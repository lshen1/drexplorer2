
#' Implementation of NCI60 method (GI50, TGI, LC50) for dose response data
#'
#' This function implements S3-style OOP. Methods such as predict, plot and lines are available.
#'  Notice that here the response is assumed to be calculated as in:
#'  http://dtp.nci.nih.gov/branches/btb/ivclsp.html. The response has range -1 to 1. 
#' 
#' The original NCI60 method applies linear interpolation to estimate GI50, TGI and LC50.
#'  Tong, F. P. improved this by using 4-parameter logistic model, which is the LL.4 model in
#'  drexplorer (from drc package). For resistant cell lines, even the LL.4 cannot be fitted since 
#'  the curve is essentially flat. Here we fit multiple models throught the model selection framework
#'  in drexplorer and use the best model to estimate GI50, TGI and LC50. This will almost give an
#'  estimate for all extreme data without giving NA.
#'
#' @param d dose original dose.
#' @param y response percentage of growth inhibition as defined by NCI60 method.
#' @param interpolation whether to use interpolation to estimat GI50, TGI and LC50.
#' @param log.d whether to return log10 transformed value for GI50, TGI and LC50
#'  (the dose concentration).
#' @return a nci60Fit object consisting different attributes. It includes a vector with GI50,
#'  TGI and LC50. 
#' @references Shoemaker, R. H. The NCI60 Human Tumour Cell line Anticancer Drug Screen.
#'  Nature Reviews, 6: 813-823, 2006.
#' @references Tong, F. P. (2010). Statistical Methods for Dose-Response Assays.
#' 
#' @export
#' @examples 
#' fit_nci60 <- nci60Fit(d = datNCI60$Dose, y = datNCI60$Growth/100)
#' fit_nci60[1:length(fit_nci60)]
nci60Fit <- function(d, y, interpolation = TRUE, log.d = FALSE){
  
  dat <- data.frame(dose = d, response = y)
  fitL <- fitOneExp(dat = dat, standardize = F, drug = '', cellLine = '', unit = '',
                    models = drModels(return = TRUE, verbose = FALSE)[['recommendedModels']],
                    alpha = 1, interpolation = TRUE)
  fitBest <- fitL$fits[[fitL$indBest]]
  fit <- fitBest
  # dose = 0 at precision 1e15 (not enough for 1e-5) would be accepted as control
  indtrt <- which(round(d, 15) != 0)
  dmin <- min(d[indtrt], na.rm = TRUE)
  dmax <- max(d[indtrt], na.rm = TRUE)
  if(log.d == FALSE){
    dmin <- log10(dmin)
    dmax <- log10(dmax)
  }
  
  if(!inherits(fit, 'try-error')){
    res <- findDoseGivenResponse(fit, response = c(0.5, 0, -0.5),
                                 interpolation = interpolation, log.d = log.d)
    if(log.d == FALSE){
      TGI <- log10(res[2])
    } else {
      TGI <- res[2]
    }
    auc_gr <- computeAUC(fit, dmin = dmin, dmax = TGI, islogd = TRUE)[1]
    auc_gr0 <- computeAUC(fit, dmin = dmin, dmax = dmax, islogd = TRUE)[2]
    aoc_gr <- 1-(auc_gr/auc_gr0)
    auc_ld <- computeAUC(fit, dmin = TGI, dmax = dmax, islogd = TRUE)[1]
    aoc_ld <- abs(auc_ld)/auc_gr0
    auc_gri <- (auc_gr + (auc_gr0 - abs(auc_ld)))/(2*auc_gr0)
    res <- c(res, aoc_gr, aoc_ld, auc_gri)
  } else {
    res <- rep(NA, 6)
  }
  base <- ifelse(log.d, 10, 1)
  names(res) <- c('GI50', 'TGI', 'LC50', "AOC_GR", "AOC_LD", "AUC_GRI")
  attr(res, 'class') <- c('numeric', 'nci60Fit')
  attr(res, 'fitDat') <- dat
  attr(res, 'fit') <- fit
  attr(res, 'base') <- base
  return(res)
}



#' predict method for nci60Fit class
#' 
#' @rdname predict
#' @method predict nci60Fit
#' @aliases predict,nci60Fit
#'
#' @export
predict.nci60Fit <- function(fit, newData = NULL){
  fitDat <- toget(fit, 'fitDat')
  fitDR <- toget(fit, 'fit') #drFit object
  if(is.null(newData)) newData <- fitDat$dose
  if(!inherits(fitDR, 'try-error')){ 
    y <- predict(fitDR, newData=newData)
  } else {
    y <- rep(NA, length(newData))
  }
  return(y)
}


getPlotDat_nci60 <- function(fit){
  fitDat <- toget(fit, 'fitDat')
  gg <- format_grid(fitDat$dose)
  top <- gg$top
  bot <- gg$bot
  xGrid <- gg$xGrid
  y <- predict(fit, newData = xGrid)
  ## too many points and leads to a large pdf: use a subset of the points
  if(length(xGrid) > 10000){
    ind1 <- which(diff(log10(xGrid)) > 1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
    ind2 <- floor(seq(max(ind1) + 1, length(xGrid), length.out = 1000))
    indSel <- c(ind1, ind2)
  } else {
    indSel <- 1:length(xGrid)
  }
  return(list(fitDat = fitDat, y = y, indSel = indSel, xGrid = xGrid))
}



#' plot method for nci60Fit class
#' 
#' @rdname plot
#' @method plot nci60Fit
#' @aliases plot,nci60Fit
#' 
#' @param cex.lab cex for label.
#' @param axes axes.
#' @param pch point size.
#' @param lty line type.
#' @param ltyh horizontal line type.
#' @param type default is set to lines and points.
#' @param h horizontal line to add indicating e.g. GI (h=0.5), TGI (h=0), LD50 (h=-0.5)
#' 
#' @export
plot.nci60Fit <- function(fit, xlab = "Dose", ylab = "Relative growth",
                          main = 'Fit of NCI60 method',
                          xlim = NULL, ylim = NULL,
                          cex.main = 1, cex.axis = 1, cex.lab = 1, axes = TRUE,
                          pcol = 'black', pch = 1, lcol = 'black', lty = 1, ltyh = 1, lwd = 2,
                          type = c('line', 'points'), h = c(-0.5, 0, 0.5), ...){

  plL <- getPlotDat_nci60(fit)
  fitDat <- plL$fitDat
  xGrid <- plL$xGrid
  indSel <- plL$indSel
  y <- plL$y
  if(is.null(ylim)) ylim <- c(-1, 1)
  if(any(fitDat$response>ylim[2])) ylim[2] <- max(fitDat$response)
  if(is.null(xlim)) xlim <- range(pretty(log10(fitDat$dose)))		
  if('points' %in% type){
    # just plot points
    # NOW IN LOG annotation
    plot(log10(fitDat$dose), fitDat$response, ylim = ylim, xlim = xlim,
         axes = axes, lty = lty, lwd = lwd,
         xlab = xlab, ylab = ylab, main = main, col = pcol, pch = pch,
         cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, xaxt = "n")
    atx <- axTicks(1)
    labels <- sapply(atx,function(i) as.expression(bquote(10^ .(i))))
    axis(1, at = atx, labels = labels)
    abline(h = h, col = 'grey', lty = ltyh)
  } else {
    # no points; so just the axis, blank plot
    plot(log10(fitDat$dose), fitDat$response, ylim = ylim, xlim = xlim,
         axes = axes, lty = lty, lwd = lwd,
         xlab = xlab, ylab = ylab, main = main, col=pcol, pch=pch,
         cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, type = 'n', xaxt = "n")
    atx <- axTicks(1)
    labels <- sapply(atx,function(i) as.expression(bquote(10^ .(i))))
    axis(1, at = atx, labels = labels)
    abline(h = h, col = 'grey', lty = ltyh)
  }
  if('line' %in% type)
    lines(log10(xGrid)[indSel], y[indSel], col = lcol, lwd = lwd)
} 



#' lines method for nci60Fit class
#' 
#' @rdname lines
#' @method lines nci60Fit
#' @aliases lines,nci60Fit
#' 
#' @param lcol line color.
#' @param lty line type.
#' @param type default is set to line.
#' @export
lines.nci60Fit <- function(fit, 
                           pcol = 'black', lcol = 'black',
                           lwd = 2, lty = 1, pch = 16, type = c('line')){
  plL <- getPlotDat_nci60(fit)
  fitDat <- plL$fitDat
  xGrid <- plL$xGrid
  indSel <- plL$indSel
  y <- plL$y
  if('line' %in% type)	
    lines(log10(xGrid)[indSel], y[indSel], col = lcol, lwd = lwd, lty = lty)
  if('points' %in% type)
    with(fitDat, points(log10(dose), response, col = pcol, pch = pch))
} 







