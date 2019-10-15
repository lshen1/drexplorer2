
#' The Hill equation function
#' 
#' Hill equation to return y = f(Dose). HS > 0 is assumed 
#'
#' @param Dose the concentration of drug (original scale).
#' @param Einf the bottom asymptote of response (e.g. at drug concentration of Inf).
#' @param E0 the top asymptote of response (e.g. at drug concentration of 0).
#' @param logEC50 natural logarithm of EC50.
#' @param HS Hill slope.
#' @return the response at dose = Dose (between 0 and 1).
#' @export
fHill <- function(Dose, Einf, E0, logEC50, HS) {
  res <- Einf + (E0-Einf) / (1+(Dose/exp(logEC50))^HS)
  return(res)
}


#' Fit Hill equation on dose response data
#'
#' This function implements S3-style OOP. Methods such as predict, plot and lines are available.
#' Notice that control data (dose = 0) is usually not fitted (supplied) in Hill equation.
#'
#' @param d dose original dose.
#' @param y response relative viability. 
#' 
#' @return it returns a hillFit object consisting a vector of estimations for 'E0', 'Einf',
#'  'EC50', 'HS', 'Emax', 'MSE' (Mean Squared Error), 'Rsq' (R squared),
#'  'RSE' (Residual Standard Error).
#' @importFrom stats nls residuals
#' @export
#' @examples 
#' sdat <- prepDRdat(drMat = ryegrass[, c(2, 1)], alpha = 0.01, fitCtr = FALSE,
#'                   standardize = TRUE)$dat
#' fit_hill <- hillFit(d = sdat$dose, y = sdat$response)
#' fit_hill[1:length(fit_hill)]
hillFit <- function(d, y){
  if(length(d) != length(y)){
    stop('Please specify equal length in d (dose) and y (response)!\n')
  }
  dat <- data.frame(dose = d, response = y)
  # this might fail due to convergence
  # need to add boundary: in the order of start
  # EC50 must between dmin and dmax
  lb <- c(0, 0, log(min(dat$dose)), 1e-5)
  ub <- c(2, 1, log(max(dat$dose)), 1e5)
  # warnOnly=TRUE: always give an result although it might be non-convergent ---> this will always return a result but Rsq can be really poor
  fit <- try(nls(y ~ Einf + (E0-Einf) /(1+(dose/exp(logEC50))^HS),
                 start = list(E0 = 1, Einf = 0, logEC50 = log(median(dat$dose)), HS = 1), 
                 algorithm = "port", trace = FALSE, lower = lb, upper = ub,
                 control = list(maxiter = 50000, warnOnly = TRUE), data = dat), silent = TRUE)
  res <- rep(NA, 8)
  names(res) <- c('E0', 'Einf', 'EC50', 'HS', 'Emax', 'MSE', 'Rsq', 'RSE')
  if(!inherits(fit, 'try-error')){
    cc <- coef(fit)
    rss <- sum(residuals(fit)^2) # residual sum of squared
    tss <- sum((y-mean(y, na.rm = T))^2, na.rm = T)
    MSE <- rss/length(y)
    Rsq <- 1-rss/tss
    df <- sum(!is.na(y))-4 # four parameters in for Hill model
    rse <- sqrt(rss/df)# residual standard error: sqrt(RSS/df)
    res['EC50'] <- exp(cc['logEC50'])
    res['E0'] <- cc['E0']
    res['Einf'] <- cc['Einf']
    res['HS'] <- cc['HS']
    res['MSE'] <- MSE
    res['Rsq'] <- Rsq
    res['RSE'] <- rse
    res['Emax'] <- predict(fit, list(dose=max(d, na.rm = T))) # response at maximum dose using nls predict method
  }	
  attr(res, 'class') <- c('numeric', 'hillFit')
  attr(res, 'fitDat') <- dat
  attr(res, 'fit') <- fit
  return(res)
}


#' predict method for hillFit class
#' 
#' @rdname predict
#' @method predict hillFit
#' @aliases predict,hillFit
#' 
#' @param fit a hillFit/nci60Fit object as returned by hillFit/nci60Fit.
#' @param newData dose in original scale to be used for prediction of response.
#' 
#' @return predicted response at specified dose.
#' @export
predict.hillFit <- function(fit, newData = NULL){
  fitDat <- toget(fit, 'fitDat')
  if(is.null(newData)) newData <- fitDat$dose 
  dose <- newData
  y <- fHill(Dose = dose, Einf = fit['Einf'], E0 = fit['E0'],
             logEC50 = log(fit['EC50']), HS = fit['HS'])
  return(y)
}



getPlotDat_hill <- function(fit){
  fitDat <- toget(fit, 'fitDat')
  gg <- format_grid(fitDat$dose)
  top <- gg$top
  bot <- gg$bot
  xGrid <- gg$xGrid
  y <- predict(fit, newData = xGrid)
  if(length(xGrid) > 10000){
    ind1 <- which(diff(log10(xGrid)) > 1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
    ind2 <- floor(seq(max(ind1) + 1, length(xGrid), length.out = 1000))
    indSel <- c(ind1, ind2)
  } else {
    indSel <- 1:length(xGrid)
  }
  return(list(fitDat = fitDat, y = y, indSel = indSel, xGrid = xGrid))
}



#' plot method for hillFit class
#' 
#' @rdname plot
#' @method plot hillFit
#' @aliases plot,hillFit
#' 
#' @param fit a hillFit object as returned by hillFit.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param xlim y limit for display.
#' @param ylim x limit for display.
#' @param main main title.
#' @param cex.main cex for main title.
#' @param cex.axis cex for axis annotation.
#' @param pcol color for points. 
#' @param lcol color for lines.
#' @param lwd line width used to specify the appearance of fitted curve.
#' @param ... additional parameters, not implemented.
#' 
#' @export
plot.hillFit <- function(fit, xlab = "Log10(Dose)", ylab = "Relative viability",
                         xlim = NULL, ylim = NULL,
                         main = 'Fit of Hill equation', 
                         cex.main = 1, cex.axis = 1,
                         pcol = 'black', lcol = 'black', lwd = 2, ...){
  plL <- getPlotDat_hill(fit)
  fitDat <- plL$fitDat
  xGrid <- plL$xGrid
  indSel <- plL$indSel
  y <- plL$y
  if(is.null(ylim)) ylim <- range(pretty(fitDat$response))
  if(is.null(xlim)) xlim <- range(pretty(log10(fitDat$dose)))		
  with(fitDat, plot(log10(dose), response, 
                    xlab = xlab, ylab = ylab, main = main,
                    xlim = xlim, ylim = ylim, cex.main = cex.main,
                    cex.axis = cex.axis))
  lines(log10(xGrid)[indSel], y[indSel], col = lcol, lwd = lwd)
}




#' lines method for hillFit class
#' 
#' @rdname lines
#' @method lines hillFit
#' @aliases lines,hillFit
#' 
#' @param fit a hillFit object as returned by hillFit.
#' @param col line color.	   
#' @param lwd line width.  
#' @param show_points whether to add points for dose-response paires (dose not 0). 
#' @param pcol color for points; only relevant if show_points is TRUE.
#' @param pch pch for points; only relevant if show_points is TRUE.	   
#' @param ... additional parametrs passed to generic lines() function.
#' 
#' @export
lines.hillFit <- function(fit, col = 1, lwd = 2,
                          show_points = FALSE, pcol = 'black', pch = 16, ...){
  plL <- getPlotDat_hill(fit)
  fitDat <- plL$fitDat
  xGrid <- plL$xGrid
  indSel <- plL$indSel
  y <- plL$y
  lines(log10(xGrid)[indSel], y[indSel], col = col, lwd = lwd)
  if(show_points){
    with(fitDat, points(log10(dose), response, col = pcol, pch = pch))
  }
}





