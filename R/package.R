
#' Dose-response Explorer
#'
#' There is a great variety of models developed for dose-response data, many of which have been
#'  implemented in the drc and DoseFinding packages. drexplorer combines both packages to aid
#'  the user to visually examine and compare how existing models perform on the data. Another
#'  important feature for drexplorer is to allow the user to identify outlier measurements and
#'  visually examine how these outliers affect the fitted model.
#' @docType package
#' @name drexplorer-package
#' @aliases drexplorer
#' @author Li Shen, Pan Tong, and Kevin R Coombes
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit}}
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @importFrom epiR epi.ccc
#' @importFrom drc drm LL.2 LL.3 LL.3u LL.4 LL.5 W1.2 W1.3 W1.4 W2.2 W2.3 W2.4 BC.4 BC.5 LL2.2 LL2.3 LL2.3u LL2.4 LL2.5 AR.2 AR.3 
#' @importFrom DoseFinding fitMod emax sigEmax exponential quadratic linear linlog logistic 
NULL


setOldClass("drc")
setOldClass("DRMod")
setClassUnion("drFit0", list("drc", "DRMod"))

#' Class \code{"drFit"}
#' 
#' Class \code{"drFit"} defines a results from drc or DoseFinding
#'
#' @name drFit-class
#' @rdname drFit-class
#' @exportClass drFit
#' @section Creating Objects: 
#' Although objects can be created directly using \code{new("drFit", ...)},
#'  the most common usage will be to pass dose-response data to the \code{drFit} function.
#'
#' @slot fit A fit object by either drc or DoseFinding package. 
#' @slot fitDat The actual data matrix (after excluding outliers) used to fit a dose-response model.
#' @slot originalDat The original dose response data, a matrix.
#' @slot alpha A scalar for significance level. This specifies the significance level to identify outliers
#' @slot fitCtr A logical value specifying whether to include the control points into the model fitting. 
#' @slot tag A string (either 'drc' or 'DoseFinding') tracking which package is used for model fitting. 
#' @slot modelName model name of the fitted model
#' @slot nPar number of parameters for the fitted model
#' @slot info A list that holds information related to the model, i.e. Residual Standard Error (rse).
#' @exportClass drFit
setClass('drFit', representation(fit = 'drFit0',
                                 fitDat = 'matrix',
                                 originalDat = 'matrix',
                                 alpha = 'numeric',
                                 fitCtr = 'logical',
                                 tag = 'character',
                                 modelName = 'character',
                                 nPar = 'numeric',
                                 info = 'list'))

#' @name predict
#' @rdname predict-methods
#' @exportMethod predict
setGeneric('predict', package = 'stats')

#' predict method for drFit object
#'
#' @param object a model object for which prediction is desired.
#' @param ... additional arguments affecting the predictions produced.
#' @param newData a vector of dose levels to predict at. Default is to predict response at observed dose levels.
#' @return a vector of predicted values
#' @rdname predict-methods
#' @aliases predict,drFit-method
setMethod('predict', signature(object = 'drFit'),
          function(object, newData) {
            if(missing(newData)) newData <- object@originalDat[, 1] 	
            fitObj <- object@fit
            tag <- object@tag
            y <- rep(NA, length(newData))
            if(tag == "drc"){ # using model fitted from drc
              y <- predict(fitObj, newdata = data.frame(dose = newData), od = TRUE)
            }
            if(tag == "DoseFinding"){ # using model fitted from DoseFinding
              # linlog always return y=-1.37; problematic for this model!
              modelName <- attributes(fitObj)$model
              func <- get(modelName) # function name, i.e. sigEmax
              coy <- as.list(fitObj$coefs)
              coy$dose <- newData
              y <- do.call(func, coy)
            }
            return(y)
          })


#' @name plot
#' @rdname plot-methods
#' @exportMethod plot
setGeneric('plot', package = 'graphics')

#' plot method for drFit object
#' 
#' @param x a drFit object.
#' @param y y is not required in a drFit object.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @param pchs pchs (a vector) specify symbols to show different data points. In particular,
#'	pchs[1] is used to show regular points. pchs[2] is used to show outliers at 0.05
#'	significance level and pchs[3] is used to show outliers at 0.01 significance level.  
#' @param cols Similarly, cols[1] is used to color regular points, cols[2] to color outliers
#'  at 0.05 significance level and cols[3] to color outliers at 0.01 significance level. 
#' @param col color used to specify the appearance of fitted curve.
#' @param lwd line width used to specify the appearance of fitted curve.
#' @param lty line type used to specify the appearance of fitted curve.
#' @param addLegend currently deprecated. Please specify through type argument. Whether to add
#'  legend indicating outlier status.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param xlim y limit for display.
#' @param ylim x limit for display.
#' @param main main title. 
#' @param type the plot type. The user can specify more flexibly by specifying a vector of
#'  entities to plot including:
#' 	line: dose response curve
#' 	points: observed data points
#' 	control: control data points
#' 	legend: legend for outlier status
#' @param style deprecated. Please specify through type argument. when style='full', observed
#'  dose-response pairs are plotted including controls, fitted curve is superimposed and legend
#'  for outliers indicated; when style='simple', only fitted curve is plotted (without dose-response
#'  points and of course, not outlier status). 
#' @param bty bty passed to legend.
#' @param h horizontal line to add indicating e.g. IC50 (h=0.5).
#' @param cex.main cex for main title.
#' @param cex.axis cex for axis annotation.
#' @param cex.lab cex for axis label.
#' @param axes whether to add axes.
#' @param return whether return data for the plotted lines; this is useful if we want to
#'  reconstruct the lines elsewhere.
#'  
#' @importFrom graphics abline axTicks axis legend lines par points
#' @importFrom plyr ddply summarise .
#'  
#' @rdname plot-methods
#' @aliases plot,drFit-method  
setMethod('plot', signature(x = 'drFit'),
          function(x, pchs = c(16, 17, 15), cols = c(1, 2, 3), col = 4, lwd = 2,
                   lty = 1, addLegend = FALSE,
                   xlab = "Dose", ylab = "Relative viability",
                   xlim = NA, ylim = NA, main,
                   type = c('line', 'points', 'control', 'legend', 'sem'),
                   style = '', bty = 'n', h = c(0.5),
                   cex.main = 1, cex.axis = 1, cex.lab = 1, axes = TRUE, return = FALSE) {
            if(missing(main)) main <- attributes(x@fit)$model
            if(style == 'full') type <- c('line', 'points', 'control', 'legend')
            if(style == 'simple') type <- c('line')
            if(addLegend) type <- c(type, 'legend')

            resPrepDR <- x@info$resPrepDR # result of prep DR data 
            datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
            isCtrl <- resPrepDR$isCtrl # global indicator if data is control
            isTrt <- resPrepDR$isTrt # global indicator if data is treatment
            hasCtrl <- resPrepDR$hasCtrl
            indicator1 <- resPrepDR$indicator1
            indicator2 <- resPrepDR$indicator2
            dose <- datAll$dose # may include 0 dose
            dose1 <- dose[isTrt]
            trtScaled <- datAll$response[isTrt] # response at nonzero dosage
            ctrScaled <- datAll$response[isCtrl] # response at nonzero dosage
            
            if(is.na(ylim[1])) ylim <- c(0, range(pretty(datAll$response))[2])

            pCols <- rep(cols[1], nrow(datAll)) # cols[1] for regular points
            pCols[indicator1] <- cols[2] # cols[2] for outliers at 5% significance level
            pCols[indicator2] <- cols[3] # cols[3] for outliers at 1% significance level
            pPchs <- rep(pchs[1], nrow(datAll))
            pPchs[indicator1] <- pchs[2]
            pPchs[indicator2] <- pchs[3]
            # setup grid
            gg <- format_grid(dose1)
            top <- gg$top
            bot <- gg$bot
            xGrid <- gg$xGrid
            # modified on 02/23/2014: use a subset of points so that the pdf figure will not be too large
            ## too many points and leads to a large pdf: use a subset of the points
            ## modified on 20150406: grid out has 5000 points, this should be reasonable; no subsetting is needed if grid length less than 10000
            if(length(xGrid) > 10000){
              ind1 <- which(diff(log10(xGrid))>1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
              ind2 <- floor(seq(max(ind1)+1, length(xGrid), length.out = 1000))
              indSel <- c(ind1, ind2)
            } else {
              indSel <- 1:length(xGrid)
            }
            y <- predict(x, newData = xGrid)
            dd <- data.frame(line_x = log10(xGrid)[indSel], line_y = y[indSel])
            res <- list(dat = dd, ylim = ylim)
            datAll_trt <- subset(datAll, dose != 0)
            mysmry <- ddply(datAll_trt, .(dose), summarise, N = length(response),
                            meanResponse = mean(response, na.rm = TRUE),
                            stdResponse = sd(response, na.rm = TRUE))
            N_rep <- mean(mysmry$N, na.rm = T) # N replicates, averaged across doses
            mysmry$ss <- mysmry$stdResponse/sqrt(N_rep)
            # where no plot is needed
            if(return) return(res)
            # initialize the plot: the scaled values for ctr and trt
            par(bg = "white")
            ## the actual data points
            if(is.na(xlim[1])) xlim <- log10(c(bot, top))
            if(!'plot' %in% type ){
              # just plot points?
              # 	plot(log10(dose1), trtScaled, ylim=ylim, xlim=xlim, axes=axes, lty=lty, lwd=lwd,
              #       xlab=xlab, ylab=ylab, main=main, col=pCols[isTrt], pch=pPchs[isTrt], cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab)
              # NOW IN LOG annotation
              plot(log10(dose1), trtScaled, ylim = ylim, xlim = xlim, axes = axes,
                   lty = lty, lwd = lwd,
                   xlab = xlab, ylab = ylab, main = main,
                   col = pCols[isTrt], pch = pPchs[isTrt], cex.main = cex.main, cex.axis = cex.axis,
                   cex.lab = cex.lab, xaxt = "n")
              atx <- axTicks(1)
              labels <- sapply(atx,function(i) as.expression(bquote(10^ .(i))))
              axis(1,at=atx,labels=labels)
            } else {
              # no points; so just the axis
              # plot(log10(dose1), trtScaled, ylim=ylim, xlim=xlim, axes=axes, lty=lty, lwd=lwd, 
              #         xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab, col=pCols[isTrt], pch=pPchs[isTrt], type='n')
              plot(log10(dose1), trtScaled, ylim = ylim, xlim = xlim, axes = axes,
                   lty = lty, lwd = lwd,
                   xlab = xlab, ylab = ylab, main = main, cex.main = cex.main, cex.axis = cex.axis,
                   cex.lab = cex.lab, col = pCols[isTrt], pch = pPchs[isTrt],
                   type = 'n', xaxt = "n")
              atx <- axTicks(1)
              labels <- sapply(atx,function(i) as.expression(bquote(10^ .(i))))
              axis(1, at = atx, labels = labels)
            }
            if('control' %in% type) {
              # add control points if available
              if(hasCtrl)	points(rep(log10(bot), length(ctrScaled)), ctrScaled, pch=8)	
            }
            if(!is.null(h))	abline(h=h, col='grey')
            ### add curve from the fitted model
            if('line' %in% type)
              lines(log10(xGrid)[indSel], y[indSel], col = col, lwd = lwd, lty = lty)
            # add outlier legend
            if('legend' %in% type)
              legend("topright", c("okay", "5%", "1%"), col = cols, pch = pchs, bty = bty)
            if('sem' %in% type)
              with(mysmry, errbar(log10(dose), meanResponse,
                                  meanResponse + ss, meanResponse - ss,
                                  add = T, errbar.col = col, lwd = 0.7, col = col, cex = 0.6))
            ## add legend
            #if(addLegend) legend("topright", c("okay", "95%", "99%"), col=cols, pch=pchs) ## use significance level to remove confusion (modified on 09/03/2013)
            #if(style!='full')
            #	addLegend <- FALSE # simple style removes the outlier status
            # if alpha=1, of course no need for legend
            #if(x@alpha==1) addLegend <- FALSE
            #if(addLegend) legend("topright", c("okay", "5%", "1%"), col=cols, pch=pchs, bty=bty)
            #browser()
          })
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("response"))



#' @name lines
#' @rdname lines-methods
#' @exportMethod lines
setGeneric('lines', package = 'graphics')

#' lines method for drFit object
#' 
#' @param x a drFit object	   
#' @param col line color.
#' @param lwd line width.
#' @param lty line type.  
#' @param show_points whether to add points for dose-response paires (dose not 0). 
#' @param pcol color for points; only relevant if show_points is TRUE.
#' @param pch pch for points; only relevant if show_points is TRUE.
#' @param type the plot type. The user can specify more flexibly by specifying a vector of entities
#'  to plot including:
#' 	line: dose response curve.
#' 	points: observed data points.	   
#' @param ... additional parametrs passed to generic lines() function.
#' 
#' @importFrom graphics lines
#' 
#' @rdname lines-methods
#' @aliases lines,drFit-method
setMethod('lines', signature(x  ='drFit'),
          function(x, col = 5, lwd = 2, lty = 1, show_points = FALSE, pcol = 'black',
                   pch = 16, type = c('line'), ...) {
            if(show_points) type <- c(type, 'points')
            resPrepDR <- x@info$resPrepDR # result of prep DR data 
            datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
            isTrt <- resPrepDR$isTrt # global indicator if data is treatment
            dose <- datAll$dose # may include 0 dose
            dose1 <- dose[isTrt]
            trtScaled <- datAll$response[isTrt] # response at nonzero dosage
            # setup grid	
            gg <- format_grid(dose1)
            top <- gg$top
            bot <- gg$bot
            xGrid <- gg$xGrid
            y <- predict(x, newData = xGrid)
            datAll_trt <- subset(datAll, dose != 0)
            mysmry <- ddply(datAll_trt, .(dose), summarise, N = length(response),
                            meanResponse = mean(response, na.rm = TRUE),
                            stdResponse = sd(response, na.rm = TRUE))
            N_rep <- mean(mysmry$N, na.rm = T) # N replicates, averaged across doses
            mysmry$ss <- mysmry$stdResponse/sqrt(N_rep)
            # modified on 02/23/2014: use a subset of points so that the pdf figure will not be too large
            ## modified on 20150406: grid out has 5000 points, this should be reasonable; no subsetting is needed if grid length less than 10000
            if(length(xGrid) > 10000){
              ind1 <- which(diff(log10(xGrid)) > 1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
              ind2 <- floor(seq(max(ind1)+1, length(xGrid), length.out = 1000))
              indSel <- c(ind1, ind2)
            } else {
              indSel <- 1:length(xGrid)
            }
            if('line' %in% type)	
              lines(log10(xGrid)[indSel], y[indSel], col = col, lwd = lwd, lty = lty, ...)
            if('points' %in% type)
              points(log10(dose1), trtScaled, col = pcol, pch = pch)
            if('sem' %in% type)
              with(mysmry, errbar(log10(dose), meanResponse, meanResponse + ss, meanResponse - ss,
                                  add = TRUE, errbar.col = col, lwd = 0.7, col = col, cex = 0.6))
            # sometimes the added lines can be accompanied with points to show original data
            #if(show_points){
            #	points(log10(dose1), trtScaled, col=pcol, pch=pch)
            #}
            #lines(log10(xGrid), y, col=col, lwd=lwd)	  
          })	











# from Hmisc
#' @importFrom graphics segments strwidth
errbar <- function (x, y, yplus, yminus, cap = 0.015, main = NULL, sub = NULL, 
                    xlab = as.character(substitute(x)),
                    ylab = if (is.factor(x) || is.character(x)) "" else as.character(substitute(y)), 
                    add = FALSE, lty = 1, type = "p", ylim = NULL, lwd = 1, pch = 16, 
                    errbar.col = par("fg"), Type = rep(1, length(y)), ...) 
{
  if (is.null(ylim)) 
    ylim <- range(y[Type == 1], yplus[Type == 1], yminus[Type == 
                                                           1], na.rm = TRUE)
  if (is.factor(x) || is.character(x)) {
    x <- as.character(x)
    n <- length(x)
    t1 <- Type == 1
    t2 <- Type == 2
    n1 <- sum(t1)
    n2 <- sum(t2)
    omai <- par("mai")
    mai <- omai
    mai[2] <- max(strwidth(x, "inches")) + 0.25
    par(mai = mai)
    on.exit(par(mai = omai))
    plot(NA, NA, xlab = ylab, ylab = "", xlim = ylim, ylim = c(1, 
                                                               n + 1), axes = FALSE, ...)
    axis(1)
    w <- if (any(t2)) 
      n1 + (1:n2) + 1
    else numeric(0)
    axis(2, at = c(seq.int(length.out = n1), w), labels = c(x[t1], 
                                                            x[t2]), las = 1, adj = 1)
    points(y[t1], seq.int(length.out = n1), pch = pch, type = type, 
           ...)
    segments(yplus[t1], seq.int(length.out = n1), yminus[t1], 
             seq.int(length.out = n1), lwd = lwd, lty = lty, col = errbar.col)
    if (any(Type == 2)) {
      abline(h = n1 + 1, lty = 2, ...)
      offset <- mean(y[t1]) - mean(y[t2])
      if (min(yminus[t2]) < 0 & max(yplus[t2]) > 0) 
        lines(c(0, 0) + offset, c(n1 + 1, par("usr")[4]), 
              lty = 2, ...)
      points(y[t2] + offset, w, pch = pch, type = type, 
             ...)
      segments(yminus[t2] + offset, w, yplus[t2] + offset, 
               w, lwd = lwd, lty = lty, col = errbar.col)
      at <- pretty(range(y[t2], yplus[t2], yminus[t2]))
      axis(side = 3, at = at + offset, labels = format(round(at, 
                                                             6)))
    }
    return(invisible())
  }
  if (add) 
    points(x, y, pch = pch, type = type, ...)
  else plot(x, y, ylim = ylim, xlab = xlab, ylab = ylab, pch = pch, 
            type = type, ...)
  xcoord <- par()$usr[1:2]
  smidge <- cap * (xcoord[2] - xcoord[1])/2
  segments(x, yminus, x, yplus, lty = lty, lwd = lwd, col = errbar.col)
  if (par()$xlog) {
    xstart <- x * 10^(-smidge)
    xend <- x * 10^(smidge)
  }
  else {
    xstart <- x - smidge
    xend <- x + smidge
  }
  segments(xstart, yminus, xend, yminus, lwd = lwd, lty = lty, 
           col = errbar.col)
  segments(xstart, yplus, xend, yplus, lwd = lwd, lty = lty, 
           col = errbar.col)
  return(invisible())
}





