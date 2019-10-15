

iaiModels <- c('sigEmax', 'linear', 'logistic', # deduced by Dianne
               'LL.3u', 'LL.5' # drc
               ) 


decomposeMixtureDose <- function(doseMixture, ratio, islog = TRUE){
  if(islog){
    dm <- 10^doseMixture
  } else {
    dm <- doseMixture
  }
  p1 <- 1/(1+ratio)
  d1 <- dm*p1
  d2 <- dm*(1-p1)
  return(data.frame(d = dm, d1 = d1, d2 = d2))
}



#' Re-implemented IAI calculation with model selection and CI
#' 
#' @param d1 dose for drug 1. 
#' @param d2 dose for drug 2. 
#' @param e corresponding percentile in the range [0, 1]; response is 1-E.
#' @param E a vector of responses (between 0 and 1) where IAI and confidence interval are to be
#'  computed from.                                           
#' @param name1 drug 1 name.
#' @param name2 drug 2 name.
#' @param ratio ratio of d1/d2.
#' 
#' @importFrom stringr str_detect
#' @importFrom plyr mutate
#' 
#' @return a list. 
fitIAI_mod0 <- function(d1, d2, e, E = seq(0.05, 0.95, 0.005),
                        name1 = 'Drug A', name2 = 'Drug B', ratio){
  # originally fitIAI_fixedRatio_MSEL_CI()
  ind1 <- d1!=0 & d2==0
  ind2 <- d1==0 & d2!=0
  ind12 <- d1!=0 & d2!=0 # where the dose differs
  ### specify models
  models <- iaiModels
  dat1 <- data.frame(Dose = d1[ind1], Response = e[ind1])
  dat2 <- data.frame(Dose = d2[ind2], Response = e[ind2])
  # notice the combined dose is used for fitting combination data
  # this is consistent to Dianne's code
  dat12 <- data.frame(Dose = d1[ind12]+d2[ind12], Response = e[ind12]) # y~A+B
  dat12A <- data.frame(Dose = d1[ind12], Response = e[ind12])  # y~A
  dat12B <- data.frame(Dose = d2[ind12], Response = e[ind12])  # y~B
  fitdat <- function(dat, drug){
    resPrepDR <- prepDRdat(dat, alpha = 1, fitCtr = FALSE, standardize = FALSE)
    info_experimentSetup <- prepDRdat2expInfo(resPrepDR)
    resL <- fitOneExp(dat = dat, drug = drug, cellLine = '', models = models,
                      percentile = E, plot = FALSE, fitCtr = FALSE,
                      transparency = 0.95, standardize = FALSE, interpolation = FALSE, unit='')
    resL$ICx <- data.frame(resL$ICx)
    resL$info_experimentSetup <- info_experimentSetup
    resL
  }
  ### fit drugs
  resL_1 <- fitdat(dat1, name1)
  resL_2 <- fitdat(dat2, name2)
  resL_12 <- fitdat(dat12, str_c(name1, name2, sep = '_')) # this is fitted with 
  # sometimes for mixute A+B, we may plot y~A, y~B, y~A+B (resL_12). so this is added on 2016/01/05
  resL_12A <- fitdat(dat12A, str_c(name1, name2, sep = '_')) # this is fitted with 
  resL_12B <- fitdat(dat12B, str_c(name1, name2, sep = '_')) # this is fitted with 
  
  ## IAI estimation requires ic value for A, B, A+B
  ## IAI variance estimation requires p (1/ratio), ic value as well as variance for A, B, A+B
  ## IAI C.I. requires number of parameters andd II
  prep4IAI <- function(resL, percentile=seq(0.1, 0.9, by=0.1)){
    # ICx here is interpolated and thus might be different from Dianne's result
    # if we set interpolation=FALSE (allow extropolation), the result is the same
    dose_log <- dose_logNA <- as.numeric(resL$ICx[1, str_detect(colnames(resL$ICx), 'IC')])
    dose_var <- as.numeric(resL$ICx_variance[1, str_detect(colnames(resL$ICx_variance), 'IC')])
    dose_logNA[is.na(dose_var)] <- NA # this leads to NA when var is NA; thus II would be NA; otherwise, II is from truncated IC value
    dfp <- resL$nParBest
    rr <- data.frame(dose = 10^dose_logNA, variance = dose_var, FA = E,
                     Response = 1-E, dfp = dfp, model = resL$models[resL$indBest]) # to original scale
    rownames(rr) <- rr$Response
    rr
  }
  dat4IAI1 <- prep4IAI(resL_1)
  dat4IAI2 <- prep4IAI(resL_2)
  dat4IAI12 <- prep4IAI(resL_12) # for combo data, dose is total dose
  p <- 1/ratio
  # this is using IC with NA; all values would be NA if there is no response that are achieved by all 3 treatment conditions! 
  # NA propogages! but this is not very useful as IAI is dependent on dose design in this case; 
  # we have a valid argument to let II being NA!
  II <- dat4IAI12$dose*(p/dat4IAI1$dose +1/dat4IAI2$dose)/(1+p) 
  II_var <- dat4IAI12$dose^2*(p^2*dat4IAI1$variance/dat4IAI1$dose^4 +
                                dat4IAI2$variance/dat4IAI2$dose^4)/((1+p)^2) +
    II^2*dat4IAI12$variance/dat4IAI12$dose^2
  n <- length(d1)
  dfp <- dat4IAI1$dfp[1] + dat4IAI2$dfp[1] + dat4IAI12$dfp[1]
  II_lb <- exp(log(II)- qt(0.975, df = n-dfp)*sqrt(II_var)/II) # 95% CI lower bound
  II_ub <- exp(log(II)+ qt(0.975, df = n-dfp)*sqrt(II_var)/II)
  datIAI <- data.frame(IAI = II, IAIvariance = II_var, IAI_lb = II_lb, IAI_ub = II_ub)
  rownames(datIAI) <- dat4IAI12$Response
  ## extract dose corresponding to each IC values
  getDose <- function(resL){
    dose_log <- as.numeric(resL$ICx[1, str_detect(colnames(resL$ICx), 'IC')])
    data.frame(dose = 10^dose_log, FA = E, Response = 1-E) # to original scale
  }
  ## extract a subset of IC values.
  getIC <- function(resL, percentile = seq(0.1, 0.9, by = 0.1), prefix = 'drug1'){
    indSel <- match(signif(percentile, 4), signif(E, 4))
    rr0 <- resL$ICx[1, str_detect(colnames(resL$ICx), 'IC')][, indSel]
    drug <- resL$drug
    rr <- data.frame(resL$info_experimentSetup, rr0)
    colnames(rr) <- str_c(prefix, colnames(rr), sep = '_')
    rr
  }
  DoseDat1 <- getDose(resL_1)
  DoseDat2 <- getDose(resL_2)
  DoseDat12 <- getDose(resL_12) # original scale
  ## IC
  ICdat1 <- getIC(resL_1, prefix = 'drug1')
  ICdat2 <- getIC(resL_2, prefix = 'drug2')
  ICdat12 <- getIC(resL_12, prefix = 'drug12')
  p1 <- ratio/(1+ratio)
  doseDat12 <- mutate(DoseDat12, d1 = dose*p1, d2 = dose*(1-p1)) # dose in the mixture
  datIAI <- data.frame(datIAI, doseDat12)
  datIAI$D1 <- DoseDat1$dose # dose for drug 1 alone
  datIAI$D2 <- DoseDat2$dose
  res <- list(resIAI_raw = datIAI, ICdat1 = ICdat1, ICdat2 = ICdat2, ICdat12 = ICdat12,
              fit1 = resL_1, fit2 = resL_2, fit12 = resL_12,
              fit12A = resL_12A, fit12B = resL_12B, ratio = ratio)
  return(res)
}



#' Re-implemented IAI calculation with model selection and CI
#' 
#' @param d1 dose for drug 1. 
#' @param d2 dose for drug 2. 
#' @param e corresponding percentile in the range [0, 1]; response is 1-E.
#' @param E a vector of responses (between 0 and 1) where IAI and confidence interval are to be
#'  computed from.                                           
#' @param name1 drug 1 name.
#' @param name2 drug 2 name.
#' @param ratio ratio of d1/d2.
#' @param tol tolerance in declaring fixed ratio.
#' @param IAIactive default is set to 0.95.
#' 
#' @return a list
#' 
#' @importFrom stringr str_split
#' 
#' @export
#' @examples 
#' fitL_mod <- fitIAI_mod(d1 = UMSCC22B[, 1], d2 = UMSCC22B[, 2], e = UMSCC22B[, 3],
#'                        name1 = 'SCH66336', name2 = '4HPR', ratio = 1)
fitIAI_mod <- function(d1, d2, e, E = seq(0.05, 0.95, 0.005), 
                       name1 = 'Drug A', name2 = 'Drug B', ratio,
                       tol = 0.1, IAIactive = 0.95){
  # originally: fitCombo_FixedRatio_MSEL_CI
  fit_fixedRay <- try(fitIAI_mod0(d1 = d1, d2 = d2, e = e,
                                  name1 = name1, name2 = name2, ratio = ratio), silent = TRUE)
  # IAI at selected response
  Esel <- round(seq(0.05, 0.95, by = 0.05), 3)
  # succeed may be unfortunate if a few of the data points are not for fixed ratio 
  # but the figure still shows a successful fit.
  succeed <- !inherits(fit_fixedRay, 'try-error')
  if(succeed){
    ## default E range is [0.05, 0.95]; E out of this range is less reliable due to saturation
    resIAIraw <- fit_fixedRay$resIAI_raw
    resIAIraw$diff <- c(0, diff(resIAIraw$FA)) # track diff in x-axis
    resIAIraw$index <- 1:nrow(resIAIraw)
    indsel <- match(Esel, round(resIAIraw$Response, 3))
    CIsel <- resIAIraw[indsel, ]
    CIsel_active <- subset(resIAIraw, IAI <= IAIactive) 
    ratio <- fit_fixedRay$ratio
    ## add IC data for decompsed drug
    decompColname <- c(str_c('decomposed_1_IC', seq(10, 90, by = 10), sep = '_'),
                       str_c('decomposed_2_IC', seq(10, 90, by = 10), sep = '_'))
    datDecomp <- data.frame(t(rep(NA, length(decompColname))))
    colnames(datDecomp) <- decompColname
    ICdat12 <- fit_fixedRay$ICdat12
    getICvec <- function(ICdat){
      isIC <- str_detect(colnames(ICdat), 'IC')
      ICvec <- as.numeric(ICdat[1, isIC])
      names(ICvec) <- sapply(str_split(colnames(ICdat)[isIC], '_'), function(x) x[[2]]) # should be consistent with getIC in fitIAI_fixedRatio_MSEL
      ICvec
    }
    ICdat1 <- fit_fixedRay$ICdat1
    ICdat2 <- fit_fixedRay$ICdat2
    IC1 <- getICvec(ICdat1)
    IC2 <- getICvec(ICdat2)
    IC12 <- getICvec(ICdat12)
    ICdat_12 <- log10(decomposeMixtureDose(IC12, ratio = ratio))	
    datDecomp[1, 1:9] <- ICdat_12$d1
    datDecomp[1, 10:18] <- ICdat_12$d2
    if(nrow(CIsel_active) >= 2){
      zz <- CIsel_active[, 'Response'] # based on response
      # WARNING: assuming active IAI is a continuous piece
      ARRsym <- abs(sum(CIsel_active$diff))
    } else {
      ARRsym <- 0
    }
  } else {
    CIsel <- data.frame(IAI = NA, FA = 1-Esel, d1 = NA, d2 = NA, D1 = NA, D2 = NA,
                        IAI_lb = NA, IAI_ub = NA)
    ARRsym <- NA
    # ratio <- str_c(exp1$d2.d1_avail, collapse = ', ')
    ratio <- ratio
  }
  resIAI <- data.frame(drug1 = name1, drug2 = name2, drug12 = str_c(name1, name2, sep = '+'), CIsel)
  IAIvec <- CIsel$IAI
  IAIvec_lb <- CIsel$IAI_lb
  IAIvec_ub <- CIsel$IAI_ub
  names(IAIvec) <- str_c('IAI_FA', 1-Esel, sep = '')
  names(IAIvec_lb) <- str_c('IAI_FA_lb', 1-Esel, sep = '')
  names(IAIvec_ub) <- str_c('IAI_FA_ub', 1-Esel, sep = '')
  # just IAI, not SE and others
  resIAIsimple0 <- data.frame(drug1 = name1, drug2 = name2,
                              drug12 = str_c(name1, name2, sep = '+'),
                              ARRsym = ARRsym, IAIactiveCutoff = IAIactive,
                              t(IAIvec), t(IAIvec_lb), t(IAIvec_ub))
  if(succeed){
    resIAIsimple <- data.frame(ratio = ratio, resIAIsimple0,
                               fit_fixedRay$ICdat1, fit_fixedRay$ICdat2,
                               fit_fixedRay$ICdat12, datDecomp)
  } else {
    resIAIsimple <- data.frame(ratio = ratio, resIAIsimple0,
                               fit_fixedRay$ICdat1, fit_fixedRay$ICdat2,
                               fit_fixedRay$ICdat12)
  }
  return(list(fit = fit_fixedRay, resIAI = resIAI, resIAIsimple = resIAIsimple))
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("IAI"))



#' Plot for the fitted object from model selection
#' 
#' @param fitL the fitted result from fitIAI_mod().
#' @param plotFA whether to plot FA value; if FALSE, dose will be used as x-axis for
#'  IAI estimates
#' @param mixtureDose options to plot mixture dose. This can be (1) 'A+B':  mixture dose is
#'  the sum of the 2 drug doses, (2) 'A': mixture dose represented by drug 1 dose, 
#'  (3) 'B': mixture dose represented by drug 2 dose, (4) 'DualAxisA' where dual axis is used
#'  for drug 1 and drug 2 dose but mixture represented by drug 1. This is useful if the two drugs
#'  have very different doses, (5) 'DualAxisB', similar to 'DualAxisA' but mixture is represented
#'  by drug 2 dose.
#' @param ylim y axis limits. Default is set to 0 to 1.2.
#' 
#' @importFrom graphics plot.new text
#' @importFrom grid plotViewport pushViewport popViewport
#' @importFrom gridBase baseViewports
#' @importFrom ggplot2 ggplot geom_line xlim ylim aes xlab ylab ggtitle scale_y_log10 geom_hline
#'    
#' @export
plotIAI_mod <- function(fitL, plotFA = F, mixtureDose = c('DualAxisA' ), ylim = c(0, 1.2)){
  res <- fitL
  fit <- res$fit
  resIAI_raw <- fit$resIAI_raw # this is finer
  resIAI <- res$resIAI # this is shorter
  resIAIsimple <- res$resIAIsimple # this is 1 row
  if(inherits(fit, 'try-error')){
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 10),
         main = 'Fitting Interaction Index failed!')
  } else {
    # library(gridBase) # mix base and ggplot2 grid figures
    op <- par(mfrow=c(2, 2))
    plot(0, type = 'n', axes = FALSE, ann = FALSE, xlim = c(0, 1), ylim = c(0, 1))
    zz <- resIAI[1, ]
    cex <- 1.3
    xx <- 0.35
    yy <- 0.84
    dy <- 0.12
    cl <- ifelse('CellLine' %in% colnames(zz), zz[, 'CellLine'], '')
    if(cl != '') text(xx, yy-dy*1, sprintf('Cell line: %s', cl), cex = cex)
    text(xx, yy-dy*2, sprintf('Drug 1: %s', zz[, 'drug1']), cex = cex)
    text(xx, yy-dy*3, sprintf('Drug 2: %s', zz[, 'drug2']), cex = cex)
    resIAI_raw <- mutate(resIAI_raw, IAItrunc = truncByLimit(IAI, Lower = 10^(-2), Upper = 10^2))
    resIAI_raw <- mutate(resIAI_raw, IAI_ub = truncByLimit(IAI_ub, Lower = 10^(-2.5), Upper = 10^2.5))
    resIAI_raw <- mutate(resIAI_raw, IAI_lb = truncByLimit(IAI_lb, Lower = 10^(-2.5), Upper = 10^2.5))
    realD <- noNA(resIAI_raw[, 1:2])
    main <- ifelse(nrow(realD) > 0, '', 'No estimate of IAI available')
    plot.new()              ## suggested by @Josh
    vps <- baseViewports()
    pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
    vp1 <- plotViewport(c(1.8,1,0,1))
    p <- ggplot(resIAI_raw, aes(x = Response, y = IAItrunc)) + geom_line()+
      xlim(c(0, 1)) + ylim(c(10^-2.5, 10^2.5))+
      geom_line(aes(x = Response, y = IAI_ub), linetype = 2, lwd = 0.5)+
      geom_line(aes(x = Response, y = IAI_lb), linetype = 2, lwd = 0.5)+
      ylab('Combination Index') + xlab('Relative viability') +
      scale_y_log10(breaks = 10^c(-2:2), labels = c(0.01, 0.1, 1, 10, 100)) +
      ggtitle(main) +
      geom_hline(yintercept = 1, linetype = 2, color = 'blue')
    print(p, vp = vp1) ## http://stackoverflow.com/questions/14124373/combine-base-and-ggplot-graphics-in-r-figure-window
    popViewport()  # this is the line to pop out figrue; critical to reset plots
    # plot dr curves
    fit1 <- fit$fit1
    fit2 <- fit$fit2
    fit12 <- fit$fit12
    fit12A <- fit$fit12A
    fit12B <- fit$fit12B
    getDoseVec <- function(fit) as.numeric(str_split(fit$info_experimentSetup$Obs_dose, ',')[[1]])
    doseVec1 <- getDoseVec(fit1)
    doseVec2 <- getDoseVec(fit2)
    doseVec12 <- getDoseVec(fit12)
    allDose <- c(doseVec1, doseVec2, doseVec12)
    if('A+B' %in% mixtureDose){
      plotOneExp(fit12, ind2plot = 'best', type = c('line', 'points', 'plot'),
                 show = 'None', pcols = 'purple', cols = 'purple', lwd = 1,
                 xlim = log10(range(allDose, na.rm = T)), main = 'Mixture using total dose')
      plotOneExp(fit1, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'red', pcols = 'red', lwd = 1, show = 'None')
      plotOneExp(fit2, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'blue', pcols = 'blue', lwd = 1, show = 'None')
      legend('bottomleft', c(fit1$drug, fit2$drug, fit12$drug),
             lty=c(1, 1, 1), col = c('red', 'blue', 'purple'), bty = 'n', cex = 0.6)
    }
    if('A' %in% mixtureDose){
      plotOneExp(fit12A, ind2plot = 'best', type = c('line', 'points', 'plot'),
                 show = 'None', pcols = 'purple', cols = 'purple', lwd = 1,
                 xlim = log10(range(allDose, na.rm = T)),
                 main = sprintf('Mixture using %s dose', fit1$drug))
      plotOneExp(fit1, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'red', pcols = 'red', lwd = 1, show = 'None')
      plotOneExp(fit2, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'blue', pcols = 'blue', lwd = 1, show = 'None')
      legend('bottomleft', c(fit1$drug, fit2$drug, fit12$drug),
             lty=c(1, 1, 1), col = c('red', 'blue', 'purple'), bty = 'n', cex = 0.6)
    }
    if('B' %in% mixtureDose){
      plotOneExp(fit12B, ind2plot = 'best', type = c('line', 'points', 'plot'),
                 show = 'None', pcols = 'purple', cols = 'purple', lwd = 1,
                 xlim = log10(range(allDose, na.rm = T)),
                 main = sprintf('Mixture using %s dose', fit2$drug))
      plotOneExp(fit1, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'red', pcols = 'red', lwd = 1, show = 'None')
      plotOneExp(fit2, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'blue', pcols = 'blue', lwd = 1, show = 'None')
      legend('bottomleft', c(fit1$drug, fit2$drug, fit12$drug), 
             lty = c(1, 1, 1), col = c('red', 'blue', 'purple'), bty = 'n', cex = 0.6)	
    } 
    if('DualAxisA' %in% mixtureDose){
      ## dual-axis
      plotOneExp(fit12A, ind2plot = 'best', type = c('line', 'points', 'plot'),
                 ylim = ylim, show = 'None', pcols = 'purple', cols = 'purple', lwd = 1,
                 xlim = log10(range(doseVec1, na.rm = T)),
                 main = sprintf('Mixture using %s dose (Dual-axis)', fit1$drug))
      plotOneExp(fit1, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'red', pcols = 'red', lwd = 1, show = 'None')
      par(new = TRUE)
      x <- fit2$fits[[fit2$indBest]]
      resPrepDR <- x@info$resPrepDR # result of prep DR data 
      datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
      isTrt <- resPrepDR$isTrt
      dose <- datAll$dose # may include 0 dose
      dose1 <- dose[isTrt]
      trtScaled <- datAll$response[isTrt] 
      ## drexplorer curves cannot be combined at this level: plotOneExp(fit2, ind2plot='best', type=c('line', 'points', 'plot'), h=NULL, col='blue', pcols='blue', lwd=1, show='None',  xlab="", ylab="", axes=FALSE, style='simple', main='')
      ll_other <- plot(x,  xlab = "", ylab = "", axes = FALSE,
                       col = "red", style = 'simple', main = '', return = T)
      dd <- ll_other$dat
      plot(dd$line_x, dd$line_y,  xlab="", ylim = ylim, ylab = "", axes = FALSE,
           type = "l", col = 'blue', lwd = 1) # WARNING: ylim assumes the previous plot is fixed with ylim=c(0, 1.2)
      points(log10(dose1), trtScaled, col = 'blue', pch = 20)
      axis(3, xlim = range(dd$line_x), lwd = 1, line = -0.2,
           col = 'blue', col.axis = 'blue', lty = 2)
      par(xpd = T) # so that outside legend is possible
      legend('bottomleft', c(fit1$drug, fit2$drug, fit12$drug),
             lty = c(1, 1, 1), col = c('red', 'blue', 'purple'), bty = 'n', cex = 0.6)	
      par(xpd = F)	
    }
    if('DualAxisB' %in% mixtureDose){
      ## dual-axis
      plotOneExp(fit12B, ind2plot = 'best', type = c('line', 'points', 'plot'),
                 show = 'None', pcols = 'purple', cols = 'purple', lwd = 1,
                 xlim = log10(range(doseVec2, na.rm = T)),
                 main = sprintf('Mixture using %s dose (Dual-axis)', fit2$drug))
      plotOneExp(fit2, ind2plot = 'best', type = c('line', 'points'),
                 cols = 'blue', pcols = 'blue', lwd = 1, show = 'None')
      par(new = TRUE)
      x <- fit1$fits[[fit1$indBest]]
      resPrepDR <- x@info$resPrepDR # result of prep DR data 
      datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
      isTrt <- resPrepDR$isTrt
      dose <- datAll$dose # may include 0 dose
      dose1 <- dose[isTrt]
      trtScaled <- datAll$response[isTrt] 
      ## drexplorer curves cannot be combined at this level: plotOneExp(fit2, ind2plot='best', type=c('line', 'points', 'plot'), h=NULL, col='blue', pcols='blue', lwd=1, show='None',  xlab="", ylab="", axes=FALSE, style='simple', main='')
      ll_other <- plot(x,  xlab = "", ylab = "", axes = FALSE,
                       col = "red", style = 'simple', main = '', return = T)
      dd <- ll_other$dat
      plot(dd$line_x, dd$line_y,  xlab = "", ylim = ylim, ylab = "", axes = FALSE,
           type = "l", col = 'red', lwd = 1) # WARNING: ylim assumes the previous plot is fixed with ylim=c(0, 1.2)
      points(log10(dose1), trtScaled, col = 'red', pch = 20)
      axis(3, xlim = range(dd$line_x), lwd = 1,line = -0.2,
           col = 'red', col.axis = 'red', lty = 2)
      par(xpd = T) # so that outside legend is possible
      legend('bottomleft', c(fit1$drug, fit2$drug, fit12$drug),
             lty = c(1, 1, 1), col = c('red', 'blue', 'purple'), bty = 'n', cex = 0.6)	
      par(xpd = F)	
    } 
    # plot IC values
    ICdat1 <- fit$ICdat1
    ICdat2 <- fit$ICdat2
    ICdat12 <- fit$ICdat12
    ratio <- fit$ratio
    getICvec <- function(ICdat){
      isIC <- str_detect(colnames(ICdat), 'IC')
      ICvec <- as.numeric(ICdat[1, isIC])
      names(ICvec) <- sapply(str_split(colnames(ICdat)[isIC], '_'), function(x) x[[2]]) # should be consistent with getIC in fitIAI_fixedRatio_MSEL
      ICvec
    }
    IC1 <- getICvec(ICdat1)
    IC2 <- getICvec(ICdat2)
    IC12 <- getICvec(ICdat12)
    ICdat_12 <- decomposeMixtureDose(IC12, ratio=ratio)
    if(plotFA){
      plot(1, xlim = c(0, max(doseVec1, na.rm = T)+0.2),
           ylim = c(0, max(doseVec2, na.rm = T)+0.2),
           xlab = sprintf('%s dose', fit1$drug), ylab = sprintf('%s dose', fit2$drug), type = "n")
      points(c(max(doseVec1, na.rm=T), 0), c(0, max(doseVec2, na.rm=T)), pch = 15)
      lines(c(10^IC1['IC30'], 0), c(0, 10^IC2['IC30']), type='l', col = 'black')
      lines(c(10^IC1['IC50'], 0), c(0, 10^IC2['IC50']), type='l', col = 'red')
      lines(c(10^IC1['IC70'], 0), c(0, 10^IC2['IC70']), type='l', col = 'blue')
      ### mixture
      td <- ICdat_12[c('IC30', 'IC50', 'IC70'), ]
      with(ICdat_12, lines(c(0, d1), c(0, d2), type = 'l', lty = 2))
      with(td, points(d1, d2, pch = 17, col = c('black', 'red', 'blue')))
      legend('top', c('FA=0.3', 'FA=0.5', 'FA=0.7'), col = c('black', 'red', 'blue'),
             lty = c(1, 1, 1), horiz = TRUE, bty = 'n', cex = 0.7,
             text.col = c('black', 'red', 'blue'))
    }
    par(op)
  }
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("IAI_lb", "IAI_ub", "IAItrunc", "Response"))







