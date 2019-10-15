

model.drc <- c('LL.2', 'LL.3', 'LL.3u', 'LL.4', 'LL.5',
               'W1.2', 'W1.3', 'W1.4', 'W2.2', 'W2.3', 'W2.4',
               'AR.2', 'AR.3', 'BC.4', 'BC.5', 
               'LL2.2', 'LL2.3', 'LL2.3u', 'LL2.4', 'LL2.5')

model.DoseFinding <- c('emax', 'sigEmax', 'exponential', 'quadratic',
                       'linear', 'linlog', 'logistic')

getPackageName <- function(modelName){
  if(modelName %in% model.drc) {
    return('drc' )
  } else if(modelName %in% model.DoseFinding) {
    return('DoseFinding')
  } else {
    stop('The specified model is not supported. Check supported models by typing drModels()!\n')
  }	
}


#' Frequently used models in drexplorer2
#'
#' This function returns frequently used dose-response models for anti-cancer drug screening.
#'
#' @return a vector of recommended models
#'
#' @export
recommendedModels <- function(){
  # LL.3 might be better than linlog; linlog removed since it is not appropriate,
  # especially predict always -1.37 due to the default off parameter in the original
  # DoseFinding package
  return(c('sigEmax', 'LL.4', 'LL.5', 'linear', 'logistic'))
}

#' Show available dose-response models with direct support. 
#' 
#' This function shows available models to be passed to drFit function as modelName.
#'
#' @param return whether to return the models in a list
#' @param verbose whether to print out the models
#' @return Following is the print out:
#'   Models implemented in drc package:
#'   LL.2
#'   LL.3
#'   LL.3u
#'   LL.4
#'   LL.5
#'   W1.2
#'   W1.3
#'   W1.4
#'   W2.2
#'   W2.3
#'   W2.4
#'   BC.4
#'   BC.5
#'   LL2.2
#'   LL2.3
#'   LL2.3u
#'   LL2.4
#'   LL2.5
#'   AR.2
#'   AR.3
#'   MM.2
#'   MM.3
#'   Models implemented in DoseFinding package:
#'   emax
#'   emaxGrad
#'   sigEmax
#'   sigEmaxGrad
#'   exponential
#'   exponentialGrad
#'   quadratic
#'   quadraticGrad
#'   betaMod
#'   betaModGrad
#'   linear
#'   linearGrad
#'   linlog
#'   linlogGrad
#'   logistic
#'   logisticGrad
#'   linInt
#'   linIntGrad
#'   More specific information can be found in the original packages.
#'
#' @examples
#' drModels()
#' @export
drModels <- function(return = FALSE, verbose = TRUE){
  if(verbose) {
    cat("Models implemented in drc package:\n")
    cat(model.drc, sep = '\n')
    cat("Models implemented in DoseFinding package:\n")
    cat(model.DoseFinding, sep  ='\n')
    cat('More specific information can be found in the original packages.\n')
    cat("Models frequently used in drug screening experiments:\n")
    cat(recommendedModels(), sep = '\n')
  }
  if(return){
    return(list(drc = model.drc,
                doseFinding = model.DoseFinding,
                recommendedModels = recommendedModels()))
  }
}


#' Identifying outliers with the range to standard deviation ratio statistic.
#' 
#' This function implements the method described by D. Newman, Biometrika, 1939 to identify outliers.
#' 
#' Given measurements from controls (no drug treated), we can compute the sample
#'  standard deviation (s). The range of responses from treated samples (w) can be computed
#'  for a given dose level. Assuming the controls to have the same variation as the drug treated case, 
#'  the distribution of ratio statistic q=w/s can be derived and used to calculate if there is
#'  outliers in the treated responses as described by D. Newman, Biometrika, 1939.  
#' 
#' Note that this function works for a single dose level. When multiple dose levels exist,
#'  one need to repeatedly call this function to identify outliers at each dose level or use the
#'  flagOutliers() function which is just a wrapper.
#'
#' @param ref a vector giving the reference values to estimate sigma.
#' @param obs a vector giving the observed values where potential outlier might exist.
#' @param alpha significance level. Only 0.01, 0.05 or 1 is allowed.
#' @param recursive whether to recursively identify outliers. Currently not implemented.
#' @return indicator a logical vector specifying if the corresponding point is flagged as outlier.
#' @references Newman, D. (1939). The distribution of range in samples from a normal population, 
#'  expressed in terms of an independent estimate of standard deviation. Biometrika, 31(1/2),
#'  20-30.\url{http://www.jstor.org/stable/2334973}
#' @importFrom stats median sd  
#' @export  
#' @examples
#' set.seed(1)
#' x <- rnorm(10, 0, 1)
#' y <- c(rnorm(5, 0, 1), rnorm(1, 0, 1) + 4)
#' # the last observation in y is an outlier
#' NewmanTest(x, y, alpha = 0.05)
NewmanTest <- function(ref, obs, alpha = 0.01, recursive = FALSE){
  check.alpha(alpha)
  indicator <- rep(FALSE, length(obs))
  if(alpha[1]==1) 
    return(indicator) # at significance level 1, no outlier at all
  f <- sum(!is.na(ref))-1 # f defined in the paper 
  n <- sum(!is.na(obs)) # n defined in the paper
  if(!(f %in% c(5:20, 24, 30, 40, 60)) | !(n %in% c(3:12, 20))) {
    warning("The scenario (sample size combination) is not available for lookup!\n Treat all measurements as non-outlier")
    return(indicator)
  }
  s <- sd(ref, na.rm = TRUE) # the standard deviation estimated from control group. It is assumed to be the same as the treatment group 
  w <- diff(range(obs, na.rm = TRUE)) # the range statistic from the treatment group
  q <- w/s # the statistic used in D. Newman, Biometrika 1939
  tab <- onePercentTab # default lookup table, table IV in the paper
  if(alpha == 0.05) tab <- fivePercentTab
  cutoff <- tab[as.character(f), as.character(n)]
  if(s != 0 & q > cutoff) {
    temp <- obs - median(obs, na.rm = TRUE)
    if(abs(min(temp, na.rm = TRUE)) > abs(max(temp, na.rm = TRUE)))
      ind <- which.min(obs) # minimum is the extreme and excluded. note that NA does not affect the index due to the design of which.min()
    else
      ind <- which.max(obs)
    indicator[ind] <- TRUE
  }
  return(indicator)
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("onePercentTab", "fivePercentTab"))


#' Identifying outliers for dose-response data
#'
#' Identifying outliers for dose-response data
#'
#' @param drMat dose-response matrix. The first column being dosage and the second column
#'  being response. controls are included by specifying dose = 0.
#' @param alpha a scalar for significance level. Only 0.01, 0.05 and 1 are allowed. 
#' 	alpha = 1 is included for comparability issue. In this case, no outliers will be identified. 
#' @return indicator a vector with length nrow(drMat) specifying if a data point is outlier. 
#'  Control points always have status FALSE
#' @export
#' @examples
#' data(ryegrass, package = 'drc') # use the ryegrass data from drc package
#' drOutlier(drMat = ryegrass[, c(2, 1)], alpha = 0.05)
drOutlier <- function(drMat, alpha = 0.01) {
  check.drMat(drMat)
  dose <- drMat[, 1]
  response <- drMat[, 2]
  ctr <- response[dose == 0] # ctr measurements have dose=0
  indicator <- rep(FALSE, length(dose))
  for(d in unique(dose[dose != 0])){
    sel <- which(dose == d)
    trt <- response[sel]
    indicator[sel] <- NewmanTest(ref = ctr, obs = trt, alpha = alpha)
  }
  return(indicator)
}


#' prepare dose-reponse data
#'
#' This function scaled the response by mean reponse in control when necessary
#'
#' The standardization means response (e.g. count) is to be scaled by mean control response so that
#'  the standardized response is relative viability, usually between 0 and 1. This function first
#'  detects outlier data points (using supplied data, either scaled or not scaled). It then scale
#'  the data and split the data into two data frames: one for dose not 0 and one for dose at 0. 
#'  Ideally, outlier detection is better at original count (or signal intensity) due to asymptotic
#'  assumption, but for data already scaled, detection outlier at this level is what we can do.
#'
#' @param drMat dose-response matrix as for drFit. The first column being dosage and the second
#'  column being response. controls are included by specifying dose = 0
#' @param alpha a scalar for significance level. This specifies the significance level to identify
#'  outliers which will be excluded from model fitting. To include all data, set alpha = 1. 
#' @param fitCtr A logic vector specifying whether to include the control points into the model fitting.
#' @param standardize whether to standardize (scale) the data based on control points.
#'  This should be disabled when no control data is supplied.
#' @return a list
#' @export
prepDRdat <- function(drMat, alpha = 0.01, fitCtr = FALSE, standardize = TRUE){
  datInput <- drMat # input dat
  drMat <- noNA(drMat)
  hasCtrl <- any(drMat[, 1] == 0) # logical indicating if there is control data by testing if any dose is 0
  dose <- drMat[, 1]
  response <- drMat[, 2]
  indCtrl <- dose == 0
  indTrt <- !indCtrl
  ctr <- response[indCtrl] # ctr measurements have dose=0
  trt <- response[indTrt] 		
  if(standardize){
    # scale the control and treatment effects with control mean
    ## datAll is scalled data including control and outlier
    datAll <- data.frame(dose = dose, response = (response/mean(ctr)))
  } else {
    # no standardization
    datAll <- data.frame(dose = dose, response = response)
  }
  #### use scaled data if asked or the input data: 20150406
  indicator <- drOutlier(drMat = datAll, alpha = alpha) # for outlier removal
  indicator1 <- suppressWarnings(drOutlier(drMat = datAll, alpha = 0.05)) # for plot purpose indicating different levels of outlier status
  indicator2 <- suppressWarnings(drOutlier(drMat=datAll, alpha = 0.01)) # for plot purpose indicating different levels of outlier status
  # dat is the data used for fitting models
  if(fitCtr) { 
    dat <- datAll[!indicator, ] # just remove outlier and include control
  } else {
    tempSel <- !indicator & indTrt
    dat <- datAll[tempSel, ]
  }
  isCtrl <- indCtrl # this include outliers; 
  isTrt <- indTrt # this include outliers; 
  return(list(dat = dat, datAll = datAll, isCtrl = isCtrl, isTrt = isTrt, isOutlier = indicator,
              indicator1 = indicator1, indicator2 = indicator2, standardize = standardize,
              fitCtr = fitCtr, alpha = alpha, hasCtrl = hasCtrl, datInput = datInput))
}



#' This is a wrapper to fit dose response models from drc and DoseFinding packages. 
#' 
#' This function can fit various dose-response models by specifying a model name and package source
#'  (either drc or DoseFinding).
#' 
#' When fit the model, the response is internally scaled with mean response in control. Therefore,
#'  the fitted response is a ratio. In visualization, i.e. dose-response curve plotting, the dose
#'  is log10 transformed so that we can see the sigmoid curve. However, in computing IC50 and doing
#'  prediction, the dose is in original scale since the model is trained in the original scale of
#'  dose. Correspondingly, predicted response is also scaled response. When the user requires a
#'  log10 transformed dose, it is done after estimating dose in original scale through the
#'  computeIC() function.
#'
#' @param drMat dose-response matrix. the first column being dosage and second column being response.
#'  controls are included by specifying dose = 0
#' @param modelName a dose-response model. For available models to be specified, see drModels().
#' @param alpha a scalar for significance level. This specifies the significance level to identify
#'  outliers which will be excluded from model fitting. To include all data, set alpha = 1. 
#' @param fitCtr A logic vector specifying whether to include the control points into the model
#'  fitting.
#' @param standardize whether to standardize (scale) the data based on control points. This should
#'  be disabled when no control data is supplied
#' @return the function returns a drFit S4 object.
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit-class}}
#' @author Li Shen (\email{lishen@@mdanderson.org}),
#'  Kevin R Coombes (\email{kcoombes@@mdanderson.org}) and
#'  Pan Tong (\email{nickytong@@gmail.com})
#' @importFrom methods new
#' @importFrom stats coef
#' @export
#' @examples
#' data(ryegrass, package = 'drc') # use the ryegrass data from drc package
#' # fit a sigmaEmax model without outlier removal. the controls are excluded from model fitting
#' fit.LL.3 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3", alpha=0.01, fitCtr=FALSE)
#' fit.LL.3u <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3u", alpha=0.01, fitCtr=FALSE)
#' fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.01, fitCtr=FALSE)
drFit <- function(drMat, modelName = "sigEmax",
                  alpha = 0.01, fitCtr = FALSE, standardize = TRUE){
  package <- getPackageName(modelName) # which package source is the model implemented
  resPrepDR <- prepDRdat(drMat, alpha = alpha, fitCtr = fitCtr, standardize = standardize)
  dat <- resPrepDR$dat
  rse <- nPar <- NA
  if(package == 'drc'){
    f <- get(modelName) ## i.e. the LL.3 or LL.3u function from drc package
    model <- try(drc::drm(response ~ dose, data = dat, fct = f())) ## pass f to drm in the drc package to fit dose-response curve
    if(inherits(model, "try-error")) {
      stop(cat("Model", modelName, "failed\n"))
    } else {
      rse <- summary(model)$rseMat[1]
      nPar <- length(model$fit$par)
    }
  } else {
    model <- try(DoseFinding::fitMod('dose', 'response', data = dat, model = modelName))
    if(inherits(model, "try-error")) {
      stop(cat("Model", modelName, "failed\n"))
    } else {
      coef.est <- coef(model)
      nPar <- length(coef.est)
      rse <- sqrt(model$RSS/model$df) # sqrt(RSS/df)
    }
  }
  res <- new('drFit', fit = model, fitDat = data.matrix(dat), originalDat = data.matrix(drMat),
             alpha = alpha, fitCtr = fitCtr, tag = package, modelName = modelName, nPar = nPar,
             info = list(RSE = rse, standardize = standardize, resPrepDR = resPrepDR))
  return(res)
}



#' Calculate AUC from a fitted object
#'
#' AUC is calculated through the integrate() function based on dose-response curve.
#'
#' @param fit usually a drFit object. However, any object with a predict method would work.
#' @param dmin minimum dose. The integral range is [dmin, dmax].
#' @param dmax maximum dose. The integral range is [dmin, dmax].
#' @param islogd whether the supplied dose dmin/dmax is in log10 scale. The user should be
#'  responsible for the consistency between the actual value of dmin, dmax and islogd. If
#'  log10 transformed dmin/dmax is supplied with islogd = TRUE, the AUC is calculated based on
#'  dose-response curve with x-axis being log10(dose); On the other hand, if dmin, dmax is original
#'  scale and islogd = FALSE, the AUC is calculated based on a dose-response curve with x-axis
#'  being dose.
#' @return a vector of AUC, AUC0 and AUCs. AUC is the area under dose response curve; AUC0 is
#'  the area under the line response = 1; AUCS is AUC/AUC0
#' @importFrom stats integrate  
#' @export
#' @examples
#' data(ryegrass, package = 'drc') # use the ryegrass data from drc package
#' fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha = 0.01, fitCtr = FALSE)
#' computeAUC(fit.sigEmax, dmin = -0.027, dmax = 1.477, islogd = TRUE)
computeAUC <- function(fit, dmin, dmax, islogd = TRUE) {
  # assumption: the predict function is trained on original dose
  # the rationale: predict response at each dose; apply integration;
  # AUC can be calculated based on log10 dr curve or original dr curve (in this case, absolute AUC is smaller, scaled AUC is also smaller usually); Notice that this only changes the integration range/scale (x) without affecting response (y)
  # However, the user should supply dmin/dmax in log10 scale if islogd==TRUE (for log10 AUC) and dmin/dmax in original scale for islogd=FALSE (for untransformed AUC)
  # according to the Nature paper, AUC is defined on the log10 scale; this is also more intuitive since it matches dr curve which is in log10 scale
  # islogd: this makes it possible to calculate AUC either on log10dose or original dose. However, the user
  #      should make sure dmin and dmax is on the same scale (take log10 correspondingly)
  f_response <- function(Dose, islogd = islogd) {
    if(islogd) {
      dd <- 10^(Dose)
    } else {
      dd <- Dose
    }
    ## currently only implements for drFit
    #if(class(fit)=='drFit'){
    #	res <- predict(fit, newData=dd)
    #}
    res <- predict(fit, newData = dd)
    res  
  }
  # reference response, e.g. the line response=1; this can be used to scale AUC
  f_responseRef <- function(Dose, reference = 1, islogd = islogd) {
    if(islogd) {
      dd <- 10^(Dose)
    } else {
      dd <- Dose
    }
    res <- rep(reference, length(dd))
    res 
  }
  
  AUC <- try(integrate(f_response, lower = dmin, upper = dmax, islogd = islogd)$value, silent = TRUE)  
  AUC0 <- try(integrate(f_responseRef, lower = dmin, upper = dmax, islogd = islogd)$value, silent = TRUE)  
  if(inherits(AUC, 'try-error')) {
    cat(AUC)
    AUC = NA
  }
  if(inherits(AUC0, 'try-error')) {
    cat(AUC0)
    AUC0 = NA
  }
  res <- c(AUC = AUC, AUC0 = AUC0, AUCs = AUC/AUC0)
  return(res)
}



#' @importFrom stats uniroot
RootFindingIC <- function(drFit, percent = 0.5, log.d = TRUE,
                          lower, upper, dmin = 0, dmax = Inf, ...) {
  #### the predicted response is bounded with theoretical lower and upper bound. When the required response specified by IC is out of the theoretical range, we need to specify the returned dose
  #### originally we specify NA, in which case we may get IC50=NA, which can be due to too sensitive or too negative. Thus, it is better to give a value.
  #### we assign dmin for the scaled response larger than theoretical response; we assign dmax for response smaller than theretical minimum response.
  #### notice this only affects the IC by prediction result
  #dmin <- 1e-30 # rediculously low dose
  #dmax <- 1e30 # rediculously high dose
  #dmin <- -Inf # rediculously low dose
  #dmax <- Inf
  tag <- drFit@tag
  res <- NA
  dose <- drFit@originalDat[, 1]
  if(missing(lower)) lower <- 1e-10 # cannot allow log(lower)=log(0) as uniroot lower bound
  if(missing(upper)) upper <- 1e10 # assume maximum dose will not exceed 1e10
  myIC <- percent #
  # what to return for NA
  getNAres <- function(percent){
    res <- NA
    names(res) <- paste('response_', percent, sep = '')
    res
  }
  if(tag == "drc"){
    objFct <- drFit@fit$fct
    pm <- drFit@fit$parmMat # par mat
    objFct$"fct"(dose, t(pm))
    fl.drc <- function(ld, parm, IC) objFct$"fct"(exp(ld), parm) - IC
    ymin <- fl.drc(ld = log(upper), parm = t(pm), IC = 0)
    ymax <- fl.drc(ld = log(lower), parm = t(pm), IC = 0)
    # ymin can be NaN when the dr function is not well behaved
    # we just return NA
    if(is.nan(ymin) | is.nan(ymax)) {
      return(getNAres(percent))
    }
    # in case the curse is increasing, ymax < ymin
    # we swap them
    if(ymax < ymin) {
      tt <- ymax
      ymax <- ymin
      ymin <- tt
    }
    if(myIC > ymax) {
      res <- dmin # required response > theoretical maximum response, throw NA
    } else if(myIC<ymin){
      res <- dmax # required response < theoretical minimum response, throw NA
    } else {	
      #res <- uniroot(f.drc, parm=t(pm), IC=myIC, lower=lower, upper=upper, ...)$root # to feed with biological IC
      res <- exp(uniroot(fl.drc, parm = t(pm), IC = myIC,
                         lower = log(lower), upper = log(upper), ...)$root) 
    }
  }
  if(tag == "DoseFinding"){
    modelName <- attributes(drFit@fit)$model
    fname <- get(modelName) # function name, i.e. sigEmax
    coy <- as.list(drFit@fit$coefs) # coefs, a list for do.call
    # define root finding function: a function of d. coef specifies 4 essential parameter: e0, emax, ed50, h
    fl.DoseFinding <- function(ld, coef, IC) {
      coef$dose <- exp(ld)
      do.call(fname, coef)-IC
    }
    ymax <- fl.DoseFinding(ld = log(lower), coef = coy, IC = 0)
    ymin <- fl.DoseFinding(ld = log(upper), coef = coy, IC = 0)
    # ymin can be NaN when the dr function is not well behaved
    ## this is due to a bug in sigEmax function: e0 + eMax * dose^h/(ed50^h + dose^h); when dose=1e60, d^h=Inf; Inf/Inf=NaN; we can manually correct for this, but it is easier to set smaller upper, i.e. 1e10
    # we just return NA
    if(is.nan(ymin) | is.nan(ymax)) {
      return(getNAres(percent))
    }
    # in case the curse is increasing, ymax < ymin
    # we swap them
    if(ymax < ymin) {
      tt <- ymax
      ymax <- ymin
      ymin <- tt
    }
    if(myIC > ymax) {
      res <- dmin # required response > theoretical maximum response, throw NA
    } else if(myIC < ymin){
      res <- dmax # required response < theoretical minimum response, throw NA
    } else {	
      # now myIC is between ymax and ymin
      res <- exp(uniroot(fl.DoseFinding, coef = coy, IC = myIC,
                         lower = log(lower), upper=log(upper), ...)$root)
    }	
  }
  if(log.d) res <- log10(res)
  names(res) <- paste('response_', percent, sep = '')
  return(res)
}


icByInterpolation <- function(perc, xv, yv, ubound = max(xv), lbound = min(xv)){
  av <- abs(yv - perc)
  wv <- which(av == min(av))[1]
  yr <- range(yv, na.rm = TRUE)
  if(perc > yr[2]){ # response too large
    ic <- lbound 
  } else if(perc < yr[1]) { # response too small
    ic <- ubound
  } else { # response in range
    ic <- xv[wv]
  }
  # truncation
  if (ic > ubound) ic <- ubound ## saturated case
  if (ic < lbound) ic <- lbound
  return(ic)
}


#' @importFrom stats approx
format_grid <- function(dose1, n = 1000, resolution = 50, stepLen = NA){
  top <- 10^ceiling(log10(max(dose1))) ######### modified: sometimes only 6 doses are observed. This modification guarantees all 7 doses are present
  bot <- 10^floor(log10(min(dose1)))
  dif <- diff(sort(dose1, decreasing=FALSE)); 
  min_d <- min(dif[dif!=0]) # add sort so that dif will always be positive
  xGrid <- 10^approx(log10(sort(unique(dose1))), n = n)$y # 5000 points equal space between each dose
  return(list(top = top, bot = bot, xGrid = xGrid, stepLen = stepLen))
}


findDoseGivenResponse <- function(drFit, response = 0.50, log.d = TRUE,
                                  interpolation = TRUE, stepLen = NA, lower, upper, ...) {
  #percent <- 1-percent ### convert to biological percent; updated on 2014/02/13---> not ok: the order is reversed!
  if(interpolation == FALSE){
    if(length(response) == 1) return(RootFindingIC(drFit, response, log.d, lower = lower, upper = upper, ...))
    if(length(response) > 1) {
      res <- rep(NA, length(response))
      for(i in 1:length(response)) {
        tm <- try(RootFindingIC(drFit, response[i], log.d, lower = lower, upper = upper, ...),
                  silent = TRUE)
        if(class(tm) != 'try-error') res[i] <- tm
      }
      names(res) <- paste('response_', response, sep = '')
    }
  } else {
    resPrepDR <- drFit@info$resPrepDR # result of prep DR data 
    datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
    isTrt <- resPrepDR$isTrt # global indicator if data is treatment
    dose <- datAll$dose # may include 0 dose
    dose1 <- dose[isTrt]
    #
    gg <- format_grid(dose1=dose1, stepLen=stepLen, resolution=100)
    top <- gg$top
    bot <- gg$bot
    xGrid <- gg$xGrid
    yv <- predict(drFit, newData = xGrid) ## predicted values. dose at the original scale
    if(length(response) == 1) {
      res <- icByInterpolation(response, xv = xGrid, yv, max(dose1), min(dose1))
    } else {
      res <- rep(NA, length(response))
      for(i in 1:length(response)) {
        tm <- try(icByInterpolation(response[i], xv = xGrid, yv, max(dose1), min(dose1)),
                  silent=TRUE)
        if(class(tm) != 'try-error') res[i] <- tm
      }
      names(res) <- paste('response_', response, sep = '')
    }
    # final presentation
    if(log.d) {
      res <- log10(res)
    }
  }
  return(res)
}


#' This function uses rootfinding with fitted curve to compute IC values.
#'
#' @param drFit A drFit object as returned by drFit() function.
#' @param percent the inhibition ratio to be searched against. A vector between 0 and 1.
#'  Corresponding response is 1 - percent
#' @param log.d whether to return log10(dose) or the raw dose. Default is set to TRUE.
#' @param interpolation whether to use interpolation to estimate IC values. In this case,
#'  the computed IC values will be bound by the observed dosages.
#' @param stepLen step length to construct equally spaced intervals during interpolation.
#'  Only used when interpolation = TRUE.
#' @param lower lower bound in root search. Default is min(c(0, min(dose))) where dose is the
#'  observed dose levels from drFit@@originalDat.
#' @param upper upper bound in root search. Default is max(dose)*1e6 where dose is the observed
#'  dose levels from drFit@@originalDat.
#' @param ... other optimization parameters to be passed to uniroot().
#' @return A named vector giving the IC values at specified percentiles.
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit}}
#' @author Kevin R Coombes (\email{kcoombes@@mdanderson.org}), Pan Tong (\email{nickytong@@gmail.com})
#' @export
#' @examples
#' data(ryegrass, package = 'drc') # use the ryegrass data from drc package
#' fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha = 0.01, fitCtr = FALSE)
#' computeIC(fit.sigEmax, percent = seq(0, 1, by = 0.1), log.d = FALSE)	
computeIC <- function(drFit, percent = 0.50, log.d = TRUE, interpolation = TRUE,
                      stepLen = NA, lower, upper, ...) {
  res <- findDoseGivenResponse(drFit, response = 1-percent, log.d = log.d,
                               interpolation = interpolation, stepLen = stepLen,
                               lower = lower, upper = upper, ...)
  names(res) <- paste('IC', percent*100, sep = '')
  return(res)
}


#' compute variance of IC value at each inhibition percentage (original dose scale)
#'
#' Notice that this does not allow log.d option for log10 dose; it also does not allow
#'  interpretation/extrapolation. When a response is not achieved, the variance estimate would
#'  be NA indicating not achieved. 
#'
#' @param drFit A drFit object as returned by drFit() function.
#' @param percent the inhibition ratio to be searched against. A vector between 0 and 1.
#'  Corresponding response is 1-percent
#' @return a vector of estimated variance of the doses achieving each response level
#' @importFrom drc ED
#' @importFrom stats vcov
#' @importFrom foreach foreach
#' @export
computeICvariance <- function(drFit, percent = 0.50){
  tag <- drFit@tag
  fitObj <- drFit@fit
  modelName <- drFit@modelName
  response <- 1-percent
  res <- tryCatch(
    {
      if(tag == "drc"){ # using model fitted from drc
        dfp <- length(fitObj$fit$par) #  number of parameters
        # two columns: dose and variance
        tmp <- drc::ED(fitObj, respLev = response, interval = "delta",
                        type = "absolute", display = FALSE)
        varEst <- (tmp[, 2])^2 
      }
      if(tag == "DoseFinding"){ 
        coef.est <- coef(fitObj)
        cov.est <- vcov(fitObj)
        dfp <- length(coef.est)
        # only work for specific DoseFinding models
        # NA can be possible when y>1 all over for all models: IAI itself is NA since no dose will achieve y<1 in this case!
        if(modelName == 'linear'){
          # iteract through each response due to matrix operation
          varEst <- foreach(i = 1:length(response), .combine = 'c') %do% {
            # ic <- (response[i]-coef.est[1])/coef.est[2]
            der.coef <- c(-1/coef.est[2], -(response[i]-coef.est[1])/coef.est[2]^2)
            as.vector(der.coef%*%cov.est%*%der.coef)
          }
        } else if(modelName == 'sigEmax') {
          varEst <- foreach(i = 1:length(response), .combine = 'c') %do% {
            ic <- coef.est[3]/((coef.est[2]/(response[i]-coef.est[1])-1)^(1/coef.est[4]))
            der.coef <- c(coef.est[2]/(coef.est[4]*(coef.est[2]-response[i]+coef.est[1])*
                                         (response[i]-coef.est[1])),
                          -1/(coef.est[4]*(coef.est[2]-response[i]+coef.est[1])), 1/coef.est[3],
                          log(coef.est[2]/(response[i]-coef.est[1])-1)/coef.est[4]^2)
            as.vector(der.coef%*%cov.est%*%der.coef*ic^2)
          }
        } else if(modelName == 'logistic') {
          varEst <- foreach(i=1:length(response), .combine='c') %do% {
            der.coef <- c(coef.est[4]*coef.est[2]/((coef.est[2]-response[i]+coef.est[1])*
                                                     (response[i]-coef.est[1])),
                          -coef.est[4]/(coef.est[2]-response[i]+coef.est[1]), 1,
                          -log(coef.est[2]/(response[i]-coef.est[1])-1))
            as.vector(der.coef%*%cov.est%*%der.coef)
          }
        } else if(modelName == 'exponential') {
          varEst <- foreach(i=1:length(response), .combine='c') %do% {
            der.coef <-c(-coef.est[3]/(response[i]-coef.est[1]+coef.est[2]), 
                         -coef.est[3]*(response[i]-coef.est[1])/(coef.est[2]*
                                                                   (response[i]-coef.est[1]+coef.est[2])),
                         log((response[i]-coef.est[1]+coef.est[2])/coef.est[2]))
            as.vector(der.coef%*%cov.est%*%der.coef)
          }
        } else {
          varEst <- rep(NA, length(response))
        }
      }
      rr <- unname(varEst)
      rr[is.na(rr)] <- NA
      names(rr) <- paste('IC', percent*100, sep = '')
      rr
    }, error = function(cond){
      #message(sprintf('Diagnostic info: tag=%s; model=%s', tag, modelName))
      #message('Variane estimation failed with:')
      #message(cond)
      rr <- rep(NA, length(response))
      names(rr) <- paste('IC', percent*100, sep = '')
      return(rr)
    }, finally = {
      # nothing
    })
  res
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("i", "%do%"))


