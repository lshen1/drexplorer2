

getFixedRatioInfo <- function(d1, d2, e, tol = 0.001){
  ind1 <- d1!=0 & d2==0
  ind2 <- d1==0 & d2!=0
  ind12 <- d1!=0 & d2!=0 # where the dose differs
  ratio_d <- d2[ind12]/d1[ind12] # pointwise ratio
  isRayDesign <- all.equal(ratio_d, rep(mean(ratio_d), length(ratio_d)), tolerance = tol) # TRUE or a string
  if(isRayDesign != TRUE) isRayDesign <- FALSE # make string as false
  if(isRayDesign != TRUE){
    d2.d1 <- NA
    warning('This is not a fixed ratio design!')
  } else {
    d2.d1 <- mean(ratio_d) # this is not necessary for fixed ratio design
  }
  return(list(drMat = data.frame(d1=d1, d2=d2, e=e), isRayDesign = isRayDesign, 
              d1 = d1[ind1], e1 = e[ind1], d2 = d2[ind2], e2 = e[ind2],
              d12 = d1[ind12] + d2[ind12], e12 = e[ind12], ind12 = ind12,
              d2.d1 = d2.d1, d2.d1_avail = ratio_d))
}


#' convert d12 = d1+d2 to projected dose sqrt(d1^2+d2^2) in the direction of (d1, d2)
#'  given fixed ratio d2/d1
#'  
#' @param d12 d12 = d1+d2 to projected dose sqrt(d1^2+d2^2) in the direction of (d1, d2).
#' @param d2.d1 fixed ratio d2/d1.
d12_to_dproj <- function(d12, d2.d1){
  r <- d2.d1
  return(d12*sqrt(r^2+1)/(1+r))
}

#' convert d1 to projected dose sqrt(d1^2+d2^2) given fixed ratio d2/d1
#' 
#' @param d1 d1 to projected dose sqrt(d1^2+d2^2).
#' @param d2.d1 fixed ratio d2/d1.
d1_to_dproj <- function(d1, d2.d1){
  return(d1*sqrt(d2.d1^2+1))
}

#' convert d2 to projected dose sqrt(d1^2+d2^2) given fixed ratio d2/d1
#'
#' @param d2 d1 to projected dose sqrt(d1^2+d2^2). 
#' @param d2.d1 fixed ratio d2/d1.
#' @param islogd loged or not. Default is set to FLASE.
d2_to_dproj <- function(d2, d2.d1, islogd = FALSE){
  if(islogd) {
    d2 <- 10^d2
    res <- d2*sqrt(d2.d1^2+1)/d2.d1 #  return original dose
  } else {
    res <- d2*sqrt(d2.d1^2+1)/d2.d1		
  }
  return(res)
}


#' convert projected dose sqrt(d1^2+d2^2) to d1 given fixed ratio d2/d1
#'
#' @param dproj projected dose sqrt(d1^2+d2^2) to d1. 
#' @param d2.d1 fixed ratio d2/d1.
#' @param islogd loged or not. Default is set to FLASE.
dproj_to_d1 <- function(dproj, d2.d1, islogd = FALSE){
  if(islogd) {
    dproj <- 10^dproj
    res <- dproj/sqrt(d2.d1^2+1)
  } else {
    res <- 	dproj/sqrt(d2.d1^2+1)
  }
  return(res)
}


#' convert projected dose sqrt(d1^2+d2^2) to d2 given fixed ratio d2/d1
#'
#' @param dproj projected dose sqrt(d1^2+d2^2) to d2. 
#' @param d2.d1 fixed ratio d2/d1.
#' @param islogd loged or not. Default is set to FLASE.
dproj_to_d2 <- function(dproj, d2.d1, islogd = FALSE){
  if(islogd) {
    dproj <- 10^dproj
    res <- dproj*d2.d1/sqrt(d2.d1^2+1)
  } else {
    res <- dproj*d2.d1/sqrt(d2.d1^2+1)
  }
  return(res)	
}



getBliss <- function(y1, y2){
  # when y1>1, the yBliss is larger than individual one, obsurd!
  if(length(y2) == 1){
    # this is for constant dose design where B@x and thus y2 has a single value
    y2 <- rep(y2, length(y1))
  }
  res <- pmin(y1, 1)*pmin(y2, 1)
  return(res)
}


getHSA <- function(y1, y2){
  # when y1>1, the yBliss is larger than individual one, obsurd!
  if(length(y2) == 1){
    # this is for constant dose design where B@x and thus y2 has a single value
    y2 <- rep(y2, length(y1))
  }
  res <- pmin(y1, y2)
  return(res)
}


boundAUC <- function(AUC, ub = 1){
  AUC[AUC > ub] <- ub
  return(AUC)
}



# This is to calculate values in synergy plot table, which is useful to calculate synergy and
#  log volume values.
# The interpretation: it is the extra kill percent where the reference is the most extreme
#  observation under a given alpha; so this is conservative extra percent killing.
getSynergyPlotValues <- function(avg, mean, sd, alpha = 0.05){
  Zalpha <- qnorm(1-alpha/2, 0, 1) # 1.96 for alpha=0.05, 2.57 for alpha=0.01
  res1 <- (avg-(mean+Zalpha*sd))
  res2 <- (avg-(mean-Zalpha*sd))
  res <- rep(0, length(res1)) 
  i1 <- which(avg > mean+Zalpha*sd)
  res[i1] <- res1[i1]
  i2 <- which(avg<mean-Zalpha*sd)
  res[i2] <- res2[i2]
  res <- forceZero(res)
  res[is.na(res1) | is.na(res2)] <- NA # 2017/0501: may contain NA since not enough data from lab as they deleted some
  return(res)
}	




fitHSA_ray_ <- function(d1, d2, e, tol = 0.05, strict = TRUE,
                        name1 = 'Drug A', name2 = 'Drug B',
                        dmin = NULL, dmax = NULL, islogd = TRUE){
  frInfo <- getFixedRatioInfo(d1, d2, e, tol = tol)
  if(!frInfo$isRayDesign & strict) {
    rs <- round(frInfo$d2.d1_avail, 2)
    ratios <- safeCSV(sort(unique(rs)))
    message <- sprintf('\nNot fixed ratio design!\n Drug1: %s; Drug2: %s\tRatios observed: %s; Ratio majority: %s\n',
                       name1, name2, ratios, names(sort(table(rs), decreasing = T))[1])
    stop(message)
  }
  d2.d1 <- frInfo$d2.d1
  d2.d1_avail <- frInfo$d2.d1_avail
  ind1 <- d1!=0 & d2==0
  ind2 <- d1==0 & d2!=0
  ind12 <- d1!=0 & d2!=0 # where the dose differs
  ### specify models
  models <- iaiModels
  dat1 <- data.frame(Dose = d1[ind1], Response = e[ind1])
  dat2 <- data.frame(Dose = d2[ind2], Response = e[ind2])
  dat12 <- data.frame(Dose = d1[ind12]+d2[ind12], Response = e[ind12]) 
  # notice the combined dose is used for fitting combination data
  # this is consistent to Dianne's code
  fitdat <- function(dat, drug){
    resPrepDR <- prepDRdat(dat, alpha = 1, fitCtr = FALSE, standardize = FALSE)
    info_experimentSetup <- prepDRdat2expInfo(resPrepDR)
    resL <- fitOneExp(dat = dat, drug = drug, cellLine = '', models = models,
                      percentile = seq(0.05, 0.95, 0.05), plot = FALSE, fitCtr = FALSE,
                      transparency = 0.95, standardize = FALSE, interpolation = FALSE, unit = '')
    resL$ICx <- data.frame(resL$ICx)
    resL$info_experimentSetup <- info_experimentSetup
    resL
  }
  dat12 <- mutate(dat12, Dose = d12_to_dproj(Dose, d2.d1 = frInfo$d2.d1_avail))
  ###
  ### fit drugs: drug1: y~d1; drug2: y~d2; combo: y~projected dose
  ###
  ###
  ### plot: combo (projected dose, f(projected dose))
  ### drug 1: (projected dose, f1(d1))
  ### drug 2: (projected dose, f2(d2))
  ###
  resL_1 <- fitdat(dat1, name1)
  resL_2 <- fitdat(dat2, name2)
  dat1 <- mutate(dat1, dproj = d1_to_dproj(d1 = Dose, d2.d1 = d2.d1_avail))
  dat2 <- mutate(dat2, dproj = d2_to_dproj(d2 = Dose, d2.d1 = d2.d1_avail))
  resL_12 <- fitdat(dat12, str_c(name1, name2, sep = '_'))
  df12 <- data.frame(dat12, d1 = d1[ind12], d2 = d2[ind12], d12 = d1[ind12]+d2[ind12],
                     dproj = dat12$Dose)
  ##
  ## AUC calculation
  ##
  # assumption: the predict function is trained on original dose
  # The AUC is always calculated at projected dose axis. So to fulfill this, given dproj, convert it to the same dose axis as fit (to d1 for A, d2 for B or dproj for combo) and get the response
  # then do integration; Note, this can be simplified: all fits are forced to on dproj, then no dose conversion needed
  # the rationale: predict response at each dose; apply integration;
  # AUC can be calculated based on log10 dr curve or original dr curve (in this case, absolute AUC is smaller, scaled AUC is also smaller usually); Notice that this only changes the integration range/scale (x) without affecting response (y)
  # However, the user should supply dmin/dmax in log10 scale if islogd==TRUE (for log10 AUC) and dmin/dmax in original scale for islogd=FALSE (for untransformed AUC)
  # according to the Nature paper, AUC is defined on the log10 scale; this is also more intuitive since it matches dr curve which is in log10 scale
  # islogd: this makes it possible to calculate AUC either on log10dose or original dose. However, the user
  #      should make sure dmin and dmax is on the same scale (take log10 correspondingly)
  ## always supply dproj, but need to convert to d1 or d2 for fit1 and fit2
  fint_12 <- function(fit, dproj, islogd = islogd) {
    if(islogd) {
      dd <- 10^(dproj)
    } else {
      dd <- dproj
    }
    res <- predict(fit, newData=dd)
    res  
  }
  fint_1 <- function(fit, dproj, islogd = islogd) {
    if(islogd) {
      dd <- 10^(dproj)
    } else {
      dd <- dproj
    }
    res <- predict(fit, newData = dproj_to_d1(dd, d2.d1 = d2.d1))
    res  
  }
  fint_2 <- function(fit, dproj, islogd = islogd) {
    if(islogd) {
      dd <- 10^(dproj)
    } else {
      dd <- dproj
    }
    res <- predict(fit, newData = dproj_to_d2(dd, d2.d1 = d2.d1))
    res  
  }
  fint_HSA <- function(fit1, fit2, dproj, islogd = islogd) {
    if(islogd) {
      dd <- 10^(dproj)
    } else {
      dd <- dproj
    }
    y1 <- predict(fit1, newData = dproj_to_d1(dd, d2.d1 = d2.d1))
    y2 <- predict(fit2, newData = dproj_to_d2(dd, d2.d1 = d2.d1))
    res <- pmin(y1, y2) 
    res
  }
  fint_Bliss <- function(fit1, fit2, dproj, islogd = islogd) {
    if(islogd) {
      dd <- 10^(dproj)
    } else {
      dd <- dproj
    }
    y1 <- predict(fit1, newData = dproj_to_d1(dd, d2.d1 = d2.d1))
    y2 <- predict(fit2, newData = dproj_to_d2(dd, d2.d1 = d2.d1))
    res <- getBliss(y1, y2)
    res
  }
  f_ref <- function(dproj, reference = 1, islogd = islogd) {
    if(islogd) {
      dd <- 10^(dproj)
    } else {
      dd <- dproj
    }
    res <- rep(reference, length(dproj))
    res 
  }
  res <- list(fit1 = resL_1$fits[[resL_1$indBest]], 
              fit2 = resL_2$fits[[resL_2$indBest]], 
              fit12 = resL_12$fits[[resL_12$indBest]],
              dat1 = dat1,
              dat2 = dat2,
              dat12 = df12,
              d2.d1 = frInfo$d2.d1,
              d2.d1_avail = frInfo$d2.d1_avail,
              drug1 = name1,
              drug2 = name2,
              drug12 = sprintf('%s+%s', name1, name2))
  datMacsyn <- ddply(res$dat12, .(dproj), summarise,
                     CPavg = mean(Response, na.rm = T),
                     CPstd = sd(Response, na.rm = T),
                     d1 = mean(d1, na.rm = T), d2 = mean(d2, na.rm = T))
  # add predicted value and bliss value
  datMacsyn <- mutate(datMacsyn,
                      y1 = predict(res$fit1, d1),
                      y2 = predict(res$fit2, d2),
                      yBliss = getBliss(y1, y2),
                      yHSA = getHSA(y1, y2))
  # add InhibitionAvg
  datMacsyn_HSA <- datMacsyn
  datMacsyn <- mutate(datMacsyn, InhibitionAvg = 1-CPavg,
                      InhibitionAvgPercent = InhibitionAvg*100, 
                      InhibitionStdPercent = CPstd*100, AdditiveInhibitionPercent = 100*(1-yBliss),
                      extra_kill_percent = InhibitionAvgPercent - AdditiveInhibitionPercent)
  datMacsyn_HSA <- mutate(datMacsyn_HSA, InhibitionAvg = 1-CPavg,
                          InhibitionAvgPercent = InhibitionAvg*100, 
                          InhibitionStdPercent = CPstd*100,
                          AdditiveInhibitionPercent = 100*(1-yHSA),
                          extra_kill_percent = InhibitionAvgPercent - AdditiveInhibitionPercent)
  # add P value
  datMacsyn <- mutate(datMacsyn, p = getNormPvec(InhibitionAvgPercent, AdditiveInhibitionPercent, InhibitionStdPercent),
                      synergy_plot = forceZero(extra_kill_percent),
                      synergy_plot_95 = getSynergyPlotValues(avg = InhibitionAvgPercent, mean = AdditiveInhibitionPercent, sd = InhibitionStdPercent, alpha = 0.05),
                      synergy_plot_99 = getSynergyPlotValues(avg = InhibitionAvgPercent, mean = AdditiveInhibitionPercent, sd  =InhibitionStdPercent, alpha = 0.01),
                      synergy_plot_99.9 = getSynergyPlotValues(avg = InhibitionAvgPercent, mean = AdditiveInhibitionPercent, sd = InhibitionStdPercent, alpha = 0.001)
  )
  datMacsyn_HSA <- mutate(datMacsyn_HSA, p = getNormPvec(InhibitionAvgPercent, AdditiveInhibitionPercent, InhibitionStdPercent),
                          synergy_plot = forceZero(extra_kill_percent),
                          synergy_plot_95 = getSynergyPlotValues(avg = InhibitionAvgPercent, mean = AdditiveInhibitionPercent, sd = InhibitionStdPercent, alpha = 0.05),
                          synergy_plot_99 = getSynergyPlotValues(avg = InhibitionAvgPercent, mean = AdditiveInhibitionPercent, sd = InhibitionStdPercent, alpha = 0.01),
                          synergy_plot_99.9 = getSynergyPlotValues(avg = InhibitionAvgPercent, mean = AdditiveInhibitionPercent, sd = InhibitionStdPercent, alpha  =0.001)
  )
  res$datMacsyn <- datMacsyn
  res$datMacsyn_HSA <- datMacsyn_HSA
  if(is.null(dmin)) {
    dmin <- min(res$dat12$dproj)
    if(islogd) dmin <- log10(dmin)
  }
  if(is.null(dmax)) {
    dmax <- max(res$dat12$dproj)
    if(islogd) dmax <- log10(dmax)
  }
  tryError2NA <- function(res){
    if(class(res) == 'try-error') res <- NA
    res
  }
  AUC_ref <- tryError2NA(try(integrate(f_ref, lower=dmin, upper=dmax, islogd=islogd)$value, silent=TRUE))  
  AUC_12 <- tryError2NA(try(integrate(fint_12, fit=res$fit12, lower=dmin, upper=dmax, islogd=islogd)$value/AUC_ref  , silent=TRUE))
  AUC_HSA <- tryError2NA(try(integrate(fint_HSA, fit1=res$fit1, fit2=res$fit2, lower=dmin, upper=dmax, islogd=islogd)$value/AUC_ref, silent=TRUE))  
  AUC_Bliss <- tryError2NA(try(integrate(fint_Bliss, fit1=res$fit1, fit2=res$fit2, lower=dmin, upper=dmax, islogd=islogd)$value/AUC_ref, silent=TRUE))
  AUC_1 <- tryError2NA(try(integrate(fint_1, fit=res$fit1, lower=dmin, upper=dmax, islogd=islogd)$value/AUC_ref, silent=TRUE))
  AUC_2 <- tryError2NA(try(integrate(fint_2, fit=res$fit2, lower=dmin, upper=dmax, islogd=islogd)$value/AUC_ref, silent=TRUE))  
  ### BOUNDING
  AUC_12 <- boundAUC(AUC_12)
  AUC_1 <- boundAUC(AUC_1)
  AUC_2 <- boundAUC(AUC_2)
  AUC_HSA <- boundAUC(AUC_HSA)
  AUC_Bliss <- boundAUC(AUC_Bliss)
  AUC_ref <- boundAUC(AUC_ref)
  deltaAUC_HSA <- AUC_12-AUC_HSA
  deltaAUC_Bliss <- AUC_12-AUC_Bliss
  res$islogd <- islogd
  res$dmin <- dmin
  res$dmax <- dmax
  res$AUC_ref <- AUC_ref
  res$AUC_12 <- AUC_12
  res$AUC_HSA <- AUC_HSA
  res$AUC_Bliss <- AUC_Bliss
  res$deltaAUC_HSA <- deltaAUC_HSA
  res$deltaAUC_Bliss <- deltaAUC_Bliss
  res$AUC_1 <- AUC_1
  res$AUC_2 <- AUC_2
  return(res)
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Dose", "dproj", "y1", "y2",
                                                        "CPavg", "CPstd",
                                                        "yBliss", "yHSA",
                                                        "InhibitionAvg",
                                                        "InhibitionAvgPercent",
                                                        "InhibitionStdPercent",
                                                        "extra_kill_percent",
                                                        "AdditiveInhibitionPercent"))


#' bootstrap of fitHSA_ray to obtain variation of deltaAUC
#' 
#' @param d1 drug 1 dose in original scale.
#' @param d2 drug 2 dose in original scale.
#' @param e response at the dose pair. 
#' @param tol tolerance in declaring fixed ratio.
#' @param strict whether to stop if not fixed ratio; this is just for visualization;
#'  the calculations of synergy is wrong since d12 = d1+d2 is meaningless for variable ratio! 
#'  At the end, this is a waste of time: it would not make sense and impossible to do this!
#' @param name1 drug 1 name.
#' @param name2 drug 2 name.
#' @param dmin start position to calcualte AUC at projected dose, in original scale,
#'  default is minimum projected dose.
#' @param dmax end position to calcualte AUC at projected dose, in original scale,
#'  default is maximum projected dose.
#' @param islogd whether to compute AUC at log10 dose; if TRUE, dmin and dmax must also be
#'  log10 scale.
#' @param bootstrapB number of boostrap iterations.
#' @param nrep how many replicates are sampled with replacement for each dose combination.
#'  nrep=3 will make it much narrower CI band.
#' @param seed random seed. Dafault is set to 100.
#' 
#' @importFrom stats quantile  
fitHSA_ray_boot <- function(d1, d2, e, tol = 0.05, strict = TRUE,
                            name1 = 'Drug A', name2 = 'Drug B',
                            dmin = NULL, dmax = NULL, islogd = TRUE,
                            bootstrapB = 10, nrep = 3, seed = 100){
  set.seed(seed)
  myd <- data.frame(d1 = d1, d2 = d2, e = e)
  myd <- mutate(myd, d1_d2 = str_c(d1, d2, sep = '_'))
  resBoot <- foreach(b = 1:bootstrapB) %do% {
    # sample data
    myd_s <- ddply(myd, .(d1_d2), function(x) {
      N <- nrow(x)
      ss <- sample(1:N, size = nrep, replace = TRUE)
      x[ss, ]
    })
    # fit
    tp <- with(myd_s, fitHSA_ray_(d1 = d1, d2 = d2, e = e, tol = tol, strict = strict,
                                  name1 = name1, name2 = name2, dmin = dmin, dmax = dmax,
                                  islogd = islogd))
    tp[16:23]
  }
  deltaAUC_Bliss_CI <- quantile(sapply(resBoot, function(x) x$deltaAUC_Bliss), c(0.025, 0.975))
  deltaAUC_HSA_CI <- quantile(sapply(resBoot, function(x) x$deltaAUC_HSA), c(0.025, 0.975))
  return(list(resBoot = resBoot,
              deltaAUC_Bliss_CI = deltaAUC_Bliss_CI,
              deltaAUC_HSA_CI = deltaAUC_HSA_CI))
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("d1_d2"))




#' Fit HSA and Bliss model given a fixed ratio design with input (d1, d2, e)
#' 
#' @param d1 drug 1 dose in original scale.
#' @param d2 drug 2 dose in original scale.
#' @param e response at the dose pair. 
#' @param tol tolerance in declaring fixed ratio.
#' @param strict whether to stop if not fixed ratio; this is just for visualization;
#'  the calculations of synergy is wrong since d12 = d1+d2 is meaningless for variable ratio! 
#'  At the end, this is a waste of time: it would not make sense and impossible to do this!
#' @param name1 drug 1 name.
#' @param name2 drug 2 name.
#' @param dmin start position to calcualte AUC at projected dose, in original scale,
#'  default is minimum projected dose.
#' @param dmax end position to calcualte AUC at projected dose, in original scale,
#'  default is maximum projected dose.
#' @param islogd whether to compute AUC at log10 dose; if TRUE, dmin and dmax must also be
#'  log10 scale.
#' @param bootstrapB number of boostrap iterations.
#' @param nrep how many replicates are sampled with replacement for each dose combination.
#'  nrep=3 will make it much narrower CI band.
#' @param seed random seed. Dafault is set to 100.
#' 
#' @return a list
#' @export
fitHSA_ray <- function(d1, d2, e, tol = 0.05, strict = TRUE,
                       name1 = 'Drug A', name2 = 'Drug B', 
                       dmin = NULL, dmax = NULL, islogd = TRUE,
                       bootstrapB = 0, nrep = 3, seed = 100){
  # this calls the core
  res <- fitHSA_ray_(d1 = d1, d2 = d2, e = e, tol = tol, strict = strict,
                     name1 = name1, name2 = name2,
                     dmin = dmin, dmax = dmax, islogd = islogd)
  if(bootstrapB > 1){
    resB <- fitHSA_ray_boot(d1 = d1, d2 = d2, e = e, tol = tol, strict = strict,
                            name1 = name1, name2 = name2,
                            dmin = dmin, dmax = dmax, islogd = islogd,
                            bootstrapB = bootstrapB, nrep = nrep, seed = seed)
    # see https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Methods_for_bootstrap_confidence_intervals
    res$deltaAUC_Bliss_CI <- 2*res$deltaAUC_Bliss - rev(resB$deltaAUC_Bliss_CI)
    res$deltaAUC_HSA_CI_CI <- 2*res$deltaAUC_Bliss - rev(resB$deltaAUC_HSA_CI) # may not cover the mean estimate!!!
    res$bootstrap <- resB
  }
  res$design <- 'fixedRatio'
  return(res)
}





#' plot the ray design from fit object returned by fitHSA_ray
#' 
#' @param fit the fitted result from fitHSA_ray().
#' @param cols color corresponding to mixture, drugA, drugB.
#' @param lwds line width.
#' @param xlim x-axis limit.
#' @param ylim y-axis limit.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param main figure title.
#' @param tag title tag.
#' @param type plot types.
#' @param legend whether to show legend.
#' @param legend.cex cex of the legend.
#' @param plotSetup whether initiate a new figure; set to False when just to add curve as lines.
#' @param digit Default is set to 3.
#' 
#' @export
plotHSA_ray <- function(fit, cols = c(1, 2, 3), lwds = 3,
                        xlim = NULL, ylim = NULL, xlab = NULL, ylab = 'Relative viability', 
                        main = NULL, tag = NULL, type = c('line', 'HSA', 'Bliss', 'se'),
                        legend = TRUE, legend.cex = 0.6, plotSetup = TRUE, digit = 3){
  design <- fit$design
  if(design == 'fixedRatio'){
    if(is.null(xlab)) xlab <- 'Projected dose'
    fit1 <- fit$fit1
    fit2 <- fit$fit2
    fit12 <- fit$fit12
    dat12 <- fit$dat12
    dat1 <- fit$dat1
    dat2 <- fit$dat2
    d2.d1 <- fit$d2.d1
    gg1 <- format_grid(dat12$d1)
    datp <- data.frame(d1 = gg1$xGrid)
    datp <- mutate(datp, d2 = d2.d1*d1, dproj = d12_to_dproj(d1+d2, d2.d1 = d2.d1),
                   y1 = predict(fit1, d1),
                   y2 = predict(fit2, d2),
                   y12 = predict(fit12, dproj),
                   yHSA = pmin(y1, y2),
                   yBliss = getBliss(y1, y2))
    if(is.null(ylim)) ylim <- c(0, range(pretty(c(dat12$Response, dat1$Response, dat2$Response)))[2])
    if(is.null(xlim)) xlim <- range(pretty(log10(datp$dproj)))
    if(is.null(main)){
      main <- sprintf('AUC: mixture=%.3f, HSA=%.3f, Bliss=%.3f\nDiff AUC: HSA=%.3f, Bliss=%.3f',
                      fit$AUC_12, fit$AUC_HSA, fit$AUC_Bliss, fit$deltaAUC_HSA, fit$deltaAUC_Bliss)
    }
    if(!is.null(tag)) main <- sprintf('%s\n%s', tag, main)
    if('line' %in% type){
      if(plotSetup){
        # plot setup; disable this if multiple plots are to be shown on one figure
        op1 <- par(mar = c(8, 7, 5, 2) + 0.1)
        plot(log10(datp$dproj), datp$y12, col = cols[1], lwd = lwds,
             ylim = ylim, xlim = xlim, type = 'l', xlab = '', ylab = ylab, main = main, xaxt = 'n')
      } else {
        lines(log10(datp$dproj), datp$y12, col = cols[1], lwd = lwds)
      }
      lines(log10(datp$dproj), datp$y1, col = cols[2], lwd = lwds)
      lines(log10(datp$dproj), datp$y2, col = cols[3], lwd = lwds)
    }
    if('points' %in% type){
      points(log10(dat12$dproj), dat12$Response, col = cols[1])
      points(log10(dat1$dproj), dat1$Response, col = cols[2])
      points(log10(dat2$dproj), dat2$Response, col = cols[3])
    }	
    getSmry <- function(dat) {
      ddply(dat, .(dproj), function(x) data.frame(Mean = mean(x$Response),
                                                  SE = sd(x$Response)/sqrt(nrow(x))))
    }
    if('se' %in% type){
      with(getSmry(dat12), points(log10(dproj), Mean, cex = 0.5, col = cols[1]))
      with(getSmry(dat1), points(log10(dproj), Mean, cex = 0.5, col = cols[2]))
      with(getSmry(dat2), points(log10(dproj), Mean, cex = 0.5, col = cols[3]))
      with(getSmry(dat12), arrows(log10(dproj), Mean-SE, log10(dproj), Mean+SE,
                                  col = cols[1], code = 3, length = 0.02, angle = 90))
      with(getSmry(dat1), arrows(log10(dproj), Mean-SE, log10(dproj), Mean+SE,
                                 col = cols[2], code = 3, length = 0.02, angle = 90))
      with(getSmry(dat2), arrows(log10(dproj), Mean-SE, log10(dproj), Mean+SE,
                                 col = cols[3], code = 3, length = 0.02, angle = 90))
    }
    if('HSA' %in% type){
      lines(log10(datp$dproj), datp$yHSA, col = 'blue', lty = 2, lwd = lwds)
    }	
    if('Bliss' %in% type){
      lines(log10(datp$dproj), datp$yBliss, col = 'purple', lty=2, lwd = lwds)
    }	
    ltext <- c(fit$drug12, fit$drug1, fit$drug2)
    legend_col <- cols
    legend_lty <- c(1, 1, 1)
    legend_lwd <- c(lwds, lwds, lwds)
    if('HSA' %in% type){
      ltext <- c(ltext, 'HSA')
      legend_col <- c(legend_col, 'blue')
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, lwds)
    }
    if('Bliss' %in% type){
      ltext <- c(ltext, 'Bliss')
      legend_col <- c(legend_col, 'purple')
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, lwds)
    }
    if(legend){
      # legend('bottomleft', legend = ltext, col = legend_col,
      #        lty = legend_lty, lwd = legend_lwd,
      #        inset = c(-0.275, 1), xpd = TRUE, horiz = FALSE, bty = 'n')
      legend('bottomright', legend = ltext, col = legend_col,
             lty = legend_lty, lwd = legend_lwd,
             inset = c(0,1), horiz = TRUE, cex = legend.cex, xpd = TRUE, bty = 'n')
    }
    # add backup axis
    if(plotSetup){
      mtext(xlab, side = 1, line = 2)
      at_proj <- axTicks(1)
      at_d1 <- dproj_to_d1(at_proj, d2.d1 = fit$d2.d1, islogd = TRUE)
      at_d2 <- dproj_to_d2(at_proj, d2.d1 = fit$d2.d1, islogd = TRUE)
      # for proj. dose
      axis(1, at = at_proj, labels = round(10^at_proj, digits = digit), tick = T,
           lwd = 1, line = 0, col = 'black', col.axis = 'black', col.ticks = 'black', lty = 1)
      # for d1
      axis(1, at = at_proj, labels = round(at_d1, digits = digit), tick = F,
           lwd = 1, line = 2.5, col = cols[2], col.axis = cols[2], col.ticks = cols[2], lty = 2)
      axis(1, at = at_proj, labels = round(at_d1, digits = digit), tick = F,
           lwd = 1, line = 2.5, col = cols[2], col.axis = cols[2], col.ticks = cols[2], lty = 2)
      # for d2
      axis(1, at = at_proj, labels = round(at_d2, digits = digit), tick = F,
           lwd = 1, line = 5, col = cols[3], col.axis = cols[3], col.ticks = cols[3], lty = 2)
      # 
      mtat <- par("usr")[1]-0.15*diff(par("usr")[1:2])
      mtext(sprintf("%s dose", fit$drug1),
            side = 1, line = 2.5+1, cex.lab = 1, las = 1, at = mtat, col = cols[2])
      mtext(sprintf("%s dose", fit$drug2),
            side = 1, line = 5+1, cex.lab = 1, las = 1, at = mtat, col = cols[3])
    }
  } else {
    #################################################################################################################################	
    # constant design
    #################################################################################################################################	
    drug1 <- fit$drug1
    fit1 <- fit$fit1
    fit12 <- fit$fit12
    dat12 <- fit$dat12
    dat1 <- fit$dat1
    dat2 <- fit$dat2
    yConstant <- fit$yConstant
    if(is.null(xlab)) xlab <- sprintf('%s dose', drug1)
    gg1 <- format_grid(dat12$Dose)
    datp <- data.frame(d1 = gg1$xGrid)
    datp <- mutate(datp,
                   y1 = predict(fit1, d1),
                   y12 = predict(fit12, d1),
                   yHSA = getHSA(y1, yConstant),
                   yBliss = getBliss(y1, yConstant)
    )
    if(is.null(ylim)) ylim <- c(0, range(pretty(c(dat12$Response, dat1$Response, dat2$Response)))[2])
    if(is.null(xlim)) xlim <- range(pretty(log10(datp$d1)))
    if(is.null(main)){
      main <- sprintf('AUC: mixture=%.3f, HSA=%.3f, Bliss=%.3f\nDiff AUC: HSA=%.3f, Bliss=%.3f',
                      fit$AUC_12, fit$AUC_HSA, fit$AUC_Bliss, fit$deltaAUC_HSA, fit$deltaAUC_Bliss)
    }
    if(!is.null(tag)) main <- sprintf('%s\n%s', tag, main)
    if('line' %in% type){
      if(plotSetup){
        # plot setup; disable this if multiple plots are to be shown on one figure
        op1 <- par(mar = c(8, 7, 4, 2) + 0.1)
        plot(log10(datp$d1), datp$y12, col = cols[1], lwd = lwds, ylim = ylim, xlim = xlim,
             type = 'l', xlab = '', ylab = ylab, main = main, xaxt = 'n')
      } else {
        lines(log10(datp$d1), datp$y12, col = cols[1], lwd = lwds)
      }
      lines(log10(datp$d1), datp$y1, col = cols[2], lwd = lwds)
      abline(h = yConstant, col = cols[3], lty = 5)
    }
    if('points' %in% type){
      points(log10(dat12$Dose), dat12$Response, col = cols[1])
      points(log10(dat1$Dose), dat1$Response, col = cols[2])
      points(rep(xlim[1], nrow(dat2)), dat2$Response, col = cols[3])
    }	
    getSmry <- function(dat) {
      ddply(dat, .(Dose), function(x) data.frame(Mean = mean(x$Response),
                                                 SE = sd(x$Response)/sqrt(nrow(x))))
    }
    if('se' %in% type){
      with(getSmry(dat12), points(log10(Dose), Mean, cex = 0.5, col = cols[1]))
      with(getSmry(dat1), points(log10(Dose), Mean, cex = 0.5, col = cols[2]))
      with(getSmry(dat12), arrows(log10(Dose), Mean-SE, log10(Dose), Mean+SE,
                                  col = cols[1], code = 3, length = 0.02, angle = 90))
      with(getSmry(dat1), arrows(log10(Dose), Mean-SE, log10(Dose), Mean+SE,
                                 col = cols[2], code = 3, length = 0.02, angle = 90))
    }
    if('HSA' %in% type){
      lines(log10(datp$d1), datp$yHSA, col = 'blue', lty = 2, lwd = lwds)
    }	
    if('Bliss' %in% type){
      lines(log10(datp$d1), datp$yBliss, col = 'purple', lty = 2, lwd = lwds)
    }	
    ltext <- c(fit$drug12, fit$drug1, sprintf('%s@%.2f', fit$drug2, dat2[1, 'Dose']))
    legend_col <- cols
    legend_lty <- c(1, 1, 5)
    legend_lwd <- c(lwds, lwds, lwds)
    if('HSA' %in% type){
      ltext <- c(ltext, 'HSA')
      legend_col <- c(legend_col, 'blue')
      legend_lty <- c(legend_lty, lwds)
      legend_lwd <- c(legend_lwd, lwds)
    }
    if('Bliss' %in% type){
      ltext <- c(ltext, 'Bliss')
      legend_col <- c(legend_col, 'purple')
      legend_lty <- c(legend_lty, lwds)
      legend_lwd <- c(legend_lwd, lwds)
    }
    if(legend){
      # legend('bottomleft', legend = ltext, col = legend_col,
      #        lty = legend_lty, lwd = legend_lwd,
      #        inset = c(-0.275, 1), xpd = TRUE, horiz = FALSE, bty = 'n')
      legend('bottomright', legend = ltext, col = legend_col,
             lty = legend_lty, lwd = legend_lwd,
             inset = c(0,1), horiz = TRUE, cex = legend.cex, xpd = TRUE, bty = 'n')
    }
    # add backup axis
    if(plotSetup){
      mtext(xlab, side = 1, line = 2)
      at_proj <- axTicks(1)
      labels <- sapply(at_proj,function(i) as.expression(bquote(10^ .(i))))
      axis(1, at = at_proj, labels = labels, tick = T,
           lwd = 1,line = 0, col = 'black', col.axis = 'black', col.ticks = 'black', lty = 1)
    }
  }
}
## In the same source file (to remind you that you did it) add:
if(getRversion() >= "2.15.1")  utils::globalVariables(c("d1", "d2"))

