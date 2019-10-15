
check.drMat <- function(drMat){
  if(missing(drMat)) stop("Input drMat missing!") 
  if(!is.matrix(drMat) & !is.data.frame(drMat)) stop('Input drMat must be a matrix or a data frame!')
  if(ncol(drMat) > 2) warning('More than 2 columns are present in drMat. Only the first 2 columns will be used!')
  if(nrow(drMat) <= 3) stop('At least 3 rows are needed for drMat!')
  if(!all(apply(drMat, 2, is.numeric))) stop('All elements in drMat must be numeric!')
}

check.alpha <- function(alpha) {
  if(length(alpha) > 1) stop('More than 1 alpha levels specified. Only 1 value can be specified!')
  if(!(alpha %in% c(0.01, 0.05, 1))) stop('Alpha can only be 0.01, 0.05 or 1!')
}

noNA <- function(dat, returnIndex = FALSE){
  sel <- stats::complete.cases(dat)
  if(returnIndex)
    return(sel)
  if(is.null(dim(dat))) { 
    res <- dat[sel]
  } else {
    res <- dat[sel, ]
  }
  return(res)
}


truncByLimit <- function(x, Lower, Upper) {
  x[x < Lower] <- Lower
  x[x > Upper] <- Upper
  return(x)
}


truncate_effect <- function(e, min = 1e-4, max = 1-1e-4){
  e[e < min] <- min
  e[e > max] <- max
  return(e)
}


forceZero <- function(x, eps = 1e-10){
  x[abs(x) < eps] <- 0
  return(x)
}


safeCSV <- function(vec, empty = '', sep = ','){
  if(length(vec) > 0){
    csv <- str_c(noNA(vec), collapse = sep)
  } else {
    csv <- empty
  }	
  return(csv)
}

#' calculate P value assuming avg ~ Normal(mean, sd)
#' 
#' @param avg average.
#' @param mean mean.
#' @param sd standard deviation.
#' 
#' @importFrom stats pnorm
getNormP <- function(avg, mean, sd){
  if(avg > mean){
    return(1-pnorm(avg, mean, sd))
  } else {
    return(pnorm(avg, mean, sd))
  }
}

#' calculate P value assuming avg ~ Normal(mean, sd)
#' 
#' @param avg average.
#' @param mean mean.
#' @param sd standard deviation.
#'
#' @importFrom stats pnorm
getNormPvec <- function(avg, mean, sd){
  ind1 <- avg > mean
  res1 <- 1-pnorm(avg, mean, sd)
  res2 <- pnorm(avg, mean, sd)
  res <- res2
  i1 <- which(ind1)
  res[i1] <- res1[i1]
  return(res)
}


#' General getter function for attributes in an object
#'
#' @param obj object with attributes
#' @param what specify what attribute to get 
#' @return attributes extracted 
toget <- function(obj, what = 'avail'){
  if(what %in% names(attributes(obj))){
    res <- attributes(obj)[[what]]
  } else if (what == 'avail'){
    res <- names(attributes(obj))
  } else {
    res <- NULL
  }
  return(res)
}