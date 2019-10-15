
#'  Table III and IV published by Newman
#'
#' This dataset contains two 20-by-12 matrices corresponding to table III and table IV published
#'  by Newman, D. In particular, onePercentTab gives the quantiles at significance level of 0.01
#'  and fivePercentTab gives the quantiles at significance level of 0.01. The degree of freedom f
#'  and sample size n are used as row and column names for each table.
#'
#' @format Each of the two objects (onePercentTab, fivePercentTab)
#' is a matrix with 20 rows (degree of freedom f) and 12 columns (sample size n).
#' @references Newman, D. (1939). The distribution of range in samples from a normal population, 
#'  expressed in terms of an independent estimate of standard deviation. Biometrika, 31(1/2), 20-30.
#' @name NewmanTables
#' @aliases onePercentTab fivePercentTab
NULL





#'  ryegrass dataset from drc package
#'
#' This dataset is from the ryegrass data in drc package. It is a 2 column data frame where the
#'  conc column is the concentration of ferulic acid in mM and rootl column is the root length of
#'  ryegrass. 
#'
#' @references Inderjit and J. C. Streibig, and M. Olofsdotter (2002) Joint action of phenolic
#'  acid mixtures and its significance in allelopathy research, Physiologia Plantarum, 114,
#'  422-428, 2002.
#' @name ryegrass
"ryegrass"

#'  Example raw data for calculating GI50, TGI and LC50
#'
#' This dataset is an example data set where the response is between -1 and 1.
#'  It is a 2 column data frame where the Dose column is the concentration of drug 
#'  and Growth column is growth percentage as defined by the NCI60 method. 
#'
#' @references Holbeck, S. L., Collins, J. M., & Doroshow, J. H. (2010). Analysis of food
#'  and drug administration approved anticancer agents in the NCI60 panel of human tumor cell lines.
#'  Molecular cancer therapeutics, 9(5), 1451-1460.
#' @name datNCI60
"datNCI60"

#' The drug combination data between SCH66336 and 4-HPR. This is an all-possible combination design.
#'
#' This dataset is published in Lee et al, 2007, Journal of Biopharmaceutical Statistics, Table 1
#'
#' \itemize{
#'   \item schd SCH66336 dose
#'   \item hpr 4-HPR dose
#'   \item y1 response
#' }
#'
#' @format A data frame with 30 rows and 3 columns
#' @source \url{https://biostatistics.mdanderson.org/SoftwareDownload/}
#' @references Lee, J. Jack, et al. Interaction index and different 
#'  methods for determining drug interaction in combination therapy. Journal of biopharmaceutical
#'  statistics 17.3 (2007): 461-480.
#' @name nl22B2
"nl22B2"


#' The drug combination data of squamous cell carcinoma cells (UMSCC22B) surviving after 72 hours
#'  of treatment by single and combination dose levels of SCH66336 and 4-HPR. This is a fixed-ratio
#'  design.
#'
#' This dataset is published in Lee et al, 2009, Statistics in biopharmaceutical research, Table 2
#'
#' \itemize{
#'   \item dose1 SCH66336 dose
#'   \item dose2 4-HPR dose
#'   \item e response
#' }
#'
#' @format A data frame with 13 rows and 3 columns
#' @source \url{https://biostatistics.mdanderson.org/SoftwareDownload/}
#' @references Lee, J. J., & Kong, M. (2009). Confidence intervals 
#'  of interaction index for assessing multiple drug interaction. Statistics in biopharmaceutical
#'  research, 1(1), 4-17.
#' @name UMSCC22B
"UMSCC22B"

