

#' common code in fgui
#' 
#' @importFrom fgui guiGetValue guiSet
parse_alpha <- function(){
  alpha <- guiGetValue("alpha")
  if(alpha == 'Ignore') {
    alpha <- 1
  } else {
    alpha <- as.numeric(alpha)
  }	
  guiSet("alpha_", alpha) # alpha is modified; so save it as a global var
}


#' setup output dir
set_dirOutput_ <- function(){
  # set output dir according to input data
  dirOutput <- file.path(getwd(), 'drexplorer_results')
  if(!file.exists(dirOutput)){
    dir.create(dirOutput, recursive=TRUE)
  }
  guiSet("dirOutput_", dirOutput)
  return(dirOutput)
}


#' @importFrom fgui guiGetSafe
GUI_1_main <- function(datFile, exampleDat, models, freqModels, alpha, genDRfigure, outputRes){
  parse_alpha()
  dat <- guiGetSafe("PERSONAL_dataset")
  if(class(dat)[1] != "data.frame") stop("Data must be loaded.")
  models <- guiGetValue("models")
  alpha <- guiGetSafe("alpha_") # already passed in genDRfigure_press 
  resL <- fitOneExp(dat, drug = '', cellLine = '', tag = '', models = models,
                    alpha = alpha, plot = TRUE, transparency = 0.95)
  sprintf('%d outliers detected; The best model is: %s, corresponding IC50=%.3f',
          sum(resL$datWithOutlierStatus$isOutlier), 
          resL$ICx[1, 'Model'], resL$ICx['IC50'])
}

GUI_1_main_callback <- function(arg){
  if( arg == "datFile" ) {
    datFile_press()
  }else if( arg == "exampleDat" ) {
    exampleDat_press()
  }else if( arg == "freqModels" ) {
    freqModels_press()
  }else if( arg == "genDRfigure" ) {
    genDRfigure_press()
  }else if( arg == "outputRes" ) {
    outputRes_press()
  }
}


#' @importFrom utils read.csv read.delim
datFile_press <- function(){
  f_name <- guiGetValue("datFile")
  f_ext <- tolower(tools::file_ext(f_name))
  if(f_ext == 'csv'){
    data <- read.csv(f_name)
  } else {
    data <- read.delim(f_name)
  }
  guiSet("PERSONAL_dataset", data)
  tt <- set_dirOutput_()
}


#' @importFrom utils write.csv
#' @importFrom fgui guiSetValue
exampleDat_press <- function() {
  # set output dir according to input data
  dirOutput <- set_dirOutput_()
  dose <- ryegrass[, 2]
  response <- ryegrass[, 1]
  ryegrass_dat <- data.frame(dose = dose, response = response)
  ff <- file.path(dirOutput, 'ryegrass.csv')
  write.csv(ryegrass_dat, file = ff, row.names = F)
  guiSetValue("datFile", ff)
  GUI_1_main_callback("datFile")
}
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ryegrass"))


#' @importFrom fgui setListElements
freqModels_press <- function(){
  recModels <- c("sigEmax", "LL.4", "LL.5", "linear", "LL.3", "LL.3u", "logistic")
  setListElements("models", recModels)
}


genDRfigure_press <- function(){
  parse_alpha()
  resL <- fitOneExp(dat = guiGetSafe("PERSONAL_dataset"),
                    drug = '', cellLine = '', tag = '',
                    models = guiGetValue("models"), alpha = guiGetSafe("alpha_"),
                    plot = TRUE, transparency = 0.95)
}


#' @importFrom grDevices pdf dev.off
outputRes_press <- function(){
  parse_alpha()
  drug <- ''
  cellLine <- ''
  dirOutput <- set_dirOutput_()
  pdf(file = file.path(dirOutput, sprintf('DoseResponseCurve.pdf')))
  resL <- fitOneExp(dat = guiGetSafe("PERSONAL_dataset"),
                    drug = '', cellLine = '', tag = '',
                    models = guiGetValue("models"), alpha = guiGetSafe("alpha_"),
                    plot = TRUE, transparency = 0.95)
  dev.off()
  ICmat <- resL$ICmat
  datWithOutlierStatus <- resL$datWithOutlierStatus
  write.csv(ICmat, file = file.path(dirOutput, sprintf('ICtable.csv')))
  write.csv(datWithOutlierStatus, file = file.path(dirOutput, sprintf('datWithOutlierStatus.csv')))
}


#' Launch Graphical User Interface (GUI) for dose response curve fitting
#' 
#' This function will launch GUI for curve fitting. The GUI works across different platforms,
#'  but the appearance would be slightly different.
#'
#' @importFrom fgui gui
#' @export
drexplorerGUI_1 <- function(){
  gui(GUI_1_main,
      argFilename = list(datFile = NULL),
      argCommand = list(exampleDat = NULL, freqModels = NULL, genDRfigure = NULL, outputRes = NULL),
      argList = list(models = c(model.drc, model.DoseFinding)),
      argOption = list(alpha = c("0.05","0.01","Ignore")), 
      argGridOrder = c(1, 1, 2, 2, 3, 4, 5, 6),
      callback = GUI_1_main_callback,
      helps = list(datFile = "A tab-delimilated file or csv file with column 1 being the drug concentration and column 2 being the cell viability",
                   exampleDat = "Use example data from drexplorer package",
                   models = "A drop-down list of models to be fitted",
                   alpha = "Significance level in detecting outlier data points. Choose Ignore to omit outlier removal"),
      title = "drexplorer GUI-Curve Fitting",
      argText = list(datFile = "Load File ...",
                     exampleDat = "Use example data",
                     models = "Select models", 
                     freqModels = "Restrict to frequently used models",
                     alpha = "Significance level of outlier detection",
                     genDRfigure = "Generate dose-response curve",
                     outputRes = "Save results"),
      verbose = FALSE)
}




