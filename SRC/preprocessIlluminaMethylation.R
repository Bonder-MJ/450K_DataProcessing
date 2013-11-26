# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function performs a complete Illumina 450K array data preprocessing (but no normalization).
#
# Args:
#  - path2data: path to sample methylation '.txt' data file
#	- path2controlData: path to control probes methylation '.txt' data file
#	- projectName: project name to use for results
#	- nbBeads.threshold: number of minimal beads for significant probes (by default =3)
#	- detectionPval.threshold: detection p-value threshold for significancy (by default, significant detection p-values <0.01)
#	- detectionPval.perc.threshold: percentage of significant probe methylation signals in a given sample (by default, >80% for "good quality" samples)
#	- sample2keep: path to a sample ID list in case of selection of a subset of samples for preprocessing.
#	- probeSNP_LIST: list of probes associated to frequent SNP to filter out from dataset
#	- XY.filtering: logical, if 'TRUE', remove all "allosomal" probes (probes located on X and Y chromosomes).
#	- colorBias.corr: logical, if 'TRUE' performs a color bias correction of methylated and unmethylated signals.
#	- bg.adjust: character string, used to specify which kind of bg correction to perform (by default = "separatecolors")
#	- PATH: character string that specifies the location of results folder (by default ="./")
#
# Returns: methylumi object.
#

preprocessIlluminaMethylation <- function(
  qcAfterMerging = FALSE,
  path2data,
  path2controlData,
  projectName,
  nbBeads.threshold=3,
  detectionPval.threshold=0.01,
  detectionPval.perc.threshold=95,
  detectionPval.perc.threshold2=1,
  sample2keep,
  probeSNP_LIST,
  XY.filtering,
  colorBias.corr=TRUE,
  average.U.M.Check = FALSE,
  minimalAverageChanelValue = minimalAverageChanelValue,
  maxratioDifference = maxratioDifference,
  bg.adjust="separatecolors",
  PATH="./",
  QCplot=TRUE)
{
  #set pipeline steps counter
  i=0
  
  #loadMethylumi and store nbBeads_A and nbBeads_B info
  cat("\n\tStart data loading...\n")
  methLumi_data <- loadMethylumi2(methylationData = path2data, controlData = path2controlData)
  cat("\tProject sample nb: ", length(sampleNames(methLumi_data)), ".\n")
  cat("\tData dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
  cat("\t...data loaded..\n\n")
  
  #plot raw data QC
  if(QCplot){
    plotQC(getMethylumiBeta(methLumi_data), figName=paste(projectName, "_beta.raw", sep=""), PATH = PATH)
  }
  
  # remove controls or unrelevant samples
  if(!is.null(sample2keep)) {
    i<-i+1
    cat(" Step ", i, ": 'pertinent' sample selection...\n")
    sample2keep <- read.table(file=sample2keep, sep="\t", header=FALSE, quote="")[[1]]
    methLumi_data <- getSamples(methLumi_data, sample2keep)
    
    if(is.null(methLumi_data)){
      cat("\t Project samples nb after sample selection: 0.\n")
      cat("\t Skipping sub project.\n")
      return(NULL)
    }
    
    cat("\t Project samples nb after sample selection: ", length(sampleNames(methLumi_data)), ".\n")
    cat("\t...done.\n\n")
  }
  
  #Filter any non annotated probes
  probeFilter <- fData(methLumi_data)$TargetID[ which(fData(methLumi_data)$CHR=="-") ]
  indexFilter <- which(is.element(featureNames(methLumi_data), probeFilter))
  if(length(indexFilter) > 0) methLumi_data <- methLumi_data[-indexFilter,]
  rm(probeFilter, indexFilter)
  
  # XY chz filtering
  if(tolower(XY.filtering)=="autosomal"){
    i <- i+1
    cat(" Step ", i, ": start elimination of X & Y chr. probes...\n")
    methLumi_data <- filterXY(methLumi_data)
    cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
    cat("\t...done.\n\n")
  }else if(tolower(XY.filtering)=="allosomal"){
    i <- i+1
    cat(" Step ", i, ": start elimination of none X & Y chr. probes...\n")
    methLumi_data <- filterNoneXY(methLumi_data)
    cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
    cat("\t...done.\n\n")
  }
  
  # SNP filtering
  if(!is.null(probeSNP_LIST)){
    
    i <- i+1
    cat(" Step ", i, ": start frequent SNP filtering...\n")
    indexProbe2remove <- which(is.element(featureNames(methLumi_data), probeSNP_LIST))
    if(length(indexProbe2remove)>0) methLumi_data <- methLumi_data[-indexProbe2remove,]
    
    cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
    cat("\t...done.\n\n")
  }
  
  # nbBeads filtering
  if(!is.null(nbBeads.threshold)){
    i<-i+1
    cat(" Step ", i, ": start nb beads/probe filtering...\n")
    methLumi_data <- nbBeadsFilter(methLumi_data, nbBeads.threshold)
    cat("\t...done.\n\n")
  }
  
  # sample QC and filtering
  if((!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold)) || average.U.M.Check || (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold2))){
    i<-i+1
    
    # remove bad U + M samples
    cat(" Step ", i, ": start sample QC & filtering...\n")
    cat("\t Project samples nb. before QC: ", length(sampleNames(methLumi_data)), ", probes nb. before QC: ", dim(unmethylated(methLumi_data))[1],".\n")
    if(average.U.M.Check) {
      methLumi_data <- AverageUandM.filter(methLumi_data, minimalAverageChanelValue, maxratioDifference, projectName, PATH=PATH)
      cat("\t Project samples nb. average chanel filtering: ", length(sampleNames(methLumi_data)), ".\n")
    }
    
    if(length(sampleNames(methLumi_data))==0){
      cat("\t Warning: during sample QC all samples where removed.\n")
      return(NULL)
    }
    
    # remove bad p-value samples
    if(!qcAfterMerging && (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold))){
      methLumi_data <- detectionPval.filter(methLumi_data, detectionPval.threshold, detectionPval.perc.threshold, projectName, PATH=PATH)
      cat("\t Project samples nb. after after P-value filtering: ", length(sampleNames(methLumi_data)), ".\n")
    }
    
    if(length(sampleNames(methLumi_data))==0){
      cat("\t Warning: during sample QC all samples where removed.\n")
      return(NULL)
    }
    
    # remove bad probes
    if(!qcAfterMerging && (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold2))){
      methLumi_data <- detectionPval.filter2(methLumi_data, detectionPval.threshold, detectionPval.perc.threshold2, projectName, PATH=PATH)
      cat("\t Project probes nb. after after P-value filtering: ", dim(unmethylated(methLumi_data))[1], ".\n")
    }
    
    if(length(sampleNames(methLumi_data))==0){
      cat("\t Warning: during sample QC all probes where removed.\n")
      return(NULL)
    }
    cat("\t...done.\n\n")
  }
  
  #plot raw data QC
  if(QCplot){
    plotQC(getMethylumiBeta(methLumi_data), figName=paste(projectName, "_beta.filter", sep=""), PATH = PATH)
  }
  
  #Color bias correction
  if(colorBias.corr){
    i <- i+1
    cat(" Step ", i, ": start color bias correction...\n")
    methLumi_data <- adjColorBias.quantile(methLumi_data)
    cat("\t...done.\n\n")
  }
  
  #BG subtraction
  if(bg.adjust=="separatecolors"){
    i <- i+1
    cat(" Step ", i, ": start background subtraction (" , bg.adjust, ") ...\n")
    methLumi_data <- lumiMethyB(methLumi_data, separateColor = TRUE)
    cat("\t...done.\n\n")
  }
  if(bg.adjust=="unseparatecolors"){
    i <- i+1
    cat(" Step ", i, ": start background subtraction (" , bg.adjust, ") ...\n")
    methLumi_data <- lumiMethyB(methLumi_data, separateColor = FALSE, verbose = FALSE)
    cat("\t...done.\n\n")
  }
  if(bg.adjust=="no"){
    i <- i+1
    cat(" Step ", i, ": no background subtraction.\n")
  }
  
  #plot raw data QC
  if(QCplot){
    plotQC(getMethylumiBeta(methLumi_data), figName=paste(projectName, "_beta.preproc", sep=""), PATH = PATH)
  }
  
  return(methLumi_data)
}


# This function performs a complete Illumina 450K array data preprocessing (but no normalization), directly from .idat files.
#
# Args:
# - initial methylumiObject
# - SampleAnnotation matrix
#	- ControlAnnotation matrix
#	- projectName: project name to use for results
#	- nbBeads.threshold: number of minimal beads for significant probes (by default =3)
#	- detectionPval.threshold: detection p-value threshold for significancy (by default, significant detection p-values <0.01)
#	- detectionPval.perc.threshold: percentage of significant probe methylation signals in a given sample (by default, >80% for "good quality" samples)
#	- sample2keep: path to a sample ID list in case of selection of a subset of samples for preprocessing.
#	- probeSNP_LIST: list of probes associated to frequent SNP to filter out from dataset
#	- XY.filtering: logical, if 'TRUE', remove all "allosomal" probes (probes located on X and Y chromosomes).
#	- colorBias.corr: logical, if 'TRUE' performs a color bias correction of methylated and unmethylated signals.
#	- bg.adjust: character string, used to specify which kind of bg correction to perform (by default = "separatecolors")
#	- PATH: character string that specifies the location of results folder (by default ="./")
#
# Returns: methylumi object.
#

preprocessIlluminaMethylationIdat <- function(
  qcAfterMerging = FALSE,
  methLumi_dataTmpData,
  sampleAnnotationInfomation,
  projectName,
  nbBeads.threshold=3,
  detectionPval.threshold=0.01,
  detectionPval.perc.threshold=95,
  detectionPval.perc.threshold2=1,
  probeSNP_LIST,
  XY.filtering,
  colorBias.corr=TRUE,
  average.U.M.Check = FALSE,
  minimalAverageChanelValue = minimalAverageChanelValue,
  maxratioDifference = maxratioDifference,
  bg.adjust="separatecolors",
  PATH="./",
  QCplot=TRUE)
{
  featureData(methLumi_dataTmpData) <- sampleAnnotationInfomation
  #set pipeline steps counter
  i=0
  
  #plot raw data QC
  if(QCplot){
    plotQC(getMethylumiBeta(methLumi_dataTmpData), figName=paste(projectName, "_beta.raw", sep=""), PATH = PATH)
  }
  
  #Filter any non annotated probes
  probeFilter <- fData(methLumi_dataTmpData)$TargetID[ which(fData(methLumi_dataTmpData)$CHR=="-") ]
  indexFilter <- which(is.element(featureNames(methLumi_dataTmpData), probeFilter))
  if(length(indexFilter) > 0) methLumi_dataTmpData <- methLumi_dataTmpData[-indexFilter,]
  rm(probeFilter, indexFilter)
  
  # XY chz filtering
  if(tolower(XY.filtering)=="autosomal"){
    i <- i+1
    cat(" Step ", i, ": start elimination of X & Y chr. probes...\n")
    methLumi_dataTmpData <- filterXY(methLumi_dataTmpData)
    cat("\t Data dimensions: ", dim(methLumi_dataTmpData)[1],"x", dim(methLumi_dataTmpData)[2], ".\n")
    cat("\t...done.\n\n")
  }else if(tolower(XY.filtering)=="allosomal"){
    i <- i+1
    cat(" Step ", i, ": start elimination of none X & Y chr. probes...\n")
    methLumi_dataTmpData <- filterNoneXY(methLumi_dataTmpData)
    cat("\t Data dimensions: ", dim(methLumi_dataTmpData)[1],"x", dim(methLumi_dataTmpData)[2], ".\n")
    cat("\t...done.\n\n")
  }
  # SNP filtering
  if(!is.null(probeSNP_LIST)){
    
    i <- i+1
    cat(" Step ", i, ": start frequent SNP filtering...\n")
    indexProbe2remove <- which(is.element(featureNames(methLumi_dataTmpData), probeSNP_LIST))
    if(length(indexProbe2remove)>0) methLumi_dataTmpData <- methLumi_dataTmpData[-indexProbe2remove,]
    
    cat("\t Data dimensions: ", dim(methLumi_dataTmpData)[1],"x", dim(methLumi_dataTmpData)[2], ".\n")
    cat("\t...done.\n\n")
  }
  
  # nbBeads filtering
  if(!is.null(nbBeads.threshold)){
    i<-i+1
    cat(" Step ", i, ": start nb beads/probe filtering...\n")
    methLumi_dataTmpData <- nbBeadsFilterIdat(methLumi_dataTmpData, nbBeads.threshold)
    cat("\t...done.\n\n")
  }
  
  # sample QC and filtering
  if((!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold)) || average.U.M.Check || (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold2))){
    i<-i+1
    cat(" Step ", i, ": start sample QC & filtering...\n")
    cat("\t Project samples nb. before QC: ", length(sampleNames(methLumi_dataTmpData)), ", probes nb. before QC: ", dim(unmethylated(methLumi_dataTmpData))[1],".\n")
    
    # remove bad U + M samples
    if(average.U.M.Check) {
      methLumi_dataTmpData <- AverageUandM.filter(methLumi_dataTmpData, minimalAverageChanelValue, maxratioDifference, projectName, PATH=PATH)
      cat("\t Project samples nb. average chanel filtering: ", length(sampleNames(methLumi_dataTmpData)), ".\n")
    }
    
    if(length(sampleNames(methLumi_dataTmpData))==0){
      cat("\t Warning: during sample QC all samples where removed.\n")
      return(NULL)
    }
    
    if(!qcAfterMerging && (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold))){
      methLumi_dataTmpData <- detectionPval.filter(methLumi_dataTmpData, detectionPval.threshold, detectionPval.perc.threshold, projectName, PATH=PATH)
      cat("\t Project samples nb. after after P-value filtering: ", length(sampleNames(methLumi_dataTmpData)), ".\n")
    }
    
    if(length(sampleNames(methLumi_dataTmpData))==0){
      cat("\t Warning: during sample QC all samples where removed.\n")
      return(NULL)
    }
    
    if(!qcAfterMerging && (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold2))){
      methLumi_dataTmpData <- detectionPval.filter2(methLumi_dataTmpData, detectionPval.threshold, detectionPval.perc.threshold2, projectName, PATH=PATH)
      cat("\t Project probes nb. after after P-value filtering: ", dim(unmethylated(methLumi_dataTmpData))[1], ".\n")
    }
    
    if(length(sampleNames(methLumi_dataTmpData))==0){
      cat("\t Warning: during sample QC all probes where removed.\n")
      return(NULL)
    }
    cat("\t...done.\n\n")
    
  }
  
  #plot raw data QC
  if(QCplot){
    plotQC(getMethylumiBeta(methLumi_dataTmpData), figName=paste(projectName, "_beta.filter", sep=""), PATH = PATH)
  }
  
  #Color bias correction
  if(colorBias.corr){
    i <- i+1
    cat(" Step ", i, ": start color bias correction...\n")
    methLumi_dataTmpData <- adjColorBias.quantile(methLumi_dataTmpData)
    cat("\t...done.\n\n")
  }
  
  #BG subtraction
  if(bg.adjust=="separatecolors"){
    i <- i+1
    cat(" Step ", i, ": start background subtraction (" , bg.adjust, ") ...\n")
    methLumi_dataTmpData <- lumiMethyB(methLumi_dataTmpData, separateColor = TRUE)
    cat("\t...done.\n\n")
  }
  if(bg.adjust=="unseparatecolors"){
    i <- i+1
    cat(" Step ", i, ": start background subtraction (" , bg.adjust, ") ...\n")
    methLumi_dataTmpData <- lumiMethyB(methLumi_dataTmpData, separateColor = FALSE, verbose = FALSE)
    cat("\t...done.\n\n")
  }
  if(bg.adjust=="no"){
    i <- i+1
    cat(" Step ", i, ": no background subtraction.\n")
  }
  
  #plot raw data QC
  if(QCplot){
    plotQC(getMethylumiBeta(methLumi_dataTmpData), figName=paste(projectName, "_beta.preproc", sep=""), PATH = PATH)
  }
  
  return(methLumi_dataTmpData)
}


qcAfterMerg <- function(
  matrix,
  detectionP,
  detectionPval.threshold=0.01,
  detectionPval.perc.threshold=95,
  detectionPval.perc.threshold2=1,
  PATH = PATH
  ){
    if((!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold)) || (!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold2))){
      cat("\n")
      cat(" Start after merging sample QC & filtering...\n")
      cat("\t Project samples nb. before QC: ", ncol(matrix), ", probes nb. before QC: ", nrow(matrix),".\n")
      
      if(!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold)){

        t <- detectionPval.filter_2(matrix, detectionP, detectionPval.threshold, detectionPval.perc.threshold, projectName, PATH=PATH)
        matrix <- t[[1]]
        detectionP <- t[[2]]
        cat("\t Project samples nb. after after P-value filtering: ", ncol(matrix), ".\n")
      }
         
       if(ncol(matrix)==0){
         cat("\t Warning: during sample QC all samples where removed.\n")
         return(NULL)
       }
       
       if(!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold2)){
         t <- detectionPval.filter2_2(matrix, detectionP, detectionPval.threshold, detectionPval.perc.threshold2, projectName, PATH=PATH)
         matrix <- t[[1]]
         detectionP <- t[[2]]
         cat("\t Project probes nb. after after P-value filtering: ", dim(matrix)[1], ".\n")
       }
       
       if(ncol(matrix)==0 && colnames(matrix)){
         cat("\t Warning: during sample QC all probes where removed.\n")
         return(NULL)
       }
       cat("\t...done.\n\n")
         
    }
    return(list(matrix,detectionP))
}