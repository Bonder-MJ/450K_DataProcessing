# Marc Jan Bonder
# m.j.bonder @ umcg.nl
#
############################
# Release data: 06-06-2014 #
############################

#
#######################
#### TO READ FIRST ####
#######################
# This script, when sourced ( source("Trans_meQTL_pipelineIlluminaMethylation.R") ), loads raw methylation data and performs a complete preprocessing and normalization of a batch of Illumina 450K data (corresponding to different plates).
# Please read and change carefully! Read up untill: "NOW YOU CAN SOURCE THIS SCRIPT !"
#
#################
# Pre-requisites:
# - install last 'lumi' and 'methylumi' bioconductor packages
# - all pacakges that are used: "lumi", "methylumi", "RPMM", "preprocessCore", "minfi", "matrixStats" and "IlluminaHumanMethylation450k.db"
# - data format (1 raw data):
#	- raw idat files
#
####################
# PIPELINE'S STEPS #
####################
#
# 1. Check that methylation data files and additionnal information files have the required format.
# 2. Open R environment (R console).
# 3. Increase the memory that R can use. In Windows, use the command 'memory.limit(20000)' for example.
# 4. set all pathes and pipeline options.
# 5. Source this script, in R environment, by using the command 'source("pathToThisScript/pipelineIlluminaMethylation.R").
# 6. Results are automatically saved, in 4 files (beta values, detection p-values and the same informations for SNP probes ("rs" probes that are removed befor data normalization)), as matrices.
# 7. Start the "real" exciting work !
#
#
###########
# OUTPUTS #
###########
#
# In the folder defined by the PATH_RES variable (see below) you will find:
#	- xxxx_signifPvalStats_threshold0.01.txt file : for each plate, per sample QC results (# and % of probes associated to a detection value lower than 0.01 (by default)
#	- xxx_beta.txt : normalized beta values for all plates (an m x n matrix where m rows represent probes and n columns represent samples)
# - xxx__beta_intermediate.txt : non-normalized beta values for all plates, values are filtered and background & color corrected (an m x n matrix where m rows represent probes and n columns represent samples)
#	- xxxx_detectionPvalue.txtv : normalized detection p-values for all plates (an m x n matrix where m rows represent probes and n columns represent samples)
#	- xxx_beta.RData : R compressed archive representing normalized beta values for all plates (an m x n matrix where m rows represent probes and n columns represent samples)
#	- xxxx_detectionPvalue.RData : R compressed archive representing normalized detection p-values for all plates (an m x n matrix where m rows represent probes and n columns represent samples)
#	- xxxx_betaSNPprobes.csv : beta values (non normalized) for SNP related probes, for all plates (an m' x n' matrix where m' rows represent probes and n' columns represent samples)
#	- xxxx_detectionPvalueSNPprobes.csv : detection p-values (non normalized) for SNP related probes, for all plates (an m' x n' matrix where m' rows represent probes and n' columns represent samples)
#	- QC plots :
#		. Hierarchical clustering plots : euclidean distance with Ward agglomerative method. 4 plots.
#		. Methylation profile plots: 4 plots.
#			# xxx_beta.raw.jpeg: one plot per plate. Raw beta values.
#			# xxx_beta.filter.jpeg: one plot per plate. Raw beta values after sample and probe filtering.
#			# xxx_beta.preproc.jpeg: one plot per plate. Corrected beta values, after color bias correction and background subtraction.
#			# xxx_beta.preproc.norm.jpeg: one plot for all plates. Normalized beta values
#
##############################
###### VARIABLES TO SET ######
##############################
#
### PATHs to files and folders
#set working directory (If necessary)
setwd("")
#
# If working on Windows set this:
#memory.limit(20000)
#
# set PATH to a folder of "projects" where each project corresponds to a folder of 450K plate extracted data. Only subfolders for plates can exist, otherwise the program will try to open any existing file as folder and crash.
PATH_PROJECT_DATA <- "./DATA/"
# set PATH to a folder were the results can be written to. 
PATH_RES <- "./RES/"
# set Project output name.
projectName = "ILLUMINA450K"

#
#
#

#################################################
#                   Defaults                    #
#################################################

	PATH_SRC <- "./SRC/"
	PATH_Annot <- "./ADDITIONAL_INFO/"
	PATH_RES <- "./RES/"
	PATH_ProbeSNP_LIST <- c(paste(PATH_Annot, "/ProbeFiltering/freq1percent/probeToFilter_450K_GoNL.hg19.ALL_alleleFreq1percent.txt", sep=""), paste(PATH_Annot, "/ProbeFiltering/ProbesBindingNonOptimal/Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt", sep=""))
	QCplot=FALSE
#
##############################################
#####  NOW YOU CAN SOURCE THIS SCRIPT !  #####
##############################################
#####
#####
##############################################
###### source scripts and load libraries #####
##############################################
require(lumi)
require(methylumi)
require(RPMM)
require(preprocessCore)
require(minfi)
require(matrixStats)
require(IlluminaHumanMethylation450k.db)

source(paste(PATH_SRC,"loadMethylumi2.R", sep=""))
source(paste(PATH_SRC,"lumiMethyR2.R", sep=""))
source(paste(PATH_SRC,"preprocessIlluminaMethylation.R", sep=""))
source(paste(PATH_SRC,"getMethylumiBeta.R", sep=""))
source(paste(PATH_SRC,"concatenateMatrices.R", sep=""))
source(paste(PATH_SRC,"normalizeIlluminaMethylation.R", sep=""))
source(paste(PATH_SRC,"nbBeadsFilter.R", sep=""))
source(paste(PATH_SRC,"detectionPval.filter.R", sep=""))
source(paste(PATH_SRC,"getSamples.R", sep=""))
source(paste(PATH_SRC,"filterXY.R", sep=""))
source(paste(PATH_SRC,"robustQuantileNorm_Illumina450K.R", sep=""))
source(paste(PATH_SRC,"coRankedMatrices.R", sep=""))
source(paste(PATH_SRC,"referenceQuantiles.R", sep=""))
source(paste(PATH_SRC,"adaptRefQuantiles.R", sep=""))
source(paste(PATH_SRC,"getQuantiles.R", sep=""))
source(paste(PATH_SRC,"normalize.quantiles2.R", sep=""))
source(paste(PATH_SRC,"robustQuantileNorm_Illumina450K.probeCategories.R", sep=""))
source(paste(PATH_SRC,"dataDetectPval2NA.R", sep=""))
source(paste(PATH_SRC,"uniqueAnnotationCategory.R", sep=""))
source(paste(PATH_SRC,"findAnnotationProbes.R", sep=""))
source(paste(PATH_SRC,"pipelineIlluminaMethylation.batch2.R", sep=""))
source(paste(PATH_SRC,"plotQC.R", sep=""))
source(paste(PATH_SRC,"plotMethylationDensity.R", sep=""))
source(paste(PATH_SRC,"hclustPlot.R", sep=""))
source(paste(PATH_SRC,"Additions/Type2_M-value_Correction.R", sep=""))
source(paste(PATH_SRC,"Additions/nasen.R", sep=""))
source(paste(PATH_SRC,"Average_U+M.filter.R", sep=""))
#
{
  if(is.character(PATH_ProbeSNP_LIST)){
    if(length(PATH_ProbeSNP_LIST)==1){
      probeSNP_LIST <- unlist(read.table(file=PATH_ProbeSNP_LIST, quote="", sep="\t", header=TRUE))  
    } else {
      probeSNP_LIST <- NULL
      for(id in 1:length(PATH_ProbeSNP_LIST)){
        probeSNP_LIST <- union(probeSNP_LIST, unlist(read.table(file=PATH_ProbeSNP_LIST[id], quote="", sep="\t", header=TRUE)))
      }
    }
  } 
  else{
    probeSNP_LIST <- NULL
  }
}
#
data.preprocess.norm <- NULL
print("DASEN normalization procedure")

  data.preprocess.norm <- pipelineIlluminaMethylation.batch2(
    PATH_PROJECT_DATA,
    PATH_Annot = PATH_Annot,
    projectName = projectName,
    qcAfterMerging = TRUE,
    nbBeads.threshold = NULL,
    detectionPval.threshold = 0.01,
    detectionPval.perc.threshold = NULL,
    detectionPval.perc.threshold2 = NULL,
    sampleSelection = FALSE,
    probeSNP_LIST = probeSNP_LIST,
    XY.filtering = "autosomal",
    colorBias.corr = FALSE,
    average.U.M.Check = FALSE,
    minimalAverageChanelValue = 0,
    maxratioDifference = 0,
    bg.adjust = "no",
    PATH = PATH_RES,
    QCplot = QCplot,
    betweenSampleCorrection = FALSE,
    includeQuantileNormOverChanel = FALSE,
    alfa = 100,
    NormProcedure = "DASEN",
    medianReplacement = FALSE,
    MvalueConv = TRUE
  )


if(is.null(data.preprocess.norm)){
  print("No samples selected")
  stop()
}

beta <- data.preprocess.norm$beta
detection.pvalue <- data.preprocess.norm$detection.pvalue



write.table(beta, file=paste(PATH_RES, projectName, "_Mval.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
write.table(detection.pvalue, file=paste(PATH_RES, projectName, "_detectionPvalue.txt", sep=""), sep="\t", col.names = NA)
save(beta, file=paste(PATH_RES, projectName, "_Mval.RData", sep=""))
save(detection.pvalue, file=paste(PATH_RES, projectName, "_detectionPvalue.RData", sep=""))



if(QCplot){
  plotQC(data.preprocess.norm$beta, figName=paste(projectName, "_beta.preproc.norm", sep=""), PATH = PATH_RES)
}
