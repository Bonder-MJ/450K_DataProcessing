# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
######################
#### TO READ FIRST ####
######################
# setwd("d:/UMCG/450K_pipeline_Nov2012release/")
# This script, when sourced ( source("pipelineIlluminaMethylation.main.R") ), loads raw methylation data and performs a complete preprocessing and normalization of a batch of Illumina 450K data (corresponding to different plates).
#
#################
# Pre-requisites:
# - install last 'lumi' and 'methylumi' bioconductor packages
# - data format:
#	- raw methylation data: methylation data have to be extracted with GenomeStudio in two text files, one corresponding to sample methylation informations and the second one to control probe informations.
#		- sample methylation informations: must contain the following columns in addition to the columns pre-selected by GenomeStudio(before exporting to '.txt' files the extracted informations).
#  		. Index
#			. TargetID
#			. INFINIUM_DESIGN_TYPE
#			. Color_Channel
#			. CHR
#			. Mapinfo (probe position for hg19)
#			. COORDINATE_36 (probe position for hg18)
#			. STRAND
#			. Probe_SNPS
#			. Probe_SNPS_10
#			. UCSC_Refgene_Name
#			. UCSC_Refgene_Accession
#			. UCSC_Refgene_Group
#			. Relation_To_UCSC_CPG_Island
#			. DetectionPval
#			. Signal_A
#			. Signal_B
#			. Avg_NBEADS_A
#			. Avg_NBEADS_B
#		- control probes methylation informations: all columns
#		- additional data:
#			- sample IDs list in case of only a subset of the loaded samples have to be processed and normalized (to remove control samples for example): simple text file, with no header and with one sample ID per line.
#			- probe IDs list for probe filtering (probe associated to frequent SNP for example): simple text file with one probe ID per line.
#				We provide 36 lists of probes associated to frequent SNP obtained from the 1000 genome project (these 36 lists correspond to 18 different human populations. See ADDITIONAL_INFO folder) :
# 				- list of probes associated to SNP with alternative allele frequency >= 5% (the most stringent list)
# 				- list of probes associated to SNP with alternative allele frequency >= 10%
#			We consider that a probe is associated to a frequent SNP if a frequent SNP is located inside the probe sequence or if it is located in +1 (for Infinium I probes) or +2 (for Infinium II probes) position from probe's query site.
#
####################
# PIPELINE'S STEPS #
####################
#
# 0. Extract methylation data with annotation informations with GenomeStudio software.
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
#
# set PATH to R pipeline's scripts for script sourcing
setwd("d:/UMCG/450K_pipeline_Nov2012release/")
memory.limit(20000)
PATH_SRC <- "./SRC/"
#
# set PATH to results folder
#PATH_RES <- "./RES_Muscle_SWan_NoQN/"
#PATH_RES <- "./KI_RES_BMIQ/"
#PATH_RES <- "./RES_KI_SQN/"
PATH_RES <- "./RES_test/"
#PATH_RES <- "./KI_RES_SWAN/Comp3"

#
# set PATH to a folder of "projects" where each project corresponds to a folder of 450K plate extracted data : control probes methylation informations, raw sample methylation informations (data extracted with GenomeStudio) and eventually a sample IDs list to select. Only subfolders for plates can exist, otherwise the program will try to open any existing file as folder and crash.
# requirement for data file naming :
#	- control probes file: file name must starts with the pattern "TableControl"
#	- raw sample methylation file: file name must starts with the pattern "TableSample"
#	- sample IDs file for sample selection (not compulsory): file name must contain the pattern "sampleList"

#PATH_PROJECT_DATA <- "./DATA_KI/"
#PATH_PROJECT_DATA <- "./DATA/"
PATH_PROJECT_DATA <- "./DATA_test/"

#
# set PATH to the file with frequent SNP informations, on which SNP filtering is based. If = NULL, no probe removed.
#PATH_ProbeSNP_LIST <- "./ADDITIONAL_INFO/ProbeFiltering/freq5percent/probeToFilter_450K_1000G_omni2.5.hg19.EUR_alleleFreq5percent_50bp_wInterroSite.txt"
PATH_ProbeSNP_LIST <- NULL
#
#######################################
### set pipeline options and parameters
#
# The name that will be given to result files
projectName = "ILLUMINA450K"
#
# Minimal bead number for considering that a probe worked. If =NULL, does not perform bead nb. based filtering.
nbBeads.threshold = 3
#
# Threshold for detection p-values significance. A signal is considered as significant if its associated detection p-values < this threshold (0.01 by default).
detectionPval.threshold = 0.01
#
# Percentage of significant probe methylation signals in a given sample (by default, >80% for "good quality" samples). This is used for samples QC and filtering. All samples that do not respect this condition will be removed from the analysis.
detectionPval.perc.threshold = 97
#
# If sampleSelection= FALSE , all loaded samples will be processed and normalized, if sampleSelection = TRUE, a sample IDs text list, with no header and with the pattern "sampleList" in file name, will be loaded and used to select the samples to preprocess and normalize.
sampleSelection = FALSE
#
# if 'TRUE', remove all "allosomal" probes (probes located on X and Y chromosomes).
XY.filtering = TRUE
#
# if 'TRUE' performs a color bias correction of methylated and unmethylated signals.
colorBias.corr = TRUE
#
# if 'TRUE' a additional sample filtering based on average M and U values is preformed. Cut-offs are designed based on a set of >2000 450k samples. 
average.U.M.Check = TRUE
minimalAverageChanelValue = 3108.038
maxratioDifference = 0.2506493
#
# If "separatecolors", performs a separate color bg adjustement (recommended), if "unseparatecolors" performs a non separate color bg adjustement , if "no" don't perform any bg adjustement.
bg.adjust = "separatecolors"
#
# Specifying which kind of probe annotation to use for probe categories construction, must be one of "relationToCpG" or "relationToSequence". We recommand the more robust option: probeAnnotationsCategory = "relationToCpG".
probeAnnotationsCategory = "relationToCpG"
#
# If QCplot==TRUE, performs and plots hierarchical clustering and plots methylation density curves.
QCplot=FALSE
#
# Select the normalization procedure. Switch between SQN / SWAN / M-Value correction  and BMIQ. If not give / write error this will default to SQN
NormProcedure = "SQN"
#NormProcedure = "SWAN"
#NormProcedure = "BMIQ"
#NormProcedure = "M-ValCor"

#Set alfa, used during transformation from U + M to Beta value
alfa = 100

# Do QN after SWAN / BMIQ / M-ValCor"
betweenSampleCorrection = FALSE

# Do M-value conversion.
MvalueConv = FALSE

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
require(compiler)

enableJIT(1)

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
source(paste(PATH_SRC,"pipelineIlluminaMethylation.batch.R", sep=""))
source(paste(PATH_SRC,"plotQC.R", sep=""))
source(paste(PATH_SRC,"plotMethylationDensity.R", sep=""))
source(paste(PATH_SRC,"hclustPlot.R", sep=""))
source(paste(PATH_SRC,"Additions/BMIQ_1.1_Pipeline.R", sep=""))
source(paste(PATH_SRC,"Additions/swan2.R", sep=""))
source(paste(PATH_SRC,"Additions/Type2_M-value_Correction.R", sep=""))
source(paste(PATH_SRC,"Average_U+M.filter.R", sep=""))
#
#
#
	{
	if(is.character(PATH_ProbeSNP_LIST)) probeSNP_LIST <- unlist(read.table(file=PATH_ProbeSNP_LIST, quote="", sep="\t", header=TRUE))
	else{probeSNP_LIST <- NULL}
	}
#
  data.preprocess.norm <- NULL
  print(paste(NormProcedure ,"normalization procedure"))
  if(NormProcedure != "SWAN"){
    data.preprocess.norm <- pipelineIlluminaMethylation.batch(
      PATH_PROJECT_DATA,
      projectName = projectName,
      nbBeads.threshold = nbBeads.threshold,
      detectionPval.threshold = detectionPval.threshold,
      detectionPval.perc.threshold = detectionPval.perc.threshold,
      sampleSelection = sampleSelection,
      probeSNP_LIST = probeSNP_LIST,
      XY.filtering = XY.filtering,
      colorBias.corr = colorBias.corr,
      average.U.M.Check = average.U.M.Check,
      minimalAverageChanelValue = minimalAverageChanelValue,
      maxratioDifference = maxratioDifference,
      bg.adjust = bg.adjust,
      PATH = PATH_RES,
      QCplot = QCplot,
      betweenSampleCorrection = betweenSampleCorrection,
      alfa,
	    NormProcedure
    )
  } else {
    data.preprocess.norm <- pipelineIlluminaMethylation.batch.SWAN(
      PATH_PROJECT_DATA,
      projectName = projectName,
      nbBeads.threshold = nbBeads.threshold,
      detectionPval.threshold = detectionPval.threshold,
      detectionPval.perc.threshold = detectionPval.perc.threshold,
      sampleSelection = sampleSelection,
      probeSNP_LIST = probeSNP_LIST,
      XY.filtering = XY.filtering,
      colorBias.corr = colorBias.corr,
      average.U.M.Check = average.U.M.Check,
      minimalAverageChanelValue = minimalAverageChanelValue,
      maxratioDifference = maxratioDifference,
      bg.adjust = bg.adjust,
      PATH = PATH_RES,
      QCplot = QCplot,
      alfa,
      betweenSampleCorrection = betweenSampleCorrection
    )
  }
  
  if(is.null(data.preprocess.norm)){
    print("No samples selected")
    return(0)
  }

	beta <- data.preprocess.norm$beta
	detection.pvalue <- data.preprocess.norm$detection.pvalue
	
  if(MvalueConv){
    
    for(i in 1:ncol(beta)){
      beta1 <- beta[,i]
      m1 <- log2(beta1/(1 - beta1))  
      beta[,i] <- m1
    }
    
    
    write.table(beta, file=paste(PATH_RES, projectName, "_Mval.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
    #write.table(detection.pvalue, file=paste(PATH_RES, projectName, "_detectionPvalue.txt", sep=""), sep="\t", col.names = NA)
    #save(beta, file=paste(PATH_RES, projectName, "_Mval.RData", sep=""))
    #save(detection.pvalue, file=paste(PATH_RES, projectName, "_detectionPvalue.RData", sep=""))
    
  } else {
    write.table(beta, file=paste(PATH_RES, projectName, "_beta.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
    write.table(detection.pvalue, file=paste(PATH_RES, projectName, "_detectionPvalue.txt", sep=""), sep="\t", col.names = NA)
    #save(beta, file=paste(PATH_RES, projectName, "_beta.RData", sep=""))
    #save(detection.pvalue, file=paste(PATH_RES, projectName, "_detectionPvalue.RData", sep=""))
    
  }
	
	if(QCplot){
		plotQC(data.preprocess.norm$beta, figName=paste(projectName, "_beta.preproc.norm", sep=""), PATH = PATH_RES)
	}
