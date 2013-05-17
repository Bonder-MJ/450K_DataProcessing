# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function performs a complete Illumina 450K array data preprocessing (but no normalization).
#
# Args:
#	- path2data: path to sample methylation '.txt' data file
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
	path2data,
	path2controlData,
	projectName,
	nbBeads.threshold=3,
	detectionPval.threshold=0.01,
	detectionPval.perc.threshold=80,
	sample2keep,
	probeSNP_LIST,
	XY.filtering,
	colorBias.corr=TRUE,
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

# nbBeads filtering
	if(!is.null(nbBeads.threshold)){
		i<-i+1
		cat(" Step ", i, ": start nb beads/probe filtering...\n")
		methLumi_data <- nbBeadsFilter(methLumi_data, nbBeads.threshold)
		cat("\t...done.\n\n")
	}
	
# sample QC and filtering
	if(!is.null(detectionPval.threshold) && !is.null(detectionPval.perc.threshold)){
		i<-i+1
		cat(" Step ", i, ": start sample QC & filtering...\n")
		methLumi_data <- detectionPval.filter(methLumi_data, detectionPval.threshold, detectionPval.perc.threshold, projectName, PATH=PATH)
		cat("\t Project samples nb. after sample QC: ", length(sampleNames(methLumi_data)), ".\n")
		cat("\t...done.\n\n")
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
# SNP filtering
	if(!is.null(probeSNP_LIST)){
		
		i <- i+1
		cat(" Step ", i, ": start frequent SNP filtering...\n")
		indexProbe2remove <- which(is.element(featureNames(methLumi_data), probeSNP_LIST))
		if(length(indexProbe2remove)>0) methLumi_data <- methLumi_data[-indexProbe2remove,]

		cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
		cat("\t...done.\n\n")
	}
	
# XY chz filtering
	if(XY.filtering){
		i <- i+1
		cat(" Step ", i, ": start elimination of X & Y chr. probes...\n")
		methLumi_data <- filterXY(methLumi_data)
		cat("\t Data dimensions: ", dim(methLumi_data)[1],"x", dim(methLumi_data)[2], ".\n")
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
