# Load raw methylation data (extracted with Genome Studio) and provide a methylumi object (see methylumi and lumi packkages).
#
# Args:
#	methylationData: string providing the path to ".txt" file with methylation signals and probe annotation informations extracted from GenomeStudio
#	controlData: string providing the path to ".txt" file with control probe signals and annotations extracted from GenomeStudio
#
# Returns: a methylumi object.
#
loadMethylumi2 <- function(methylationData, controlData=NULL){

	methLumi_data <- lumiMethyR2(methylationData)
	
	if(!is.null(controlData)) methLumi_data <- addControlData2methyLumiM(controlData, methLumi_data)

	return(methLumi_data)
}
