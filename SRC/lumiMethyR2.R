# Load raw methylation data (extracted with Genome Studio) and provide a methylumi object (see methylumi and lumi packkages).
# This function corresponds to a modified version of the 'lumiMethyR' provided by lumi package in order to get bead number information for each signal (methylated & unmethylated.
#
# Args: please refer to the help for 'lumiMethyR' function provided by lumi package
#
# Returns: a methylumi object.
#

lumiMethyR2 <- function (data, lib = NULL, controlData = NULL){

    methyLumiSet <- methylumiR(data, sep="\t")

	#sort by sample names
	indexSamples <- sort(sampleNames(methyLumiSet), index.return=TRUE)$ix
	methyLumiSet <- methyLumiSet[,indexSamples]
	
	methyLumiM <- as(methyLumiSet, "MethyLumiM")
	# 'A' for unmethylated and 'B' for methylated
	assayDataElement(methyLumiM, "Avg_NBEADS_A") <- assayDataElement(methyLumiSet, "Avg_NBEADS_A")
	assayDataElement(methyLumiM, "Avg_NBEADS_B") <- assayDataElement(methyLumiSet, "Avg_NBEADS_B")
	rm(methyLumiSet)
	
    if (!is.null(lib)) methyLumiM <- addAnnotationInfo(methyLumiM, lib = lib)
    if (!is.null(controlData)){
        if (is.character(controlData)) controlData <- methylumiR(controlData)
        if (is(controlData, "MethyLumiQC")) controlData(methyLumiM) <- controlData
        else { cat("Provided controlData is not supported!\n") }
    }
	
    return(methyLumiM)
}
