# Filter non significant probe signals according to probe associated bead number.
# The filtering is based on probe detection p-value updating.
# Non significant probes (bead nb < threshold) are associated to a detection p-value equal to 1.
# 
# Args:
#	- methyLumi: a methylumi object (see lumi & methylumi packages)
#	- nbBeads.threshold: minimal bead number per probe for probe selection. By default set to 3 following Illumina company recommandation.
#
# Returns: a methylumi object.

nbBeadsFilter <- function(methyLumi, nbBeads.threshold = 3){

	#get nBeads for A and B signals (unmethylated and methylated signals)
	nbBeadA <- assayDataElement(methyLumi, "Avg_NBEADS_A")
	nbBeadB <- assayDataElement(methyLumi, "Avg_NBEADS_B")
	
	#set nBeadsA and B to NULL
	assayDataElement(methyLumi, "Avg_NBEADS_A") <- NULL
	assayDataElement(methyLumi, "Avg_NBEADS_B") <- NULL
	
	#identify nbBeads < 3
	indexA <- which(nbBeadA < 3)
	indexB <- which(nbBeadB < 3)
	indexAB <- union(indexA, indexB)
	rm(indexA, indexB)

	#get detection pvalues
	detectPval <- assayDataElement(methyLumi, "detection")

	#set detection p-values associated to non significant signals to 1
	detectPval[indexAB] <- 1
	assayDataElement(methyLumi, "detection") <- detectPval
	rm(detectPval, indexAB)

	return(methyLumi)
}
