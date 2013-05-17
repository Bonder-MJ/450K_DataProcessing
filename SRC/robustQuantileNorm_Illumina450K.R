# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function launches the subset based quantile normalization for annotation based probe categories.
#
# Args: 
#	- data: beta values matrix with columns representing samples and rows representing probes.
#	- infiniumI: vector of characters representing InfiniumI probe IDs
#	- infiniumII: vector of characters representing InfiniumII probe IDs
#	- detect.pval: detection p-values matrix with columns representing samples and rows representing probes.
#	- detect.pval.threshold: numeric, represents the threshold for significant detection p-values (significant if < detect.pval.threshold). Set to 0.01 by default.
#	- annotations: Illumina probe annotations.
#
# Returns: a list of two matrices, the first one corresponding to normalized beta-values and the second one to updated detection p-values.
#
robustQuantileNorm_Illumina450K <- function(
	data,
	infiniumI,
	infiniumII,
	detect.pval,
	detect.pval.threshold,
	annotations,
	verbose = TRUE){
	
	sampleNames <- colnames(data)
		
	indexInfiniumI <- which(is.element(rownames(data), infiniumI))
	indexInfiniumII <- which(is.element(rownames(data), infiniumII))
	data.infiniumI <- data[indexInfiniumI,]
	data.infiniumII <- data[indexInfiniumII,]
	probeID.infiniumI <- rownames(data.infiniumI)
	probeID.infiniumII <- rownames(data.infiniumII)
	rm(data)
				
	#print(paste("Separate, robust and filtered quantile normalization per probe categories: reference quantiles for Infinium I and II probes are based on Infinium I signals associated to p-values < ", detect.pval.threshold ,". Reference quantiles are computed separately for each kind of probe annotation.", sep=""))

	data.norm <- robustQuantileNorm_Illumina450K.probeCategories(data.infiniumI, data.infiniumII, detect.pval.infiniumI = detect.pval[indexInfiniumI,], annotations = annotations, threshold = detect.pval.threshold)

	rm(data.infiniumI, data.infiniumII, indexInfiniumI, indexInfiniumII)
		
	if(verbose) {
		cat("\tNormalized beta values example:\n")
		cat("\t", data.norm[1:2,1:2],"\n")
	}
	
	t <- coRankedMatrices(data.norm, detect.pval)
	#print(str(t))
	return(t)
}
