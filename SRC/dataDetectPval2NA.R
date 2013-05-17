# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# Replace with 'NA' beta values entries associated to detection p-values < threshold.
#
# Args:
#	- data: data matrix of beta values, with variables in rows and samples in columns
#	- detect.pval: data matrix of detection p-values, with variables in rows and samples in columns
#	- threshold: detection p-value threshold
#
# Returns: a vector of quantiles of size = sizeNew
#
dataDetectPval2NA <- function(data, detect.pval, threshold){

	indexNA <- which(detect.pval > threshold, arr.ind=TRUE)

	data[indexNA] <- NA
					
	return(data)
}

