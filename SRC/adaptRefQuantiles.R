# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function adapt a Reference Quantile vector to a larger or lower number of quantiles by respecting the initial distribution of quantile values.
#
# Args:
#	- Reference.Quantiles: vector of original reference quantiles.
#	- sizeNew: size of the new vector of quantiles
#
# Returns: a vector of quantiles of size = sizeNew
#
adaptRefQuantiles <- function(Reference.Quantiles, sizeNew, verbose=TRUE){
	#print(paste("Length sizeNew2: ", sizeNew, sep=""))

	M <- length(Reference.Quantiles)

	if(sizeNew == M){
		newReference.Quantiles <- Reference.Quantiles
		if(verbose) cat("\t\t adaptRefQuantiles:  sizeNew = length(Reference.Quantiles)")
	}
	else{ newReference.Quantiles <- getQuantiles(na.exclude(Reference.Quantiles), sizeNew)
		if(verbose){
			if(sizeNew < M) cat("\t\t adaptRefQuantiles:  sizeNew < length(Reference.Quantiles)")
			if(sizeNew > M) cat("\t\t adaptRefQuantiles:  sizeNew > length(Reference.Quantiles)")
		}
	}
	return(sort(newReference.Quantiles))
}
