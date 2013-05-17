# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# Adapts a Reference Quantile vector to a larger or lower number of quantiles by respecting the initial distribution of quantile values.
#
# Args:
#	- source: vector of original quantiles
#	- targetLength: size of the new vector of quantiles
#
# Returns: a vector of quantiles of size = targetLength
#
getQuantiles <- function(source, targetLength){

	if(targetLength < 2){
		print("WARNING! in getQuantiles.R, 'targetLenght' must be > 1")
		return("WARNING! in getQuantiles.R, 'targetLenght' must be > 1")
	}

	p <- 1/(targetLength-1)

	probs <- vector(, targetLength)
	probs[1]=0

	for(i in 2:targetLength) probs[i] <- probs[i-1]+p
	
	probs[targetLength] <- 1

	targetQuantiles <- quantile(x=sort(source), probs=probs, na.rm=TRUE)

	return(targetQuantiles)
}
