# Build a sub-methylumi object from a methylumi object corresponding to a subset of samples.
#
# Args:
#	- methLumi_data: methylumi object (see lumi & methylumi package)
#	- sample2keep: vector of sample IDs to select
#
# Returns an updated methylumi object limitted to the sample list provided in sample2keep.

getSamples <- function(methLumi_data, sample2keep){

	samples <- sampleNames(methLumi_data)
	index <- which(is.element(samples, sample2keep))

	if(length(index)>0){
		methLumi_data <- methLumi_data[, index]
	}
	else{
		print("WARNING ! No samples to select !")
		methLumi_data<-NULL
	}
	
	return(methLumi_data)
}
