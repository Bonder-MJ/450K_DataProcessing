# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function provides all probe IDs associated to a given annotation value.
#
# Args:
#	- annotations: a two columns array with probe IDs in the 1st column and a given annotation in the 2nd one
#	- annotationValue: a given annotation value for which we would like to identify the associated probes list
#	- uniqueAnnot: logical, if TRUE, provide probe IDs associated to a unique annotation value, else, provide all matching probe ID's.
#
# Returns: a vector of strings corresponding to probe IDs
#
findAnnotationProbes <- function(annotationValue, annotations, uniqueAnnot=TRUE){

	temp <- function(x){
		if(uniqueAnnot){
			res <- unique(as.vector(unlist(strsplit(x, split=";"))))
			if(length(res)>1) return("multipleAnnot")
			if(length(res)==1) return(res)
			if(length(res)==0) return(NA)
		}
		else{return(as.vector(unlist(strsplit(x, split=";")))[1])}
	}

	searchList <- apply(as.matrix(annotations[,2]), 1, temp)

	indexNA <- which(is.na(searchList))

	if(length(indexNA)>0) searchList[indexNA] <- "NA"

	if(length(searchList) != dim(annotations)[1]) {
		print("WARNING: in 'findAnnotationProbes.R' search list for annotation probes identification do not have same length than annotations")
	}

	indexProbes <- which(searchList == annotationValue)
	
	if(length(indexProbes)>0) probes <- annotations[indexProbes,1]
	else{print(paste("WARNING: no probes associated to annotation category:", annotationValue, sep="")) ; probes <- NULL}

	return(probes)
}
