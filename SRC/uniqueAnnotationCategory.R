# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function parses a list of probe annotations and provide a list of preprocessed unique annotation values.
#
# Args:
#	- list: vector of strings, represent probe annotations.
#
# Returns: a vector of strings, representing unique and preprocessed probe annotations
#
uniqueAnnotationCategory <- function(list){
		temp <- function(x){
			res <- unique(as.vector(unlist(strsplit(x, split=";"))))
			if(length(res)>1) return("multipleAnnot")
			if(length(res)==1) return(res)
			if(length(res)==0) return(NA)
		}

	searchList <- apply(as.matrix(list), 1, temp)

	indexNA <- which(is.na(searchList))

	if(length(indexNA)>0) searchList[indexNA] <- "NA"

	searchList <- unique(searchList)
		
	return(searchList)
}
