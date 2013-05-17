# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# Compute a set of reference quantiles from a data matrix. Variables are in rows and samples in columns so one quantile per row.
#
# Args: a data matrix with variables in rows and samples in columns.
#
# Returns: a vector of quantiles of size = length(rownames(data))
#
referenceQuantiles <- function(data){	
	
	sortedData <- apply(data, 2, sort, na.last=TRUE)
		
	Reference.Quantiles <- rowMeans(sortedData, na.rm=TRUE)
	
	indexNA <- which(is.na(Reference.Quantiles))
	if(length(indexNA)>0) Reference.Quantiles <- Reference.Quantiles[-indexNA]
	
	return(Reference.Quantiles)
}
