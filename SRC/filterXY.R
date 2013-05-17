# Removes probes associated to allosomes (chr X and Y)
# 
# Args: methylumi object
#
# Returns: methylumi object

filterXY <- function(methLumi_data){

	probeXY <- fData(methLumi_data)$TargetID[ which(fData(methLumi_data)$CHR=="Y" | fData(methLumi_data)$CHR=="X") ]
	indexXY <- which(is.element(featureNames(methLumi_data), probeXY))
	rm(probeXY)
	if(length(indexXY) > 0) methLumi_data <- methLumi_data[-indexXY,]
	
	print(paste("Nb probes on X & Y chr to remove: ", length(indexXY), sep=""))
	
	return(methLumi_data)
}
