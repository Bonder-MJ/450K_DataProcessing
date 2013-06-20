# Performs a sample QC on the basis of average U and M signlan. Eliminate "bad quality" samples if their Average U and M signal are out of normal distribution.
#
# Args: 
#	- methLumi_data: methylumi object
#	- PATH: a string setting the path for the folder where the report will be saved. Set to current folder by default.
#
# Returns: an updated methylumi object.



AverageUandM.filter <- function(methLumi_data, minimalAverageChanelValue, maxratioDifference, projectName = NULL, PATH="./"){

	#get detection methylation matrix
	unmethAveragePerSample <- colMeans(unmethylated(methLumi_data))
	methAveragePerSample <- colMeans(methylated(methLumi_data))
  rationMethUnmeth = methAveragePerSample / unmethAveragePerSample
  
	#get sample names
	samples <- names(unmethAveragePerSample)

  removalId1 <- which(unmethAveragePerSample < minimalAverageChanelValue)
	removalId2 <- which(methAveragePerSample < minimalAverageChanelValue)
	removalId3 <- union(which(rationMethUnmeth < 1/maxratioDifference), which(rationMethUnmeth > maxratioDifference))

	#get "bad" samples indices
	index2remove <- union(union(removalId1, removalId2), removalId3)
	rm(unmethAveragePerSample,methAveragePerSample,rationMethUnmeth,removalId1,removalId2,removalId3)
	
	#remove "bad" samples from methylumi object
	if(length(index2remove)>0) methLumi_data <- methLumi_data[,-index2remove]
  
	samples <- as.vector(colnames(unmethylated(methLumi_data)))
  
	unmethAveragePerSample <- colMeans(unmethylated(methLumi_data))
	methAveragePerSample <- colMeans(methylated(methLumi_data))
	rationMethUnmeth = methAveragePerSample / unmethAveragePerSample
  
	averageUSignal <- as.vector(unmethAveragePerSample)
	averageMSignal <- as.vector(methAveragePerSample)
	ratioUM <- as.vector(rationMethUnmeth)
  
	# save sample quality report as text file
	if(!is.null(projectName)) write.table(list(samples=samples, averageUSignal=averageUSignal, averageMSignal=averageMSignal, ratioUM=ratioUM), file=paste(PATH, projectName, "_U-M_Signal_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
	
	return(methLumi_data)
}
