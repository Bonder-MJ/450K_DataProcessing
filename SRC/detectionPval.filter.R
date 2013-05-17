# Performs a sample QC on the basis of detection p-values. Eliminate "bad quality" samples if their proportion (%) of probes associated to significant detection p-values < detectionPval.perc.threshold (80% by default).
#
# Args: 
#	- methLumi_data: methylumi object
#	- detectionPval.threshold: detection p-value threshold, set the threshold of significant signals (ok if p-value < detectionPval.threshold). Set to 0.01 by default.
#	- detectionPval.perc.threshold: threshold for the percentage of significant probes per sample. If a sample has less than this threshold it is coonsidered as bad quality. Set to 80% by defult.
#	- projectName: a string used to name the report file about sample quality. If NULL (by default) no report provided.
#	- PATH: a string setting the path for the folder where the report will be saved. Set to current folder by default.
#
# Returns: an updated methylumi object.



detectionPval.filter <- function(methLumi_data, detectionPval.threshold=0.01, detectionPval.perc.threshold=80, projectName = NULL, PATH="./"){

	#get detection p-values
	detectPval <- assayDataElement(methLumi_data, "detection")
	#get sample names
	samples <- colnames(detectPval)

	nbSignifPval <- vector()
	percentSignifPval <- vector()

	#for each sample compute the number and % or relevant signal (detection p-value < detectionPval.threshold)
	for(i in 1:length(samples)){
		nb <- length(which(detectPval[,i] <= detectionPval.threshold))
		percent <- nb*100/length(detectPval[,i])
		nbSignifPval <- c(nbSignifPval,nb)
		percentSignifPval <- c(percentSignifPval, percent)
	}

	rm(detectPval)
	#get "bad" samples indices
	index2remove <- which(percentSignifPval < detectionPval.perc.threshold)
	
	#remove "bad" samples from methylumi object
	if(length(index2remove)>0) methLumi_data <- methLumi_data[,-index2remove]
	
	index <- sort(percentSignifPval, index.return=TRUE)$ix
	
	# save sample quality report as text file
	if(!is.null(projectName)) write.table(list(samples=samples[index], nbSignifPval=nbSignifPval[index], percentSignifPval=percentSignifPval[index]), file=paste(PATH, projectName, "_signifPvalStats_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)

	return(methLumi_data)
}
