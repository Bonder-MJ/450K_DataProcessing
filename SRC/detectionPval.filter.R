# Performs a sample QC on the basis of detection p-values. Eliminate "bad quality" samples if their proportion (%) of probes associated to significant detection p-values < detectionPval.perc.threshold (95% by default).
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
	nbSignifPval <- as.vector(colSums(detectPval <= detectionPval.threshold))
	percentSignifPval <- as.vector((colSums(detectPval <= detectionPval.threshold)/nrow(methLumi_data))*100)
  
	#for each sample compute the number and % or relevant signal (detection p-value < detectionPval.threshold)
	nrMin <- (nrow(methLumi_data)/100)*detectionPval.perc.threshold
	#get "bad" samples indices
	index2remove <- which(colSums(detectPval <= detectionPval.threshold) < nrMin)
  
  rm(detectPval)
  
	if(!is.null(projectName)) write.table(list(samples=samples[index2remove], nbSignifPval=nbSignifPval[index2remove], percentSignifPval=percentSignifPval[index2remove]), file=paste(PATH, projectName, "_non-signifPvalStats-Sample_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
	#remove "bad" samples from methylumi object
	if(length(index2remove)>0) methLumi_data <- methLumi_data[,-index2remove]
	
  
	detectPval <- assayDataElement(methLumi_data, "detection")
	samples <- colnames(detectPval)
	index <- order(colSums(detectPval))
	
	nbSignifPval <- as.vector(colSums(detectPval <= detectionPval.threshold))
	percentSignifPval <- as.vector((colSums(detectPval <= detectionPval.threshold)/nrow(methLumi_data))*100)
  
	# save sample quality report as text file
	if(!is.null(projectName)) write.table(list(samples=samples[index], nbSignifPval=nbSignifPval[index], percentSignifPval=percentSignifPval[index]), file=paste(PATH, projectName, "_signifPvalStats-Sample_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)

	return(methLumi_data)
}

# Performs a probe QC on the basis of detection p-values. Eliminate "bad quality" probes if their proportion (%) of probes associated to significant detection p-values < detectionPval.perc.threshold (1% by default).
#
# Args: 
#  - methLumi_data: methylumi object
#	- detectionPval.threshold: detection p-value threshold, set the threshold of significant signals (ok if p-value < detectionPval.threshold). Set to 0.01 by default.
#	- detectionPval.perc.threshold: threshold for the percentage of significant probes per sample. If a sample has less than this threshold it is coonsidered as bad quality. Set to 80% by defult.
#	- projectName: a string used to name the report file about sample quality. If NULL (by default) no report provided.
#	- PATH: a string setting the path for the folder where the report will be saved. Set to current folder by default.
#
# Returns: an updated methylumi object.


detectionPval.filter2 <- function(methLumi_data, detectionPval.threshold=0.01, detectionPval.perc.threshold=1, projectName = NULL, PATH="./"){
  
  #get detection p-values
  detectPval <- assayDataElement(methLumi_data, "detection")
  #get sample names
  probeNames <- rownames(detectPval)
  nbSignifPval <- as.vector(rowSums(detectPval <= detectionPval.threshold))
  percentSignifPval <- as.vector((rowSums(detectPval <= detectionPval.threshold)/ncol(methLumi_data))*100)
    
  #for each probe compute the number and % or relevant signal (detection p-value < detectionPval.threshold)
  
  nrMax <- (ncol(methLumi_data)/100)*detectionPval.perc.threshold
  #get "bad" samples indices
  index2remove <- which(rowSums(detectPval > detectionPval.threshold) > nrMax)
  
  rm(detectPval)
  
  if(!is.null(projectName)) write.table(list(probeNames=probeNames[index2remove], nbSignifPval=nbSignifPval[index2remove], percentSignifPval=percentSignifPval[index2remove]), file=paste(PATH, projectName, "_non-signifPvalStats-Probe_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
  #remove "bad" samples from methylumi object
  if(length(index2remove)>0) methLumi_data <- methLumi_data[-index2remove,]
  
  detectPval <- assayDataElement(methLumi_data, "detection")
  probeNames <- rownames(detectPval)
  index <- order(rowSums(detectPval), decreasing=T)
  
  nbSignifPval <- as.vector(rowSums(detectPval <= detectionPval.threshold))
  percentSignifPval <- as.vector((rowSums(detectPval <= detectionPval.threshold)/ncol(methLumi_data))*100)
  
  # save sample quality report as text file
  if(!is.null(projectName)) write.table(list(probeNames=probeNames[index], nbSignifPval=nbSignifPval[index], percentSignifPval=percentSignifPval[index]), file=paste(PATH, projectName, "_signifPvalStats-Probe_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
  return(methLumi_data)
}


# Performs a sample QC on the basis of detection p-values. Eliminate "bad quality" samples if their proportion (%) of probes associated to significant detection p-values < detectionPval.perc.threshold (95% by default).
#
# Args: 
# - beta: beta values
# - detectionP: detectionP values
#	- detectionPval.threshold: detection p-value threshold, set the threshold of significant signals (ok if p-value < detectionPval.threshold). Set to 0.01 by default.
#	- detectionPval.perc.threshold: threshold for the percentage of significant probes per sample. If a sample has less than this threshold it is coonsidered as bad quality. Set to 80% by defult.
#	- projectName: a string used to name the report file about sample quality. If NULL (by default) no report provided.
#	- PATH: a string setting the path for the folder where the report will be saved. Set to current folder by default.
#
# Returns: an updated methylumi object.


detectionPval.filter_2 <- function(beta, detectionP, detectionPval.threshold=0.01, detectionPval.perc.threshold=80, projectName = NULL, PATH="./"){
  
  #get detection p-values
  detectPval <- detectionP
  #get sample names
  samples <- colnames(detectionP)
  nbSignifPval <- as.vector(colSums(detectionP <= detectionPval.threshold))
  percentSignifPval <- as.vector((colSums(detectionP <= detectionPval.threshold)/nrow(beta))*100)
  
  #for each sample compute the number and % or relevant signal (detection p-value < detectionPval.threshold)
  nrMin <- (nrow(beta)/100)*detectionPval.perc.threshold
  #get "bad" samples indices
  index2remove <- which(colSums(detectPval <= detectionPval.threshold) < nrMin)
  
  if(!is.null(projectName)) write.table(list(samples=samples[index2remove], nbSignifPval=nbSignifPval[index2remove], percentSignifPval=percentSignifPval[index2remove]), file=paste(PATH, projectName, "_non-signifPvalStats-Sample_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
  #remove "bad" samples from methylumi object
  if(length(index2remove)>0){
    beta <- beta[,-index2remove]
    detectPval <- detectPval[,-index2remove]
  }

  samples <- colnames(detectPval)
  index <- order(colSums(detectPval))
  
  nbSignifPval <- as.vector(colSums(detectPval <= detectionPval.threshold))
  percentSignifPval <- as.vector((colSums(detectPval <= detectionPval.threshold)/nrow(beta))*100)
  
  # save sample quality report as text file
  if(!is.null(projectName)) write.table(list(samples=samples[index], nbSignifPval=nbSignifPval[index], percentSignifPval=percentSignifPval[index]), file=paste(PATH, projectName, "_signifPvalStats-Sample_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
  return(list(beta, detectPval))
}

# Performs a probe QC on the basis of detection p-values. Eliminate "bad quality" probes if their proportion (%) of probes associated to significant detection p-values < detectionPval.perc.threshold (1% by default).
#
# Args: 
# - beta: beta values
# - detectionP: detectionP values
#	- detectionPval.threshold: detection p-value threshold, set the threshold of significant signals (ok if p-value < detectionPval.threshold). Set to 0.01 by default.
#	- detectionPval.perc.threshold: threshold for the percentage of significant probes per sample. If a sample has less than this threshold it is coonsidered as bad quality. Set to 80% by defult.
#	- projectName: a string used to name the report file about sample quality. If NULL (by default) no report provided.
#	- PATH: a string setting the path for the folder where the report will be saved. Set to current folder by default.
#
# Returns: an updated methylumi object.


detectionPval.filter2_2 <- function(beta, detectionP, detectionPval.threshold=0.01, detectionPval.perc.threshold=1, projectName = NULL, PATH="./"){
  
  #get detection p-values
  detectPval <- detectionP
  #get sample names
  probeNames <- rownames(detectPval)
  nbSignifPval <- as.vector(rowSums(detectPval <= detectionPval.threshold))
  percentSignifPval <- as.vector((rowSums(detectPval <= detectionPval.threshold)/ncol(beta))*100)
  
  #for each probe compute the number and % or relevant signal (detection p-value < detectionPval.threshold)
  
  nrMax <- (ncol(beta)/100)*detectionPval.perc.threshold
  #get "bad" samples indices
  index2remove <- which(rowSums(detectPval > detectionPval.threshold) > nrMax)
  
  if(!is.null(projectName)) write.table(list(probeNames=probeNames[index2remove], nbSignifPval=nbSignifPval[index2remove], percentSignifPval=percentSignifPval[index2remove]), file=paste(PATH, projectName, "_non-signifPvalStats-Probe_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
  #remove "bad" samples from methylumi object
  if(length(index2remove)>0){
    beta <- beta[-index2remove,]
    detectPval <- detectPval[-index2remove,]
  }
  
  probeNames <- rownames(detectPval)
  index <- order(rowSums(detectPval), decreasing=T)
  
  nbSignifPval <- as.vector(rowSums(detectPval <= detectionPval.threshold))
  percentSignifPval <- as.vector((rowSums(detectPval <= detectionPval.threshold)/ncol(beta))*100)
  
  # save sample quality report as text file
  if(!is.null(projectName)) write.table(list(probeNames=probeNames[index], nbSignifPval=nbSignifPval[index], percentSignifPval=percentSignifPval[index]), file=paste(PATH, projectName, "_signifPvalStats-Probe_threshold", detectionPval.threshold, ".txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
  
  return(list(beta, detectPval))
}