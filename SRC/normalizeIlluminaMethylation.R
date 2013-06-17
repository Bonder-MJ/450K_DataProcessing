# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function manages and launchs the subset based quantile normalization for annotation based probe categories.
#
# Args: 
#	- beta: beta values matrix with columns representing samples and rows representing probes.
#	- detect.pval: detection p-values matrix with columns representing samples and rows representing probes.
#	- quantile.norm.pvalThreshold: numeric, represents the threshold for significant detection p-values (significant if < quantile.norm.pvalThreshold). Set to 0.01 by default.
#	- probeAnnotations: list or matrix of Illumina probe annotations.
#	- probeAnnotationsCategory: a string specifying which kind of probe annotation to use for probe categories construction, must be one of "relationToCpG" or "relationToSequence".
#
# Returns: a list of two matrices, 'beta' matrix with normalized beta values and 'detection.pvalue' with updated detection p-values.
#
normalizeIlluminaMethylationSQN <- function(
	beta,
	detect.pval,
	quantile.norm.pvalThreshold = 0.01,
	probeAnnotations,
	probeAnnotationsCategory = "relationToCpG"
	)
{
	
	#get probe IDs for each kind of Illumina bead type (Infinium I & Infinium II
	infiniumI <- probeAnnotations$TargetID[which(probeAnnotations$INFINIUM_DESIGN_TYPE == "I")]
	infiniumII <- probeAnnotations$TargetID[which(probeAnnotations$INFINIUM_DESIGN_TYPE =="II")]
	
	cat("\t Quantile normalization of samples: separated and 'robust' quantile normalization for Infinium probes I and II through probe categories (reference quantiles computed from filtered Infinium I probes only and for different categories of probe annotations). \n")
	if(!is.null(probeAnnotationsCategory)){
		if(probeAnnotationsCategory == "relationToCpG") {
			index <- which(is.element(colnames(probeAnnotations), c("TargetID", "RELATION_TO_UCSC_CPG_ISLAND")))
			probeAnnotations <- probeAnnotations[,index]
		}
		if(probeAnnotationsCategory == "relationToSequence"){
			index <- which(is.element(colnames(probeAnnotations), c("TargetID", "UCSC_REFGENE_GROUP")))
			probeAnnotations <- probeAnnotations[,index]
		}
		if(!is.element(probeAnnotationsCategory, c("relationToCpG", "relationToSequence"))){
			print("WARNINGS: probe annotation category must be one of 'relationToCpG' or 'relationToSequence'.")
			return("WARNINGS: probe annotation category must be on of 'relationToCpG' or 'relationToSequence'.")
		}
	}
	else{
		print("WARNING ! You have to specify an annotation type for probe categories based normalization ('relationToCpG' or 'relationToSequence').")
		return("WARNING ! You have to specify an annotation type for probe categories based normalization ('relationToCpG' or 'relationToSequence').")
	}

	#start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
	data.norm <- robustQuantileNorm_Illumina450K(data = beta, infiniumI = infiniumI, infiniumII = infiniumII, detect.pval = detect.pval, detect.pval.threshold=quantile.norm.pvalThreshold, annotations = probeAnnotations)
	
	names(data.norm) <- c("beta", "detection.pvalue")
	rm(beta, detect.pval)
	rm(infiniumI, infiniumII, probeAnnotations)
	
	cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
	cat("Dimension of normalized detection p-values matrix: ", dim(data.norm$detection.pvalue)[1], "x", dim(data.norm$detection.pvalue)[2],"\n")
		
	return(data.norm)
}

normalizeIlluminaMethylationSWAN <- function(
	detect.pval,
	unMeth,
	meth,
	qc,
	alfa,
	betweenSampleCorrection = TRUE
	)
{	
	normalizedBeta <- swan2(unMeth, meth, qc, alfa)
	
	if(any(is.na(normalizedBeta))){
		for(i in 1: ncol(normalizedBeta)){
			if(any(is.na(normalizedBeta[,i]))){
				ids <- which(!is.numeric(normalizedBeta[,i]))
				med <- median(normalizedBeta[-ids,i])
				normalizedBeta[ids,i] <- med
			}
		}
	}
	
	##Do QN over samples after SWAN?
	if(betweenSampleCorrection==TRUE){
		normalizedBeta2 <- normalize.quantiles(as.matrix(normalizedBeta))
		rownames(normalizedBeta2) <- rownames(normalizedBeta)
		colnames(normalizedBeta2) <- colnames(normalizedBeta)
		normalizedBeta <- normalizedBeta2
		rm(normalizedBeta2)
	}
	#start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
	data.norm <- list(normalizedBeta, detect.pval)
	
	names(data.norm) <- c("beta", "detection.pvalue")
	
	cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
		
	return(data.norm)
}

normalizeIlluminaMethylationBMIQ <- function(
	beta,
	detect.pval,
	quantile.norm.pvalThreshold = 0.01,
	probeAnnotations,
	probeAnnotationsCategory = "relationToCpG",
	QCplot = TRUE,
	betweenSampleCorrection = FALSE
	)
{	
	
	bmiqAnnotation <- cbind(rownames(probeAnnotations), as.character(probeAnnotations$INFINIUM_DESIGN_TYPE))
	
	bmiqAnnotation <- bmiqAnnotation[which( bmiqAnnotation[,1] %in% rownames(beta)),]

	beta <- beta[order(rownames(beta)),]
	bmiqAnnotation <- bmiqAnnotation[order(bmiqAnnotation[,1]),]
	
	#loop through all samples
	rm(probeAnnotations)
	normalizedBeta <- beta;
	
	##Something strange with sample 34, 41
	##166/278
	
	for(k in 1: ncol(beta)){
	#foreach(i=1:ncol(beta)) %dopar% {
		#source(paste(PATH_SRC,"Additions\\BMIQ_1.1_Pipeline.R", sep=""))
		dataOut <- BMIQ(as.vector(unlist(beta[,k])), as.vector(unlist(bmiqAnnotation[,2])), 3, FALSE, nfit=10000, th1.v=c(0.2,0.75), th2.v=NULL, niter=5, tol=0.001, plots=QCplot, sampleID=k)
		normalizedBeta[,k] <- as.array(dataOut[[1]])
		rm(dataOut)

	}
	rm(beta);
	##Do QN over samples after BMIQ?
	if(betweenSampleCorrection==TRUE){
		normalizedBeta2 <- normalize.quantiles(as.matrix(normalizedBeta))
		rownames(normalizedBeta2) <- rownames(normalizedBeta)
		colnames(normalizedBeta2) <- colnames(normalizedBeta)
		normalizedBeta <- normalizedBeta2
		rm(normalizedBeta2)
	}
	
	#start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
	data.norm <- list(normalizedBeta, detect.pval)
	
	names(data.norm) <- c("beta", "detection.pvalue")
	cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")

	return(data.norm)
}

normalizeIlluminaMethylationDASEN <- function(
	detect.pval,
	unMeth,
	meth,
	qc,
	annotation,
	alfa,
	MvalueConv,
	betweenSampleCorrection = TRUE
	)
{	

	ann <- as.array((annotation$INFINIUM_DESIGN_TYPE))
	rownames(ann) <- annotation$TargetID
  keepList <- which(rownames(ann) %in% rownames(unMeth))
	ann <- ann[keepList]
  
  unMeth <- unMeth[order(rownames(unMeth)),]
	meth <- meth[order(rownames(meth)),]
	ann <- ann[order(rownames(ann))]
  
  rm(annotation, keepList)
  
	normalizedBeta <- dasen(uns=unMeth, mns=meth, onetwo=ann, MvalueConv=MvalueConv)
	
	if(any(is.na(normalizedBeta))){
		for(i in 1: ncol(normalizedBeta)){
			if(any(is.na(normalizedBeta[,i]))){
				ids <- which(is.na(normalizedBeta[,i]))
				med <- median(normalizedBeta[-ids,i])
				normalizedBeta[ids,i] <- med
			}
		}
	}
	
	##Do QN over samples after DASEN?
	if(betweenSampleCorrection==TRUE){
		normalizedBeta2 <- normalize.quantiles(as.matrix(normalizedBeta))
		rownames(normalizedBeta2) <- rownames(normalizedBeta)
		colnames(normalizedBeta2) <- colnames(normalizedBeta)
		normalizedBeta <- normalizedBeta2
		rm(normalizedBeta2)
	}
	data.norm <- list(normalizedBeta, detect.pval)
	
	names(data.norm) <- c("beta", "detection.pvalue")
	
	cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
		
	return(data.norm)
}

normalizeIlluminaMethylationMValCor <- function(
  beta,
  detect.pval,
  quantile.norm.pvalThreshold = 0.01,
  probeAnnotations,
  probeAnnotationsCategory = "relationToCpG",
  QCplot = TRUE,
  PATH_RES,
  betweenSampleCorrection = FALSE,
  medianReplacement = TRUE
)
{	
  
  MValAnnotation <- cbind(rownames(probeAnnotations), as.character(probeAnnotations$INFINIUM_DESIGN_TYPE))
  colnames(MValAnnotation) <- c("TargetID", "INFINIUM_DESIGN_TYPE")
  #loop through all samples
  rm(probeAnnotations)
  
  normalizedBeta <- MvalType2Cor(data1=beta, type1=MValAnnotation, PATH_RES=PATH_RES, QCplot=QCplot, medianReplacement=medianReplacement)
  rownames(normalizedBeta) <- normalizedBeta[,1]
  normalizedBeta <- normalizedBeta[,2:ncol(normalizedBeta)]
  colnames(normalizedBeta) <- colnames(beta)
  
  rm(beta);
  ##Do QN over samples after normalization?
  if(betweenSampleCorrection==TRUE){
    normalizedBeta2 <- normalize.quantiles(as.matrix(normalizedBeta))
    rownames(normalizedBeta2) <- rownames(normalizedBeta)
    colnames(normalizedBeta2) <- colnames(normalizedBeta)
    normalizedBeta <- normalizedBeta2
    rm(normalizedBeta2)
  }
  
  #start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
  data.norm <- list(normalizedBeta, detect.pval)
  
  names(data.norm) <- c("beta", "detection.pvalue")
  cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
  
  return(data.norm)
}

normalizeIlluminaMethylationMValCor2 <- function(
  unMeth,
  meth,
  detect.pval,
  probeAnnotations,
  QCplot = TRUE,
  PATH_RES,
  alfa,
  betweenSampleCorrection = FALSE,
  medianReplacement = TRUE
)
{	
  
  MValAnnotation <- cbind(rownames(probeAnnotations), as.character(probeAnnotations$INFINIUM_DESIGN_TYPE))
  colnames(MValAnnotation) <- c("TargetID", "INFINIUM_DESIGN_TYPE")
  #loop through all samples
  rm(probeAnnotations)
  
  normalizedBeta <- MvalType2Cor2(dataU=unMeth , dataM = meth, type1=MValAnnotation, alfa=alfa, PATH_RES=PATH_RES, QCplot=QCplot, medianReplacement=medianReplacement)
  
  rownames(normalizedBeta) <- normalizedBeta[,1]
  normalizedBeta <- normalizedBeta[,2:ncol(normalizedBeta)]
  colnames(normalizedBeta) <- colnames(beta)
  
  rm(unMeth, meth);
  ##Do QN over samples after M-val2?
  if(betweenSampleCorrection==TRUE){
    normalizedBeta2 <- normalize.quantiles(as.matrix(normalizedBeta))
    rownames(normalizedBeta2) <- rownames(normalizedBeta)
    colnames(normalizedBeta2) <- colnames(normalizedBeta)
    normalizedBeta <- normalizedBeta2
    rm(normalizedBeta2)
  }
  
  #start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
  data.norm <- list(normalizedBeta, detect.pval)
  
  names(data.norm) <- c("beta", "detection.pvalue")
  cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
  
  return(data.norm)
}

normalizeIlluminaMethylationSWAN2 <- function(
  detect.pval,
  unMeth,
  meth,
  qc,
  alfa,
  betweenSampleCorrection = TRUE
)
{	
  options(warn=-1)
  normalizedBeta <- swan2_M(unMeth, meth, qc, alfa)
  options(warn=0)
  if(any(is.na(normalizedBeta))){
    for(i in 1: ncol(normalizedBeta)){
      if(any(is.na(normalizedBeta[,i]))){
        ids <- which(!is.numeric(normalizedBeta[,i]))
        med <- median(normalizedBeta[-ids,i])
        normalizedBeta[ids,i] <- med
      }
    }
  }
  
  ##Do QN over samples after SWAN?
  if(betweenSampleCorrection==TRUE){
    normalizedBeta2 <- normalize.quantiles(as.matrix(normalizedBeta))
    rownames(normalizedBeta2) <- rownames(normalizedBeta)
    colnames(normalizedBeta2) <- colnames(normalizedBeta)
    normalizedBeta <- normalizedBeta2
    rm(normalizedBeta2)
  }
  #start subset quantile normalization, this function returns a list of 2 matrices (beta values & detection p-values)
  data.norm <- list(normalizedBeta, detect.pval)
  
  names(data.norm) <- c("beta", "detection.pvalue")
  
  cat("\nDimension of normalized beta values matrix: ", dim(data.norm$beta)[1], "x", dim(data.norm$beta)[2],"\n")
  
  return(data.norm)
}