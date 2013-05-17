# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function performs a complete Illumina 450K array data preprocessing and subset quantile normalization for a data batch (set of 450K "plates").
#
pipelineIlluminaMethylation.batch <- function(
	PATH_PROJECT_DATA,
	projectName,
	nbBeads.threshold,
	detectionPval.threshold,
	detectionPval.perc.threshold,
	sampleSelection,
	probeSNP_LIST,
	XY.filtering,
	colorBias.corr,
	bg.adjust,
	PATH,
	QCplot=TRUE,
	betweenSampleCorrection = betweenSampleCorrection,
	alfa=100,
	NormProcedure
 ){

	####################################
	# Get project pathes and load data #
	####################################

	subProjects <- dir(PATH_PROJECT_DATA)
	
	beta <- NULL
	
	#for all subprojects
	for(i in 1:length(subProjects)){
		projectName_batch <- subProjects[i]
		sampleTable <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern="TableSample")
		controlTable <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern="TableControl")
		cat("\n# ")
		cat("Processing sub-project: ", projectName_batch, "\n")

		if(length(sampleTable) > 1){
			warning <- "\tWARNING ! Sample table: too many files matching with pattern 'TableSample' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(sampleTable) < 1){
			warning <- "\tWARNING ! Sample table: no file matching with pattern 'TableSample' ! \n"
			cat(warning)
			return(NULL)
		}
		if(length(controlTable) > 1){
			warning <- "\tWARNING ! Control table: too many files matching with pattern 'TableControl' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(controlTable) < 1){
			warning <- "\tWARNING ! Control table: no file matching with pattern 'TableControl' ! \n"
			cat(warning)
			return(warning)	
		}

		if(sampleSelection){
			sampleList <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern = "sampleList")
			path2sampleList <- paste(PATH_PROJECT_DATA, projectName_batch, "/", sampleList, sep="")

			if(length(sampleList) > 1){
				warning <- "\tWARNING ! List for sample selection: too many files matching with pattern 'SampleList' ! \n"
				cat(warning)
				return(warning)
			}
			if(length(sampleList) < 1){
				warning <- "\tWARNING ! List for sample selection: no file matching with pattern 'SampleList' ! \n"
				path2sampleList <- NULL
				cat(warning)
			}
		}
		else{path2sampleList <- NULL}

		path2data <- paste(PATH_PROJECT_DATA, projectName_batch, "/", sampleTable, sep="")
		path2controlData <- paste(PATH_PROJECT_DATA, projectName_batch, "/", controlTable, sep="")
		
		cat("\tSample table: ", path2data, "\n")
		cat("\tControl table: ", path2controlData, "\n")
		cat("\tSample list (for sample selection): ", path2sampleList, "\n")
			
		#############################
		# starts data preprocessing #
		#############################

		methLumi_data <- preprocessIlluminaMethylation(
			path2data = path2data,
			path2controlData = path2controlData,
			projectName = projectName_batch,
			nbBeads.threshold = nbBeads.threshold,
			detectionPval.threshold = detectionPval.threshold,
			detectionPval.perc.threshold = detectionPval.perc.threshold,
			sample2keep = path2sampleList,
			probeSNP_LIST,
			XY.filtering = XY.filtering,
			colorBias.corr = colorBias.corr,
			bg.adjust = bg.adjust,
			PATH = PATH_RES,
			QCplot = QCplot
		)
		
		if(is.null(methLumi_data)){
			next;
		}

		################################################
		# Sub-project data & information concatenation #
		################################################

		if(i == 1){
			beta <- getMethylumiBeta(methLumi_data, alfa)
			cat("\t beta plate", i, " ok (", dim(beta)[1], "x", dim(beta)[2], ").\n")
			detectionPval <- assayDataElement(methLumi_data, "detection")
			cat("\t detection p-values plate", i, " ok (", dim(detectionPval)[1], "x", dim(detectionPval)[2], ").\n")
			#select "useful" probe annotations
			annotation <- fData(methLumi_data) ; rm(methLumi_data)
			index <- which(is.element(colnames(annotation), c("TargetID", "INFINIUM_DESIGN_TYPE", "RELATION_TO_UCSC_CPG_ISLAND", "UCSC_REFGENE_GROUP")))
			annotation <- annotation[,index]
		} else {
			#concatenate 'betas'
			beta_i <- getMethylumiBeta(methLumi_data, alfa)
			cat("\t beta_", i, " ok (", dim(beta_i)[1], "x", dim(beta_i)[2], ").\n")
			detectionPval_i <- assayDataElement(methLumi_data, "detection")
			cat("\t For all sub-projects: beta matrices concatenation & detection p-value matrices concatenation.\n")

			beta <- concatenateMatrices(beta, beta_i) ; rm(beta_i)
			detectionPval <- concatenateMatrices(detectionPval, detectionPval_i) ; rm(detectionPval_i)
			annotation <- annotation[ which(is.element(annotation$TargetID, rownames(beta))),]

			cat("\t beta ok (", dim(beta)[1], "x", dim(beta)[2], ").\n")
			cat("\t detection p-values ok (", dim(detectionPval)[1], "x", dim(detectionPval)[2], ").\n")
		}
	}
	
	if(is.null(beta)){
		return(NULL)
	}
	
	
	############################################################################################
	# Extraction of SNP probes ("rs" probes)
	############################################################################################
	
	indexSNP <- grep(pattern="rs*", x=rownames(beta))
	if(length(indexSNP)>0){
		betaSNP <- beta[indexSNP,]
		detectionPvalSNP <- detectionPval[indexSNP,]
		
		beta <- beta[-indexSNP,]
		detectionPval <- detectionPval[-indexSNP,]
		
		write.table(betaSNP, file=paste(PATH_RES, projectName, "_betaSNPprobes.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
		write.table(detectionPvalSNP, file=paste(PATH_RES, projectName, "_detectionPvalueSNPprobes.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
	}
	
	############################################################################################
	# start data normalization (subset quantile normalization per probe annotation categories) #
	############################################################################################
	write.table(beta, file=paste(PATH_RES, projectName, "_beta_intermediate.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
	
	if(NormProcedure == "BMIQ"){
		data.preprocess.norm <- normalizeIlluminaMethylationBMIQ(
			beta = beta,
			detect.pval = detectionPval,
			quantile.norm.pvalThreshold = detectionPval.threshold,
			probeAnnotations = annotation,
			probeAnnotationsCategory = probeAnnotationsCategory,
			QCplot = QCplot,
			betweenSampleCorrection = betweenSampleCorrection
		)
	} else if(NormProcedure == "M-ValCor"){
		data.preprocess.norm <- normalizeIlluminaMethylationMValCor(
			beta = beta,
			detect.pval = detectionPval,
			quantile.norm.pvalThreshold = detectionPval.threshold,
			probeAnnotations = annotation,
			probeAnnotationsCategory = probeAnnotationsCategory,
			QCplot = QCplot,
			PATH_RES,
			betweenSampleCorrection = betweenSampleCorrection
		)
	} else {
		data.preprocess.norm <- normalizeIlluminaMethylationSQN(
			beta = beta,
			detect.pval = detectionPval,
			quantile.norm.pvalThreshold = detectionPval.threshold,
			probeAnnotations = annotation,
			probeAnnotationsCategory = probeAnnotationsCategory
		)
	}

	return(data.preprocess.norm)
}

# This function performs a complete Illumina 450K array data preprocessing and swan for a data batch (set of 450K "plates").
#
pipelineIlluminaMethylation.batch.SWAN <- function(
	PATH_PROJECT_DATA,
	projectName,
	nbBeads.threshold,
	detectionPval.threshold,
	detectionPval.perc.threshold,
	sampleSelection,
	probeSNP_LIST,
	XY.filtering,
	colorBias.corr,
	bg.adjust,
	PATH,
	QCplot=TRUE,
	alfa=100,
	betweenSampleCorrection = TRUE
 ){

	####################################
	# Get project pathes and load data #
	####################################

	subProjects <- dir(PATH_PROJECT_DATA)
	
	beta <- NULL
	
	#for all subprojects
	for(i in 1:length(subProjects)){
		projectName_batch <- subProjects[i]
		sampleTable <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern="TableSample")
		controlTable <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern="TableControl")
		cat("\n# ")
		cat("Processing sub-project: ", projectName_batch, "\n")

		if(length(sampleTable) > 1){
			warning <- "\tWARNING ! Sample table: too many files matching with pattern 'TableSample' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(sampleTable) < 1){
			warning <- "\tWARNING ! Sample table: no file matching with pattern 'TableSample' ! \n"
			cat(warning)
			return(NULL)
		}
		if(length(controlTable) > 1){
			warning <- "\tWARNING ! Control table: too many files matching with pattern 'TableControl' ! \n"
			cat(warning)
			return(warning)
		}
		if(length(controlTable) < 1){
			warning <- "\tWARNING ! Control table: no file matching with pattern 'TableControl' ! \n"
			cat(warning)
			return(warning)	
		}

		if(sampleSelection){
			sampleList <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern = "sampleList")
			path2sampleList <- paste(PATH_PROJECT_DATA, projectName_batch, "/", sampleList, sep="")

			if(length(sampleList) > 1){
				warning <- "\tWARNING ! List for sample selection: too many files matching with pattern 'SampleList' ! \n"
				cat(warning)
				return(warning)
			}
			if(length(sampleList) < 1){
				warning <- "\tWARNING ! List for sample selection: no file matching with pattern 'SampleList' ! \n"
				path2sampleList <- NULL
				cat(warning)
			}
		}
		else{path2sampleList <- NULL}

		path2data <- paste(PATH_PROJECT_DATA, projectName_batch, "/", sampleTable, sep="")
		path2controlData <- paste(PATH_PROJECT_DATA, projectName_batch, "/", controlTable, sep="")
		
		cat("\tSample table: ", path2data, "\n")
		cat("\tControl table: ", path2controlData, "\n")
		cat("\tSample list (for sample selection): ", path2sampleList, "\n")
			
		#############################
		# starts data preprocessing #
		#############################
		
		methLumi_data <- preprocessIlluminaMethylation(
			path2data = path2data,
			path2controlData = path2controlData,
			projectName = projectName_batch,
			nbBeads.threshold = nbBeads.threshold,
			detectionPval.threshold = detectionPval.threshold,
			detectionPval.perc.threshold = detectionPval.perc.threshold,
			sample2keep = path2sampleList,
			probeSNP_LIST,
			XY.filtering = XY.filtering,
			colorBias.corr = colorBias.corr,
			bg.adjust = bg.adjust,
			PATH = PATH_RES,
			QCplot = QCplot
		)
		
		if(is.null(methLumi_data)){
			next;
		}
		

		################################################
		# Sub-project data & information concatenation #
		################################################

		if(i == 1){
			beta <- getMethylumiBeta(methLumi_data, alfa)
			
			unMeth <- unmethylated(methLumi_data)
			meth <- methylated(methLumi_data)
			qc <- intensitiesByChannel(QCdata(methLumi_data)) 
			
			cat("\t beta plate", i, " ok (", dim(beta)[1], "x", dim(beta)[2], ").\n")
			detectionPval <- assayDataElement(methLumi_data, "detection")
			cat("\t detection p-values plate", i, " ok (", dim(detectionPval)[1], "x", dim(detectionPval)[2], ").\n")
			#select "useful" probe annotations
			annotation <- fData(methLumi_data) ; rm(methLumi_data)
			index <- which(is.element(colnames(annotation), c("TargetID", "INFINIUM_DESIGN_TYPE", "RELATION_TO_UCSC_CPG_ISLAND", "UCSC_REFGENE_GROUP")))
			annotation <- annotation[,index]
		} else{
			#concatenate 'betas'
			beta_i <- getMethylumiBeta(methLumi_data, alfa)
			cat("\t beta_", i, " ok (", dim(beta_i)[1], "x", dim(beta_i)[2], ").\n")
			detectionPval_i <- assayDataElement(methLumi_data, "detection")
			cat("\t For all sub-projects: beta matrices concatenation & detection p-value matrices concatenation.\n")
			qc_i <- intensitiesByChannel(QCdata(methLumi_data))
			
			qc[[1]] <- cbind(qc[[1]], qc_i[[1]])
			qc[[2]] <- cbind(qc[[2]], qc_i[[2]])
			
			beta <- concatenateMatrices(beta, beta_i) ; 
			
			detectionPval <- concatenateMatrices(detectionPval, detectionPval_i) ; rm(detectionPval_i)
			annotation <- annotation[ which(is.element(annotation$TargetID, rownames(beta))),]
			
			if(length(rownames(meth)) == length(rownames(methylated(methLumi_data)))){
				if(all(rownames(meth) == rownames(methylated(methLumi_data)))){
					meth <- cbind(meth, methylated(methLumi_data))
					unMeth <- cbind(unMeth, unmethylated(methLumi_data))
				} else {
					unMeth_i <- unmethylated(methLumi_data)
					meth_i <- methylated(methLumi_data)
					
					unMeth_i<- as.matrix(unMeth_i[which(rownames(unMeth_i) %in% rownames(unMeth)),])
					meth_i<- as.matrix(meth_i[which(rownames(meth_i) %in% rownames(meth)),])
					
					unMeth<- as.matrix(unMeth[which(rownames(unMeth) %in% rownames(unMeth_i)),])
					meth<- as.matrix(meth[which(rownames(meth) %in% rownames(meth_i)),])
					
					meth <- cbind(meth[order(rownames(meth)),], meth_i[order(rownames(meth_i)),])
					unMeth <- cbind(unMeth[order(rownames(unMeth)),], unMeth_i[order(rownames(unMeth_i)),])
				}
			} else {
				unMeth_i <- unmethylated(methLumi_data)
				meth_i <- methylated(methLumi_data)
				
				unMeth_i<- as.matrix(unMeth_i[which(rownames(unMeth_i) %in% rownames(unMeth)),])
				meth_i<- as.matrix(meth_i[which(rownames(meth_i) %in% rownames(meth)),])
				
				unMeth<- as.matrix(unMeth[which(rownames(unMeth) %in% rownames(unMeth_i)),])
				meth<- as.matrix(meth[which(rownames(meth) %in% rownames(meth_i)),])
				
				meth <- cbind(meth[order(rownames(meth)),], meth_i[order(rownames(meth_i)),])
				unMeth <- cbind(unMeth[order(rownames(unMeth)),], unMeth_i[order(rownames(unMeth_i)),])
			}
			rm(beta_i)
			rm(qc_i)
			rm(unMeth_i)
			rm(meth_i)
			
			cat("\t beta ok (", dim(beta)[1], "x", dim(beta)[2], ").\n")
			cat("\t detection p-values ok (", dim(detectionPval)[1], "x", dim(detectionPval)[2], ").\n")
		}
	}
	
	if(is.null(beta)){
		return(NULL)
	}
	
	############################################################################################
	# Extraction of SNP probes ("rs" probes)
	############################################################################################
	
	indexSNP <- grep(pattern="rs*", x=rownames(beta))
	if(length(indexSNP)>0){
		betaSNP <- beta[indexSNP,]
		detectionPvalSNP <- detectionPval[indexSNP,]
		
		beta <- beta[-indexSNP,]
		detectionPval <- detectionPval[-indexSNP,]
		
		write.table(betaSNP, file=paste(PATH_RES, projectName, "_betaSNPprobes.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
		write.table(detectionPvalSNP, file=paste(PATH_RES, projectName, "_detectionPvalueSNPprobes.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
	}
	
	indexSNP <- grep(pattern="rs*", x=rownames(unMeth))
	if(length(indexSNP)>0){
		unMeth <- unMeth[-indexSNP,]
		meth <- meth[-indexSNP,]
	}
	
	
	
	############################################################################################								
	# start data normalization (subset quantile normalization per probe annotation categories) #
	############################################################################################
	
	write.table(beta, file=paste(PATH_RES, projectName, "_beta_intermediate.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
	
	data.preprocess.norm <- normalizeIlluminaMethylationSWAN(
		detect.pval = detectionPval,
		unMeth,
		meth,
		qc,
		betweenSampleCorrection = betweenSampleCorrection
	)

	return(data.preprocess.norm)
}

