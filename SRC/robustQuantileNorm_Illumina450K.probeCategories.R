# 2011-2012
# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# Function that performs the subset based quantile normalization of beta values for separate probe categories (according to Illumina annotation)
#
# Args:
#	- data.infiniumI: beta values matrix for InfiniumI probes with columns representing samples and rows representing probes.
#	- data.infiniumII: beta values matrix for InfiniumII probes with columns representing samples and rows representing probes.
#	- detect.pval.infiniumI: detection p-values matrix for InfiniumI probes with columns representing samples and rows representing probes.
#	- annotations: probe annotations matrix
#	- threshold: numeric, detection p-values significance threshold.
#
# Returns: a list of two matrices, corresponding to InfiniumI and InfiniumII normalized beta values.
#

robustQuantileNorm_Illumina450K.probeCategories <- function(data.infiniumI, data.infiniumII, detect.pval.infiniumI, annotations, threshold, verbose=TRUE){

	res.norm <- list()
	#build filtred dataI
	data.detectPvalRemoved.infiniumI <- dataDetectPval2NA(data.infiniumI, detect.pval.infiniumI, threshold)

	#get list of all possible category values
	category <- uniqueAnnotationCategory(annotations[,2])
	if(verbose) cat("\tProbe categories values: ", paste(category, sep="", collapse=", "), sep="")

	#for each category value
	for(i in 1:length(category)){
		if(verbose) cat("\n\tFor ", category[i], ":\n")

		#get list of probes associated to category value i
		probes_i <- findAnnotationProbes(annotationValue = category[i], annotations, uniqueAnnot=TRUE)
		robustProbes_i <- findAnnotationProbes(annotationValue = category[i], annotations, uniqueAnnot=TRUE)
		if(verbose) cat("\t\tNb probes for category", category[i],":", length(probes_i), "\n")

		#build indices of probesI', probesI and probeII for category value i
		probeIndexI_i <- which(is.element(rownames(data.infiniumI), probes_i))
			if(verbose) cat("\t\tNb Inf_I probes for category ", category[i],": ", length(probeIndexI_i), "\n")
		probeIndexII_i <- which(is.element(rownames(data.infiniumII), probes_i))
			if(verbose) cat("\t\tNb Inf_II probes for category ", category[i],": ", length(probeIndexII_i), "\n")
		probe.detectPvalRemoved.IndexI_i <- which(is.element(rownames(data.detectPvalRemoved.infiniumI), robustProbes_i))
			if(verbose) cat("\t\tNb filtred Inf_I' probes for category ", category[i],": ", length(probe.detectPvalRemoved.IndexI_i), "\n")

		#build associated beta value data
		data.infiniumI_i <- data.infiniumI[probeIndexI_i,]
		data.infiniumII_i <- data.infiniumII[probeIndexII_i,]
		data.detectPvalRemoved.infiniumI_i <- data.detectPvalRemoved.infiniumI[probe.detectPvalRemoved.IndexI_i,]

		probeID.infiniumI_i <- rownames(data.infiniumI_i)
		probeID.infiniumII_i <- rownames(data.infiniumII_i)

		#compute ref.quantiles from dataI'
		ref.quantiles.detectPvalRemoved.infI_i <- referenceQuantiles(data.detectPvalRemoved.infiniumI_i)
		indexRefQuantilesNA <- which(is.na(ref.quantiles.detectPvalRemoved.infI_i))
		if(length(indexRefQuantilesNA)>0 && length(indexRefQuantilesNA)<length(ref.quantiles.detectPvalRemoved.infI_i)){
			ref.quantiles.detectPvalRemoved.infI_i <- ref.quantiles.detectPvalRemoved.infI_i[-indexRefQuantilesNA]
		}
		if(length(indexRefQuantilesNA)==length(ref.quantiles.detectPvalRemoved.infI_i)) {
			print(paste("WARNING ! nb of 'NA' in ref.quantiles from dataI' for ", category[i], "equals the nb of quantiles (", length(indexRefQuantilesNA),") !!!!", sep=""))
			return(paste("WARNING ! nb of 'NA' in ref.quantiles from dataI' for ", category[i], "equals the nb of quantiles (", length(indexRefQuantilesNA),") !!!!", sep=""))
		}
		if((length(indexRefQuantilesNA)/length(ref.quantiles.detectPvalRemoved.infI_i))>0.8) {
			print(paste("WARNING !! nb of 'NA' in ref.quantiles from dataI' for ", category[i], "represents more than 80% of the ref.quantiles (", (length(indexRefQuantilesNA)/length(ref.quantiles.detectPvalRemoved.infI_i)*100),") !!!!", sep=""))
		}
			
		#compute ref.quantiles for dataI
		ref.quantiles.infI_i <- adaptRefQuantiles(ref.quantiles.detectPvalRemoved.infI_i, dim(data.infiniumI_i)[1])
		
		#compute ref.quantiles for dataII
		ref.quantiles.infII_i <- adaptRefQuantiles(ref.quantiles.detectPvalRemoved.infI_i, dim(data.infiniumII_i)[1])
			
		#normalize dataI
		data.infiniumI.norm_i <- normalize.quantiles2(data.infiniumI_i, ref.quantiles.infI_i)
		rownames(data.infiniumI.norm_i) <- probeID.infiniumI_i
		
		#normalize dataII
		data.infiniumII.norm_i <- normalize.quantiles2(data.infiniumII_i, ref.quantiles.infII_i)
		rownames(data.infiniumII.norm_i) <- probeID.infiniumII_i
		
		#bind category dataI.norm and category dataII.norm to other data
		res.norm$data.infiniumI.norm <- rbind(res.norm$data.infiniumI.norm, data.infiniumI.norm_i)
		res.norm$data.infiniumII.norm <- rbind(res.norm$data.infiniumII.norm, data.infiniumII.norm_i)
	}
	#if(verbose) print(paste("Dim norm. data_Inf_I: ", dim(res.norm$data.infiniumI.norm), sep="", collapse=" ; "))
	#if(verbose) print(paste("Dim norm. data_Inf_II: ", dim(res.norm$data.infiniumII.norm), sep="", collapse=" ; "))

	return(rbind(res.norm$data.infiniumI.norm, res.norm$data.infiniumII.norm))
}
