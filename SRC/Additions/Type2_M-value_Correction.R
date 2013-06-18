###################################################################
# R script for typeII correction.                                 #
# Based on the method of Dedeurwaerder et al.                     #
# www.futuremedicine.com/doi/pdfplus/10.2217/epi.11.105           #
# Implemented by Pan and Chen, et al.                             #
# www.landesbioscience.com/journals/epigenetics/2012EPI0192R1.pdf #
###################################################################

MvalType2Cor <- function(data1, type1, idNameY = 1, PATH_RES, QCplot=FALSE, medianReplacement=TRUE){

	if(medianReplacement){
		medianD1 <- median(which(is.numeric(data1)))
		data1[which(!is.numeric(data1))] <- medianD1
	}
	
	data1 <- cbind(rownames(data1), as.data.frame(data1))
	colnames(data1)[1] <- as.character(colnames(type1)[idNameY])
	id1 = which(colnames(type1)=="INFINIUM_DESIGN_TYPE")
	NameY = as.character(colnames(type1)[idNameY])
	type1 <- as.data.frame(type1[,c(idNameY, id1)])
	
	data1 <- base::merge(data1, type1, by.x=NameY, by.y= NameY)
	data1 <- as.data.frame(data1)
	rm(type1, NameY, id1);
	
	head(data1); dim(data1)
  
	naidx <- apply(data1[2:(length(data1)-1)], 1, function(x){ return(all(is.numeric(x)))})
	data1 <- data1[which(naidx==TRUE),]
  
	row_size <- dim(data1)[1]
	col_size <- dim(data1)[2] - 1
	idx1 <- which(data1$INFINIUM_DESIGN_TYPE=='I')
	idx2 <- which(data1$INFINIUM_DESIGN_TYPE=='II')
  
	beta1_typeI_probe <- as.character(data1[idx1,1])
	beta1_typeII_probe <- as.character(data1[idx2,1])
	len <- length(colnames(data1))
	beta1_corrected <- c(as.character(beta1_typeI_probe), as.character(beta1_typeII_probe))

	for (i in 2:col_size){
		beta1 <- data1[,i] #one sample
		cat(colnames(data1)[i], "; ")
		beta_typeI <- beta1[idx1]
		beta_typeII <- beta1[idx2]
		m1 <- log2(beta_typeI/(1 - beta_typeI))
		m2 <- log2(beta_typeII/(1 - beta_typeII))
		probe_typeII <- length(m2)
		if(QCplot){
			pdf(paste(PATH_RES, "./QC_Plot_" , (i-1), ".pdf", sep=""))
			#plot the density of typeI and typeII beta value.
			#dev.new(1)
			
			par(mfrow=c(2,2))
			plot(stats::density(beta_typeI, bw=0.05), main='Fig1: Beta value for typeI and typeII(red)')
			lines(stats::density(beta_typeII, bw=0.05), col='red')
			plot(stats::density(m2, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE), col='red', main='Fig2: M value for typeI and typeII(red)')
			lines(stats::density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE))
		}
		#typeI M-value as the base.
		dm1 <- stats::density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE)
		sigma_m1 <- dm1$x[which.max(dm1$y[dm1$x >= 0])+ length(dm1$x[dm1$x < 0])]
		sigma_u1 <- dm1$x[which.max(dm1$y[dm1$x < 0])]
		cat(sigma_m1, "; ", sigma_u1, "\n")
		
		#######adjust typeII M value according to typeI
		dm2 <- stats::density(m2, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE)
		sigma_m2 <- dm2$x[which.max(dm2$y[dm2$x >= 0])+ length(dm2$x[dm2$x < 0])]
		sigma_u2 <- dm2$x[which.max(dm2$y[dm2$x < 0])]
		cat(sigma_m2, "; ", sigma_u2, "\n")
		mm2 <- rep(0, probe_typeII)
		idx11 <- which(m2 >= 0)
		idx22 <- which(m2 < 0)
    
		if(length(idx11)>0 && (length(sigma_m2)==0 || length(sigma_m1)==0)){
		  cat("\tWarning: removed (bad) sample: ",colnames(data1)[i], "\n")
		  cat("\tChange alfa to overcome this problem.\n")
		  next
		}
		
		mm2[idx11]<- (m2[idx11]/sigma_m2)*sigma_m1
		
		if(length(idx22)>0 && (length(sigma_u2)==0 || length(sigma_u1)==0)){
		  cat("\tWarning: removed (bad) sample: ",colnames(data1)[i], "\n")
		  cat("\tChange alfa to overcome this problem.\n")
		  next
		}
		   
		   mm2[idx22]<- (m2[idx22]/sigma_u2)*sigma_u1
    
    
		beta_typeII_corrected <- 2^mm2/(2^mm2+1)
		
		if(QCplot){
			plot(stats::density(mm2,bw=0.5, kernel="gaussian", n=200, na.rm=TRUE), col='red', main="Fig3. M value in typeI and adjust typeII(red)")
			lines(stats::density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE))
			plot(stats::density(beta_typeI, bw=0.05), main='Fig4: Beta value in typeI and corrected typeII(red)')
			lines(stats::density(beta_typeII_corrected, bw=0.05), col='red')
			dev.off()
		}
		
		temp <- c(beta_typeI, beta_typeII_corrected)
		beta1_corrected <- data.frame(beta1_corrected, temp)
	}
	
	return(beta1_corrected)
}

MvalType2Cor2 <- function(dataU, dataM, type1, idNameY = 1, alfa=10, PATH_RES, QCplot=FALSE, medianReplacement=FALSE){
	options(warn=-1)
  
	data1 <- log2(dataM+alfa/dataU+alfa)
	
	if(medianReplacement){
		medianD1 <- median(which(is.numeric(data1)))
		data1[which(!is.numeric(data1))] <- medianD1
	}
	
	rm(dataU, dataM)
	
	data1 <- cbind(rownames(data1), as.data.frame(data1))
	colnames(data1)[1] <- as.character(colnames(type1)[idNameY])
	id1 = which(colnames(type1)=="INFINIUM_DESIGN_TYPE")
	NameY = as.character(colnames(type1)[idNameY])
	type1 <- as.data.frame(type1[,c(idNameY, id1)])
	
	data1 <- base::merge(data1, type1, by.x=NameY, by.y= NameY)
	
	rm(type1, NameY, id1);
	
	head(data1); dim(data1)

	naidx <- apply(data1[2:(length(data1)-1)], 1, function(x){ return(all(is.numeric(x)))})
	data1 <- data1[which(naidx==TRUE),]
	
	row_size <- dim(data1)[1]
	col_size <- dim(data1)[2] - 1
	idx1 <- which(data1$INFINIUM_DESIGN_TYPE=='I')
	idx2 <- which(data1$INFINIUM_DESIGN_TYPE=='II')
	M_typeI_probe <- as.character(data1[idx1,1])
	M_typeII_probe <- as.character(data1[idx2,1])
	len <- length(colnames(data1))
	M_corrected <- c(as.character(M_typeI_probe), as.character(M_typeII_probe))

	for (i in 2:col_size){
		
		Mval1 <- data1[,i] #one sample
		
		cat(colnames(data1)[i], "; ")
		m1 <- Mval1[idx1]
		m2 <- Mval1[idx2]
		probe_typeII <- length(m2)
		
		if(QCplot){
			pdf(paste(PATH_RES, "./QC_Plot_" , (i-1), ".pdf", sep=""))
			#plot the density of typeI and typeII beta value.
			#dev.new(1)
			
			par(mfrow=c(2,1))
			plot(stats::density(m2, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE), col='red', main='Fig1: M value for typeI and typeII(red)')
			lines(stats::density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE))
		}

		#typeI M-value as the base.
		dm1 <- stats::density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE)
		sigma_m1 <- dm1$x[which.max(dm1$y[dm1$x >= 0])+ length(dm1$x[dm1$x < 0])]
		sigma_u1 <- dm1$x[which.max(dm1$y[dm1$x < 0])]
		cat(sigma_m1, "; ", sigma_u1, "\n")

		#######adjust typeII M value according to typeI
		dm2 <- stats::density(m2, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE)
		sigma_m2 <- dm2$x[which.max(dm2$y[dm2$x >= 0])+ length(dm2$x[dm2$x < 0])]
		sigma_u2 <- dm2$x[which.max(dm2$y[dm2$x < 0])]
		cat(sigma_m2, "; ", sigma_u2, "\n")
		mm2 <- rep(0, probe_typeII)
		idx11 <- which(m2 >= 0)
		idx22 <- which(m2 < 0)
		
		if(length(idx11)>0 && (length(sigma_m2)==0 || length(sigma_m1)==0)){
		  cat("\tWarning: removed (bad) sample: ",colnames(data1)[i], "\n")
		  cat("\tChange alfa to overcome this problem.\n")
		  next
		}
    
		mm2[idx11]<- (m2[idx11]/sigma_m2)*sigma_m1
    
		if(length(idx22)>0 && (length(sigma_u2)==0 || length(sigma_u1)==0)){
      cat("\tWarning: removed (bad) sample: ",colnames(data1)[i], "\n")
      cat("\tChange alfa to overcome this problem.\n")
		  next
		}
    
		mm2[idx22]<- (m2[idx22]/sigma_u2)*sigma_u1
		
		if(QCplot){
			plot(stats::density(mm2,bw=0.5, kernel="gaussian", n=200, na.rm=TRUE), col='red', main="Fig2. M value in typeI and adjust typeII(red)")
			lines(stats::density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE))
			dev.off()
		}
		temp <- c(m1, mm2)
		M_corrected <- data.frame(M_corrected, temp)
	}
	options(warn=0)

	return(M_corrected)
}

