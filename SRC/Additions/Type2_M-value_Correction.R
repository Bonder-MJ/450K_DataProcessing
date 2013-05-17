###################################################################
# R script for typeII correction.                                 #
# Based on the method of Dedeurwaerder et al.                     #
# www.futuremedicine.com/doi/pdfplus/10.2217/epi.11.105           #
# Implemented by Pan and Chen, et al.                             #
# www.landesbioscience.com/journals/epigenetics/2012EPI0192R1.pdf #
###################################################################

MvalType2Cor <- function(data1, type1, idNameY = 1, PATH_RES, QCplot=FALSE){

	data1 <- cbind(rownames(data1), as.data.frame(data1))
	colnames(data1)[1] <- as.character(colnames(type1)[idNameY])
	id1 = which(colnames(type1)=="INFINIUM_DESIGN_TYPE")
	NameY = as.character(colnames(type1)[idNameY])
	type1 <- as.data.frame(type1[,c(idNameY, id1)])
	
	data1 <- merge(data1, type1, by.x=NameY, by.y= NameY)
	
	rm(type1, NameY, id1);
	
	head(data1); dim(data1)

	naidx <- apply(is.na(data1), 1, sum)
	data1 <- data1[which(naidx == 0),]
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
			plot(density(beta_typeI, bw=0.05), main='Fig1: Beta value for typeI and typeII(red)')
			lines(density(beta_typeII, bw=0.05), col='red')
			plot(density(m2, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE), col='red', main='Fig2: M value for typeI and typeII(red)')
			lines(density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE))
		}
		#typeI M-value as the base.
		dm1 <- density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE)
		sigma_m1 <- dm1$x[which.max(dm1$y[dm1$x >= 0])+ length(dm1$x[dm1$x < 0])]
		sigma_u1 <- dm1$x[which.max(dm1$y[dm1$x < 0])]
		cat(sigma_m1, "; ", sigma_u1, "\n")
		
		#######adjust typeII M value according to typeI
		dm2 <- density(m2, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE)
		sigma_m2 <- dm2$x[which.max(dm2$y[dm2$x >= 0])+ length(dm2$x[dm2$x < 0])]
		sigma_u2 <- dm2$x[which.max(dm2$y[dm2$x < 0])]
		cat(sigma_m2, "; ", sigma_u2, "\n")
		mm2 <- rep(0, probe_typeII)
		idx11 <- which(m2 >= 0)
		idx22 <- which(m2 < 0)
		mm2[idx11]<- (m2[idx11]/sigma_m2)*sigma_m1
		mm2[idx22]<- (m2[idx22]/sigma_u2)*sigma_u1
		beta_typeII_corrected <- 2^mm2/(2^mm2+1)
		
		if(QCplot){
			plot(density(mm2,bw=0.5, kernel="gaussian", n=200, na.rm=TRUE), col='red', main="Fig3. M value in typeI and adjust typeII(red)")
			lines(density(m1, bw=0.5, kernel="gaussian", n=200, na.rm=TRUE))
			plot(density(beta_typeI, bw=0.05), main='Fig4: Beta value in typeI and corrected typeII(red)')
			lines(density(beta_typeII_corrected, bw=0.05), col='red')
			dev.off()
		}
		
		temp <- c(beta_typeI, beta_typeII_corrected)
		beta1_corrected <- data.frame(beta1_corrected, temp)
	}
	
	return(beta1_corrected)
}
