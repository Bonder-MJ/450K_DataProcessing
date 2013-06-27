
# code from minfi minimally modified to not require minfi objects

swan2 <- function (unMeth, meth, qc, alfa) {

#LS#  make a MethylSet containing the data and an RGChannelSet just for the annotation
#    if (is.null(mSet)) 
#        mSet <- preprocessRaw(rgSet)
   if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load minfi package')
    }
   if(!library(IlluminaHumanMethylation450kmanifest, logical.return=TRUE, quietly=TRUE)){
         stop('can\'t load IlluminaHumanMethylation450kmanifest package')
    }

    mSet   <-  new("MethylSet", Meth = meth, Unmeth = unMeth)
	rgSet  <- new('RGChannelSet', annotation="IlluminaHumanMethylation450k")
	
	#IlluminaHumanMethylation450kmanifest@data$TypeI$nCpG<-da[IlluminaHumanMethylation450kmanifest@data$TypeI$Name,'nCpG']  

	#LS#
		
	typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
	typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]

	CpG.counts <- rbind(data.frame(typeI@listData), data.frame(typeII@listData))
	CpG.counts$Name <- as.character(CpG.counts$Name)
	CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
	names(CpG.counts)[2] <- "CpGs"
	#browser()
	counts <- CpG.counts[CpG.counts$Name %in% featureNames(mSet),]
	subset <- min(table(counts$CpGs[counts$Type == "I" & counts$CpGs %in% 1:3]), table(counts$CpGs[counts$Type == "II" & counts$CpGs %in% 1:3]))
	# see email from Jovana Maksimovic to Leo 29 oct 2012

	#LS#
    bg <- bgIntensitySwan.methylumi(qc)
    #LS#

    methData <- getMeth(mSet)
    unmethData <- getUnmeth(mSet)
    normMethData <- NULL
    normUnmethData <- NULL
    xNormSet <- vector("list", 2)
    xNormSet[[1]] <- minfi:::getSubset(counts$CpGs[counts$Type == "I"], 
        subset)
    xNormSet[[2]] <- minfi:::getSubset(counts$CpGs[counts$Type == "II"], 
        subset)
    for (i in 1:ncol(mSet)) {
        message(sprintf("Normalizing array %d of %d\n", i, ncol(mSet)))
        normMethData <- cbind(normMethData, minfi:::normaliseChannel(methData[rownames(methData) %in% 
            counts$Name[counts$Type == "I"], i], methData[rownames(methData) %in% 
            counts$Name[counts$Type == "II"], i], xNormSet, bg[i]))
        normUnmethData <- cbind(normUnmethData, minfi:::normaliseChannel(unmethData[rownames(unmethData) %in% 
            counts$Name[counts$Type == "I"], i], unmethData[rownames(unmethData) %in% 
            counts$Name[counts$Type == "II"], i], xNormSet, bg[i]))
    }
    colnames(normMethData) <- sampleNames(mSet)
    colnames(normUnmethData) <- sampleNames(mSet)
	
    normSet <- mSet
    assayDataElement(normSet, "Meth") <- normMethData
    assayDataElement(normSet, "Unmeth") <- normUnmethData
    normSet@preprocessMethod <- c(sprintf("SWAN (based on a MethylSet preprocesses as '%s'", mSet@preprocessMethod[1]), as.character(packageVersion("minfi")), as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
	
    return(getBeta(normSet, offset= alfa))
}


# code from minfi minimally modified to not require minfi objects

swan2_M <- function (unMeth, meth, qc, alfa) {
  
  #LS#  make a MethylSet containing the data and an RGChannelSet just for the annotation
  #    if (is.null(mSet)) 
  #        mSet <- preprocessRaw(rgSet)
  if(!library(minfi, logical.return=TRUE, quietly=TRUE)){
    stop('can\'t load minfi package')
  }
  if(!library(IlluminaHumanMethylation450kmanifest, logical.return=TRUE, quietly=TRUE)){
    stop('can\'t load IlluminaHumanMethylation450kmanifest package')
  }
  
  mSet   <-  new("MethylSet", Meth = meth, Unmeth = unMeth)
  rgSet  <- new('RGChannelSet', annotation="IlluminaHumanMethylation450k")
  
  #IlluminaHumanMethylation450kmanifest@data$TypeI$nCpG<-da[IlluminaHumanMethylation450kmanifest@data$TypeI$Name,'nCpG']  
  
  #LS#
  
  typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
  typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]
  
  CpG.counts <- rbind(data.frame(typeI@listData), data.frame(typeII@listData))
  CpG.counts$Name <- as.character(CpG.counts$Name)
  CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
  names(CpG.counts)[2] <- "CpGs"
  #browser()
  counts <- CpG.counts[CpG.counts$Name %in% featureNames(mSet),]
  subset <- min(table(counts$CpGs[counts$Type == "I" & counts$CpGs %in% 1:3]), table(counts$CpGs[counts$Type == "II" & counts$CpGs %in% 1:3]))
  # see email from Jovana Maksimovic to Leo 29 oct 2012
  
  #LS#
  bg <- bgIntensitySwan.methylumi(qc)
  #LS#
  
  methData <- getMeth(mSet)
  unmethData <- getUnmeth(mSet)
  normMethData <- NULL
  normUnmethData <- NULL
  xNormSet <- vector("list", 2)
  xNormSet[[1]] <- minfi:::getSubset(counts$CpGs[counts$Type == "I"], 
                                     subset)
  xNormSet[[2]] <- minfi:::getSubset(counts$CpGs[counts$Type == "II"], 
                                     subset)
  for (i in 1:ncol(mSet)) {
    message(sprintf("Normalizing array %d of %d\n", i, ncol(mSet)))
    normMethData <- cbind(normMethData, minfi:::normaliseChannel(methData[rownames(methData) %in% 
                                                                            counts$Name[counts$Type == "I"], i], methData[rownames(methData) %in% 
                                                                                                                            counts$Name[counts$Type == "II"], i], xNormSet, bg[i]))
    normUnmethData <- cbind(normUnmethData, minfi:::normaliseChannel(unmethData[rownames(unmethData) %in% 
                                                                                  counts$Name[counts$Type == "I"], i], unmethData[rownames(unmethData) %in% 
                                                                                                                                    counts$Name[counts$Type == "II"], i], xNormSet, bg[i]))
  }
  colnames(normMethData) <- sampleNames(mSet)
  colnames(normUnmethData) <- sampleNames(mSet)
  
  normSet <- mSet
  assayDataElement(normSet, "Meth") <- normMethData
  assayDataElement(normSet, "Unmeth") <- normUnmethData
  normSet@preprocessMethod <- c(sprintf("SWAN (based on a MethylSet preprocesses as '%s'", mSet@preprocessMethod[1]), as.character(packageVersion("minfi")), as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
  
  uns <- getUnmeth(mSet)
  mns <- getMeth(mSet)
  
  indexNegU <- which(!is.numeric(uns), arr.ind=TRUE)
  indexNegM <- which(!is.numeric(mns), arr.ind=TRUE)
  if(length(indexNegU)>0 || length(indexNegM)>0){
    cat("\tWarning: NA values introduced. Value is replaced by 0.\n")
    
    uns[indexNegU] <- 0
    mns[indexNegM] <- 0
  }
  
  #check and "correct" for negative values
  indexNegU <- which(uns < 0, arr.ind=TRUE)
  indexNegM <- which(mns < 0, arr.ind=TRUE)
  uns[indexNegU] <- 0
  mns[indexNegM] <- 0


  return(log2(mns+alfa/uns+alfa))
}


bgIntensitySwan.methylumi  <- function (rg) { # rg is a list of 
                                              # red and green matrices
    # rg <- intensitiesByChannel(QCdata(methylumi_obj))
    neg <- grep("NEGATIVE",rownames(rg$Cy3)) 

    grnMed <- colMedians(rg$Cy3[neg, ])
    redMed <- colMedians(rg$Cy5[neg, ])
    return(rowMeans(cbind(grnMed, redMed)))
}

