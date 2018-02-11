# Marc Jan Bonder
# m.j.bonder @ umcg.nl
#
############################
# Release data: 18-04-2016 #
############################

##############################
###### VARIABLES TO SET ######
##############################
#
### PATHs to files and folders
#set working directory (If necessary)
setwd("/hps/nobackup/stegle/users/mjbonder/tools/450K_Processing/")
#
# If working on Windows set this:
#memory.limit(200000)
#
# set PATH to results folder
PATH_RES <- "/hps/nobackup/stegle/users/mjbonder/tools/450K_Processing/Res_EPiC/"
#
# set PATH to a folder of "projects" where each project corresponds to a folder of 450K plate extracted data, in .idat format.
PATH_PROJECT_DATA <- "/hps/nobackup/hipsci/scratch/trans_eqtls/IPS_Bulk_Methylation/EPiC/"
#
## set PATH to the file with frequent SNP informations, on which SNP filtering is based. If = NULL, no probe removed. Can handle arrays of filenames.
PATH_ProbeSNP_LIST <- c("./ADDITIONAL_INFO/ProbeFiltering/freq5percent/probeToFilter_450K_1000G_omni2.5.hg19.EUR_alleleFreq5percent_50bp_wInterroSite.txt", "./ADDITIONAL_INFO/ProbeFiltering/ProbesBindingNonOptimal/Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt")
#PATH_ProbeSNP_LIST <- c("./ADDITIONAL_INFO/ProbeFiltering/freq1percent/probeToFilter_450K_GoNL.hg19.ALL_alleleFreq1percent.txt", "./ADDITIONAL_INFO/ProbeFiltering/ProbesBindingNonOptimal/Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt")
#PATH_ProbeSNP_LIST <- c("./ADDITIONAL_INFO/ProbeFiltering/freq1percent/probeToFilter_450K_EnsembleV70&GoNL.hg19.ALL_alleleFreq1percent.txt", "./ADDITIONAL_INFO/ProbeFiltering/AutosomeAndSnpProbes.txt", "./ADDITIONAL_INFO/ProbeFiltering/ProbesBindingNonOptimal/Source&BSProbesMappingMultipleTimesOrNotBothToBSandNormalGenome.txt")
#PATH_ProbeSNP_LIST <- NULL
#
# The name that will be given to result files
projectName = "ILLUMINA450K_IPSc_EPiC"
#
#Set number of cores to use in multi-threading.
coresMultiThread = 4;
#
##############################################
#####  NOW YOU CAN SOURCE THIS SCRIPT !  #####
##############################################
#####
#####
##############################################
###### source scripts and load libraries #####
##############################################
require(lumi)
require(methylumi)
require(minfi)
require(wateRmelon)
#require(RPMM)
#require(preprocessCore)
#require(matrixStats)
#require(IlluminaHumanMethylation450k.db)

source("./SRC/minimalFunctions.R")

#
#
#
{
  if(is.character(PATH_ProbeSNP_LIST)){
    if(length(PATH_ProbeSNP_LIST)==1){
      probeSNP_LIST <- unlist(read.table(file=PATH_ProbeSNP_LIST, quote="", sep="\t", header=TRUE))  
    } else {
      probeSNP_LIST <- NULL
      for(id in 1:length(PATH_ProbeSNP_LIST)){
        probeSNP_LIST <- union(probeSNP_LIST, unlist(read.table(file=PATH_ProbeSNP_LIST[id], quote="", sep="\t", header=TRUE)))
      }
    }
  } 
  else{
    probeSNP_LIST <- NULL
  }
}
#

####################################
# Get project pathes and load data #
####################################
subProjects <- dir(PATH_PROJECT_DATA)

unMeth <- NULL
meth <- NULL

#for all subprojects
for(i in 1:length(subProjects)){
  
  projectName_batch <- subProjects[i]
  sampleTable <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern="TableSample")
  controlTable <- dir(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern="TableControl")
  cat("\n# ")
  cat("Processing sub-project: ", projectName_batch, "\n")
  
  #####
  if(length(sampleTable) < 1 && length(controlTable) < 1 && length(list.files(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern=".idat"))>0){
    
    barcode<- list.files(paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), pattern=".idat")
    
    barcode <- gsub("_Grn","",x=barcode)
    barcode <- gsub("_Red","",x=barcode)
    barcode <- unique(barcode)
    
    if(length(barcode)<2){
      cat("\n\tSkipped folder: ",projectName_batch,"\n")
      cat("\t to little samples")
      next;
    }
    
    cat("\n\tStart data loading...\n")      
    #methLumi_dataTmpData <- methylumIDAT(barcode, idatPath=paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), parallel=T, mc.cores = coresMultiThread, n=T)
    methLumi_dataTmpData <- readEPIC(idatPath=paste(PATH_PROJECT_DATA, projectName_batch, "/", sep=""), barcodes=barcode,force=TRUE)
    #methLumi_dataTmpData <- as(methLumi_dataTmpData, 'MethyLumiM')
    
    cat("\tProject sample nb: ", length(barcode), ".\n")
    cat("\tData dimensions: ", dim(methLumi_dataTmpData)[1],"x", dim(methLumi_dataTmpData)[2], ".\n")
    cat("\t...data loaded..\n\n")
    
    #############################
    # starts data preprocessing #
    #############################
    
    indexProbe2remove <- which(is.element(featureNames(methLumi_dataTmpData), probeSNP_LIST))
    if(length(indexProbe2remove)>0) methLumi_dataTmpData <- methLumi_dataTmpData[-indexProbe2remove,]
    
    if(is.null(methLumi_dataTmpData)){
      next;
    }
    
  } else {
    next;
  }
  
  ################################################
  # Sub-project data & information concatenation #
  ################################################
  
  if(is.null(unMeth) && length(sampleNames(methLumi_dataTmpData))>0){
    unMeth <- unmethylated(methLumi_dataTmpData)
    meth <- methylated(methLumi_dataTmpData)
    
    cat("\t Chanels plate", i, " ok (", dim(unMeth)[1], "x", dim(unMeth)[2], ").\n")
    
    #select "useful" probe annotations
    rm(methLumi_dataTmpData)
    
  } else if(length(sampleNames(methLumi_dataTmpData))>0){
    print(dim(methLumi_dataTmpData))
    #concatenate 'chanels'
    
    if(length(rownames(meth)) == length(rownames(methylated(methLumi_dataTmpData))) && all(rownames(meth) == rownames(methylated(methLumi_dataTmpData)))){
      meth <- cbind(meth, methylated(methLumi_dataTmpData))
      unMeth <- cbind(unMeth, unmethylated(methLumi_dataTmpData))
    } else {
      unMeth_i <- unmethylated(methLumi_dataTmpData)
      meth_i <- methylated(methLumi_dataTmpData)
      
      unMeth_i<- as.matrix(unMeth_i[which(rownames(unMeth_i) %in% rownames(unMeth)),])
      meth_i<- as.matrix(meth_i[which(rownames(meth_i) %in% rownames(meth)),])
      
      unMeth<- as.matrix(unMeth[which(rownames(unMeth) %in% rownames(unMeth_i)),])
      meth<- as.matrix(meth[which(rownames(meth) %in% rownames(meth_i)),])
      
      meth <- cbind(meth[order(rownames(meth)),], meth_i[order(rownames(meth_i)),])
      unMeth <- cbind(unMeth[order(rownames(unMeth)),], unMeth_i[order(rownames(unMeth_i)),])
      
      rm(unMeth_i)
      rm(meth_i)
    }

    cat("\t Chanels ok (", dim(unMeth)[1], "x", dim(unMeth)[2], ").\n")
  }
}

if(is.null(unMeth) || is.null(meth)){
  stop()
}

#split out over: Type I U / Type I M /  Type II U / Type II M

#Epic
epic.ordering <- epic.ordering[which(rownames(epic.ordering)%in%rownames(unMeth)),]
epic.ordering <- epic.ordering[order(rownames(epic.ordering)),]

write.table(unMeth[intersect(which(epic.ordering$DESIGN=="I"),which(epic.ordering$COLOR_CHANNEL=="Red")),], file=paste(PATH_RES, projectName, "_T1_Red_U_Signal.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
write.table(meth[intersect(which(epic.ordering$DESIGN=="I"),which(epic.ordering$COLOR_CHANNEL=="Red")),], file=paste(PATH_RES, projectName, "_T1_Red_M_Signal.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
write.table(unMeth[intersect(which(epic.ordering$DESIGN=="I"),which(epic.ordering$COLOR_CHANNEL=="Grn")),], file=paste(PATH_RES, projectName, "_T1_Grn_U_Signal.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
write.table(meth[intersect(which(epic.ordering$DESIGN=="I"),which(epic.ordering$COLOR_CHANNEL=="Grn")),], file=paste(PATH_RES, projectName, "_T1_Grn_M_Signal.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
write.table(unMeth[which(epic.ordering$DESIGN=="II"),], file=paste(PATH_RES, projectName, "_T2_U_Signal.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
write.table(meth[which(epic.ordering$DESIGN=="II"),], file=paste(PATH_RES, projectName, "_T2_M_Signal.txt", sep=""), quote=FALSE, sep="\t", col.names = NA)
