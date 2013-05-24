#beadcount function, creates matrix with NAs representing probes with beadcount <3 from Extended RG Channel Set
beadcountMJ<-function(x){
  #select out bead count dataframe
  assayDataElement(x, "NBeads")->nb
  
  #match rownames of beadcount dataframe to addresses
  getProbeInfo(x,type="I")->typeIadd
  match(typeIadd$AddressA,rownames(nb))->typeImatchA
  match(typeIadd$AddressB,rownames(nb))->typeImatchB
  
  #match rownames of beadcount dataframe to addresses
  getProbeInfo(x,type="II")->typeIIadd
  match(typeIIadd$Address,rownames(nb))->typeIImatch
  
  nb->nbcg
  
  locusNames <- getManifestInfo(x, "locusNames")
  bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                    dimnames = list(locusNames, sampleNames(x)))
  
  TypeII.Name <- getProbeInfo(x, type = "II")$Name
  bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$Address,]
  
  TypeI <- getProbeInfo(x, type = "I")
  
  bc_temp->bcB
  bc_temp->bcA
  
  bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB,]
  bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA,]
  
  which(bcB<0)->bcB3
  which(bcA<0)->bcA3
  bcA->bcA2
  bcB->bcB2
  
  bcA2[bcA3]<-NA
  bcB2[bcB3]<-NA
  bcA2[bcB3]<-NA
  bcB2[bcA3]<-NA
  
  data.frame(bcA2)->bcM
  data.frame(bcA2)->bcU
  data.frame(bcA2)->bc
  bc
}