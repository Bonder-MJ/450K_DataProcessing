##################################
#                                #
#  Nasen is extracted from the   #
#       Watermelon package       #
#                                #
#  Procedure by Pidsley et al.   #
#                                #
##################################

nasen <- function(mns, uns, onetwo, alfa=100, MvalueConv=TRUE, ...){
  
  mns[onetwo=='I' ,] <- normalize.quantiles(mns[onetwo=='I', ])
  uns[onetwo=='I' ,] <- normalize.quantiles(uns[onetwo=='I', ])
  
  mns[onetwo=='II',] <- normalize.quantiles(mns[onetwo=='II',])
  uns[onetwo=='II',] <- normalize.quantiles(uns[onetwo=='II',])
  
  
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
  
  
  if(MvalueConv){
    return(log2((mns+alfa)/(uns + alfa)))
  } else{
    return(mns/( mns + uns + alfa ))
  }
  
}