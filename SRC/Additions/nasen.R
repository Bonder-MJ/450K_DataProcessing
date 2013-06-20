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
  
  if(MvalueConv){
    return(log2((mns+alfa)/(uns + alfa)))
  } else{
    return(mns/( mns + uns + alfa ))
  }
  
}