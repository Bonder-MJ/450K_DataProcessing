##################################
#                                #
#  Dasen is extracted from the   #
#       Watermelon package       #
#                                #
#  Procedure by Pidsley et al.   #
#                                #
##################################

dasen <- function(mns, uns, onetwo, alfa=100, MvalueConv=TRUE, ...){
   
   mnsc <- dfsfit(mns,  onetwo, ...)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   
   mnsc[onetwo=='I' ,] <- normalize.quantiles(mnsc[onetwo=='I', ])
   unsc[onetwo=='I' ,] <- normalize.quantiles(unsc[onetwo=='I', ])
   
   mnsc[onetwo=='II',] <- normalize.quantiles(mnsc[onetwo=='II',])
   unsc[onetwo=='II',] <- normalize.quantiles(unsc[onetwo=='II',])
   
   indexNegU <- which(!is.numeric(unsc), arr.ind=TRUE)
   indexNegM <- which(!is.numeric(mnsc), arr.ind=TRUE)
   if(length(indexNegU)>0 || length(indexNegM)>0){
     cat("\tWarning: NA values introduced. Value is replaced by 0.\n")
     
     unsc[indexNegU] <- 0
     mnsc[indexNegM] <- 0
   }
   
   #check and "correct" for negative values
   indexNegU <- which(unsc < 0, arr.ind=TRUE)
   indexNegM <- which(mnsc < 0, arr.ind=TRUE)
   unsc[indexNegU] <- 0
   mnsc[indexNegM] <- 0
   
   if(MvalueConv){
     return(log2((mnsc+alfa)/(unsc + alfa)))
   } else{
     return(mnsc/( mnsc + unsc + alfa ))
   }

}

dfsfit <- function(
  mn, 
  onetwo,
  roco=unlist(data.frame(strsplit(colnames(mn), '_'), stringsAsFactors=FALSE)[2,])   
  ){
    mdf<-apply(mn,2,dfs2,onetwo)
        
    if (! is.null(roco) ) {
      scol  <- as.numeric(substr(roco,6,6))
      srow  <- as.numeric(substr(roco,3,3))
      fit   <- try(  lm( mdf ~ srow + scol ), silent=TRUE) 
      if (! inherits (fit, "try-error") ) {mdf   <- fit$fitted.values}
      else { message ('Sentrix position model failed, skipping') }
    }
    
    otcor <-  matrix(
      rep(mdf, sum(onetwo=='I')),
      byrow=T, 
      nrow=sum(onetwo=='I')
    )
    
    mn[onetwo=='I',] <- mn[onetwo=='I',] - otcor
    mn
  }

# x is a matrix of intensities
# onetwo is a character vector 
# of same order and length 
# indicating assay I or II

dfs2 <- function(x, onetwo){ 
    one <- density(x[onetwo=='I'], na.rm=T, n = 2^15, from = 0, to = 5000)
    two <- density(x[onetwo=='II'],na.rm=T, n = 2^15, from = 0, to = 5000)
    one$x[which.max(one$y)] - two$x[which.max(two$y)]
  }

