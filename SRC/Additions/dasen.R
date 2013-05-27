dasen <-
function(mns, uns, onetwo, alfa=100, ...){
   mnsc <- dfsfit(mns,  onetwo, ...)  
   unsc <- dfsfit(uns,  onetwo, roco=NULL)
   mnsc[onetwo=='I' ,] <- normalize.quantiles(mnsc[onetwo=='I', ])
   mnsc[onetwo=='II',] <- normalize.quantiles(mnsc[onetwo=='II',])

   unsc[onetwo=='I' ,] <- normalize.quantiles(unsc[onetwo=='I', ])
   unsc[onetwo=='II',] <- normalize.quantiles(unsc[onetwo=='II',])

   mnsc/( mnsc + unsc + alfa )
}


dfsfit <-
function(
   mn, 
   onetwo,
   roco=unlist(
      data.frame(
         strsplit(
            colnames(mn), 
            '_'
         ), 
         stringsAsFactors=FALSE
      )[2,] 
   )
){

   mdf<-apply(mn,2,dfs2,onetwo)

   if (! is.null(roco) ) {
      scol  <- as.numeric(substr(roco,6,6))
      srow  <- as.numeric(substr(roco,3,3))
      fit   <- lm( mdf ~ srow + scol )
      mdf   <- fit$fitted.values
   }
   otcor <-  matrix(
      rep( mdf, sum(onetwo=='I')),
      byrow=T, 
      nrow=sum(onetwo=='I')
   )
   mn[onetwo=='I',] <- mn[onetwo=='I',] - otcor
   mn
}
