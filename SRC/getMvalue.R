getBetaMj <- function(u, m, alfa=100){
  
  fNames <- rownames(u)
  
  #check and "correct" for na values
  
  indexNegU <- which(is.na(u), arr.ind=TRUE)
  indexNegM <- which(is.na(m), arr.ind=TRUE)
  if(length(indexNegU)>0 || length(indexNegM)>0){
    cat("\tWarning: NA values introduced. Value is replaced by 0.\n")
    
    u[indexNegU] <- 0
    m[indexNegM] <- 0
  }
  
  
  #check and "correct" for negative values
  indexNegU <- which(u < 0, arr.ind=TRUE)
  indexNegM <- which(m < 0, arr.ind=TRUE)
  
  if(length(indexNegU)>0 || length(indexNegM)>0){
    cat("\tWarning: negative values introduced. Value is replaced by 0.\n")
    
    u[indexNegU] <- 0
    m[indexNegM] <- 0
  }
  
  
  #check and "correct" for non-numeric values
  
  indexNegU <- which(!is.numeric(u), arr.ind=TRUE)
  indexNegM <- which(!is.numeric(m), arr.ind=TRUE)
  if(length(indexNegU)>0 || length(indexNegM)>0){
    cat("\tWarning: non-numeric values introduced. Value is replaced by 0.\n")
    
    u[indexNegU] <- 0
    m[indexNegM] <- 0
  }
  
  rm(indexNegU, indexNegM)
  beta <- m/(m + u + alfa)
  rm(u, m)
  
  rownames(beta) <- fNames
  return(beta)
}
