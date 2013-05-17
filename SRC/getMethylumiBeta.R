require(lumi)
require(methylumi)

#compute beta values from methylated & unmethylated signals
#
# Args: methylumi object (see lumi & methylumi packages)
#
# Returns: a n rows x m columns matrix of beta values with n as probes and m as samples
#

getMethylumiBeta <- function(methylumi, alfa=100){

	fNames <- featureNames(methylumi)
	u <- unmethylated(methylumi)
	m <- methylated(methylumi)

	rm(methylumi)

	#check and "correct" for negative values
	indexNegU <- which(u < 0, arr.ind=TRUE)
	indexNegM <- which(m < 0, arr.ind=TRUE)
	u[indexNegU] <- 0
	m[indexNegM] <- 0

	rm(indexNegU, indexNegM)

	beta <- m/(m + u + alfa)
	rm(u, m)

	rownames(beta) <- fNames

	return(beta)
}
