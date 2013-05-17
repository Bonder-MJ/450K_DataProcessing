# Nizar TOULEIMAT
# nizar.touleimat @ cng.com
#
# This function performs a classical sample quantile normalization.
#
# Args:
#	- X: data matrix with variables in rows and samples in columns
#	- Reference.Quantiles: vector of reference quantiles of length = dim(data)[1]
#
# Returns: a matrix of normalized data with same dimension as the original data matrix.
#
normalize.quantiles2 <- function(X, Reference.Quantiles){

    apply(X, 2, function(x, y) y[rank(x)], Reference.Quantiles)
}



