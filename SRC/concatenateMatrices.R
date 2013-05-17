# Concatenate 2 matrices, by rows in order to increase column number
# If matrices have different row names (or different row numbers) it limits the concatenation to common and sorted (by names) rows.
#
# Args: 
#	- matrix1: matrix
#	- matrix2: matrix (may have different sizes and different row names than matrix1)

concatenateMatrices <- function(matrix1, matrix2){
		
		commonRows <- intersect(rownames(matrix1), rownames(matrix2))
		if(length(commonRows) > 0){
			matrix1 <- matrix1[which(is.element(rownames(matrix1), commonRows)), , drop=FALSE]
			matrix2 <- matrix2[which(is.element(rownames(matrix2), commonRows)), , drop=FALSE]
		
			matrix1 <- matrix1[sort(rownames(matrix1), index.return=TRUE)$ix, , drop=FALSE]
			matrix2 <- matrix2[sort(rownames(matrix2), index.return=TRUE)$ix, , drop=FALSE]
		
			return(cbind(matrix1, matrix2))
		}
		else{
			print("WARNING ! in 'concatenateMatrices', matrix1 and matrix2 must share at least a common row name.")
			return("WARNING ! in 'concatenateMatrices', matrix1 and matrix2 must share at least a common row name")
		}
}
