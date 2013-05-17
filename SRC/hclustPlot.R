hclustPlot <- function(data, filename, figName, xlab, h=500, l=1000, PATH="./"){

	data.dist <- dist(data, method = "euclidean")

	fit <- hclust(data.dist, method="ward")

	jpeg(paste(PATH, "hclust_", filename, ".jpeg", sep=""), height=h, width=l)
	plot(fit, main=figName, xlab=xlab)
	dev.off()
}
