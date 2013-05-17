

plotQC <- function(beta, figName, PATH="./"){

	plotMethylationDensity(beta, PATH=PATH, figName, xname="beta", xmin=-0.1, xmax=1.1, ymin=0, ymax=6, bw=0.05)

	if(dim(beta)[2]>30) L = 15*dim(beta)[2]
	else{L = 500}
	hclustPlot(t(beta), filename=figName, figName, xlab="beta-value", h=500, l=L, PATH=PATH)
}