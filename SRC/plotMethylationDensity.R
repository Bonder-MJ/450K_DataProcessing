plotMethylationDensity <-function(data, PATH, figName, xname, xmin=-0.1, xmax=1.1, ymin=0, ymax=NULL, LEGEND=FALSE, bw="nrd0", grayScale=TRUE, color=NULL){

	xlegend = xmax - (xmax/3)*2
	
	if(is.matrix(data)){
		nbSamples <- dim(data)[2]
		sampleNames <- colnames(data)
	}
	if(is.list(data)){
		nbSamples <- dim(data)[2]
		sampleNames <- colnames(data)
	}
	#print(sampleNames)

	if(nbSamples > 10) sizeLegend <- 0.8
	else{sizeLegend = 3}

	samples.sort <- sort(sampleNames, index.return=TRUE)

	data <- data[,samples.sort$ix]

	if(is.null(color)){
		if(grayScale) color <- gray(seq(0.7,0, length=nbSamples))
		else{ color <- c("blue", "red", "green", "black", "orange", "pink", "cyan", "purple", "grey", "brown", "yellow", "violet")}
	}
	plotType = "l"
	plotSize <- c(3.5, 4, 4.5, 5, 5.5)

	jpeg(file = paste(PATH, figName, ".jpeg", sep=""),width=1700, height=900, units="px", pointsize=16, quality=100, bg="white")

	ic <- 1
	it <- 1
	legendColor <- vector()
	legendPlotSize <- vector()


	if(is.null(ymax)){
	for(i in 1:nbSamples){
		d <- density(as.numeric(data[,i]), na.rm = TRUE, bw=bw)
		if(i==1) {
			ymax <- max(d$y)
		}
		else { if(ymax < max(d$y)) ymax <- max(d$y)}
	}
	}

	for(i in 1:nbSamples){
		#print(paste("Process of sample: ", i, sep=""))
		d <- density(as.numeric(data[,i]), na.rm = TRUE, bw=bw)

		legendColor <- c(legendColor, as.character(color[ic]))
		legendPlotSize <- c(legendPlotSize, plotSize[it])
		
		if(i==1) {
			plot(d, col=color[ic], type=plotType, main = "", xlab=xname, ylab="", xlim=c(xmin,xmax), ylim=c(0,ymax), lwd=plotSize[it], cex.axis=2.5, cex.lab=3)
		}
		else{ 
			lines(d, col=color[ic], type=plotType, lwd=plotSize[it])
		}
		
		if(ic == length(color)) { ic <- 1 ; it <- it + 1}
		else{ ic <- ic+1 }
		
		if(it > length(plotSize)) it <- 1

	}
	if(LEGEND) legend(x = xlegend, y = ymax, samples.sort$x, cex=sizeLegend, col=legendColor, bty="n", lwd=legendPlotSize, lty=1)

	dev.off()
}
