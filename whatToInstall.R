source("http://bioconductor.org/biocLite.R")
biocLite()

##Minimal pipeline (only works on R 2 & R 3)
biocLite("lumi")
biocLite("methylumi")

##Extended pipeline (only works on R 2)
biocLite("preprocessCore")
biocLite("minfi")
biocLite("matrixStats")
biocLite("IlluminaHumanMethylation450k.db")
install.packages("RPMM")
