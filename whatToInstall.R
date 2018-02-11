source("http://bioconductor.org/biocLite.R")
biocLite()

##Minimal pipeline (only works on R 2 & R 3)
biocLite("matrixStats")
biocLite("wateRmelon")
install.packages("RPMM")

biocLite("minfi")
biocLite("methylumi")
biocLite("lumi")
biocLite("preprocessCore")




##Extended pipeline (only works on R 2)
biocLite("IlluminaHumanMethylation450k.db")

