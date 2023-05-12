#=============================  Fig.S3 Corrplot  ================================
##  Figure Paper
#install.packages("corrplot")
library(corrplot)
library(RColorBrewer)
#setwd("C:\\Users\\Think\\Desktop\\aDist\\Figure_AutoGscore\\Figure2-Corrplot")
svm <- read.table("plot_Corrplot_Paper_PRF.L2.BC-Tax.Seq-30-F1_23.0425.txt", header = T,row.names = 1,sep="\t")
#shape3 just color
corrplot(as.matrix(svm), is.corr = FALSE, method = 'color',	#method: "pie", "square"
         col=colorRampPalette(brewer.pal(11, "Spectral"))(100),outline="white", 
         #addCoef.col = 'grey80',
         addgrid.col="grey80", tl.col="black")