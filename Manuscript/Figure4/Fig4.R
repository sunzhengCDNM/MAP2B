#=============================  Figure4 Corrplot  =========================================
#install.packages("corrplot")
library(corrplot)
library(RColorBrewer)
setwd("C:\\Users\\Think\\Desktop\\aDist\\Figure_AutoGscore\\Figure2-Corrplot")
##  MOCK
svm <- read.table("plot_Corrplot_MOCK_PRF.Dist_MPA4.mOTUs3-23.0428.txt", header = T,row.names = 1,sep="\t")
#shape4 label
corrplot(as.matrix(svm), is.corr = FALSE, method = 'color',
         col=colorRampPalette(brewer.pal(11, "Spectral"))(100),outline="white",
         addCoef.col = 'grey80',
         addgrid.col="grey80", tl.col="black")





#=============================Figure MOCK  =========================================
library(plot3D)
library(viridis)
setwd("C:\\Users\\Think\\Desktop\\aDist\\Figure_AutoGscore\\Figure4-MOCK 3D")
feature_table <- read.table("pred_MOCK_1002.tag.xls",header = T, row.names = 1, sep="\t")
#change color
library(RColorBrewer)
scatter3D(-feature_table$Tag_Depth_log_norm, 
          -feature_table$Theoretical_Reads_DTR_log_norm,
          feature_table$Coverage_log_norm, bty="b2", 
          col=colorRampPalette(brewer.pal(11, "Spectral"))(100), 
          colvar = feature_table$G_Score_log_norm, phi =0, 
          theta=50, pch = (feature_table$Tag*12+4), #, theta=50
          cex = 1.2, ticktype = "detailed")