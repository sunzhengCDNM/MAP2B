#=============================Figure1d-f 3D scatter plot==============================
#install.packages("plot3D")
library(plot3D)
library(viridis)
setwd("C:\\Users\\Think\\Desktop")
feature_table <- read.table("pred_marmgCAMI2_short_read_sample_0.tag.xls",header = T, row.names = 1, sep="\t") 
feature_table <- read.table("pred_rhimgCAMI2_sample_0.tag.xls",header = T, row.names = 1, sep="\t") 
feature_table <- read.table("pred_strmgCAMI2_short_read_sample_0.tag.xls",header = T, row.names = 1, sep="\t")
scatter3D(-feature_table$Tag_Depth_log_norm, 
          -feature_table$Theoretical_Reads_DTR_log_norm,
          feature_table$Coverage_log_norm, 
          bty="b2", colvar = feature_table$G_Score_log_norm, 
          phi =0, pch = feature_table$Tag+15, cex = 1.2, 
          ticktype = "detailed")