  library(ggplot2)
  setwd("C:\\Users\\sunzh\\OneDrive\\×ÀÃæ\\STREAM\\FalsePositives\\Final_format\\supFigs\\Fig1")
  ph <- read.table("Theoretical_2b_tags.txt",header=T,sep="\t")
  ph$type <- factor(ph$type,levels = c("CjePI","CjeI","HaeIV","Hin4I","BslFI","BcgI","BsaXI","Bsp24I","AlfI","FalI","BplI","BaeI","AloI","PpiI","PsrI","CspCI"))
  ggplot(ph,aes(x=type,y=log10(Theoretical)))+geom_violin(aes(fill=type))+
    geom_boxplot(width = 0.2,outlier.colour = NA)+theme_bw()+
    scale_y_continuous(breaks = c(-1,0,1,2,3,4,5,6))
  ggsave("Theoretical.pdf",width = 12, height = 3)
