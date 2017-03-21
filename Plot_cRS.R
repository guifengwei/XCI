


data <- read.table("Chr11_CellLines_comm2_RS_Pvalue.txt", header = TRUE, sep = "\t")


data$site = ceiling(data$TSSStart/1000000)

#for(i in 1:max(data$site)){
#	print(c(i, median(data[data$site==i,]$MiniXist_cRS),  median(data[data$site==i,]$DelBMiniXist_cRS)), quote=F)
#}

library(ggplot2)
data$site2=ifelse(data$site<=61,"L", "R")
p<-ggplot(data, aes(x=log2(MiniXist_cRS+0.001), y=log2(DelB_MiniXist_cRS+0.001), size = Sig, colour=site2)) +geom_point(alpha=0.65) + scale_color_manual('Position', values=c("red","cornflowerblue")) + geom_abline(intercept=0,slope=1)

