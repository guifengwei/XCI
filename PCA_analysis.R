




data <- read.table("Chr11_CellLines_comm2_RS.txt", header=T, sep="\t")
library(FactoMineR)
x<-c(5:14,17,18,19,20)
data_f <- data[,x]
data_pca <- PCA(t(data_f[330:1010,]), quanti.sup=1:14)
plot(data_pca)
