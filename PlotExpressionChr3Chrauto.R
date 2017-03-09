
data <- read.table('GeneExpressionTable.Normalized.txt2', header = F, sep='\t')

Sum <- apply(data[,5:10],1,sum) 

DATA <- cbind(data, Sum)


datachr3 <- DATA[which(DATA$V1 == 'chr3' & DATA$Sum > 8),]

datachrAuto <- DATA[which( DATA$V1 != 'chr3' & DATA$V1 != 'chrX' & DATA$V1 != 'chrY' & DATA$Sum > 10), ]

boxplot(datachr3$V5, datachr3$V6, datachr3$V7, datachr3$V8, datachr3$V9, datachr3$V10,  outline=F, col='cyan', ylim=c(0,150))

boxplot(datachrAuto$V5, datachrAuto$V6, datachrAuto$V7, datachrAuto$V8, datachrAuto$V9, datachrAuto$V10, 
        outline=F, col='pink', ylim=c(0,150))


boxplot(datachr3$V5, datachrAuto$V5, datachr3$V6, datachrAuto$V6,datachr3$V7, datachrAuto$V7, datachr3$V8, datachrAuto$V8,  datachr3$V9, datachrAuto$V9, datachr3$V10, datachrAuto$V10, outline=F, col=rep(c('cyan', 'pink'), 6) , ylim=c(0,150) ) 
