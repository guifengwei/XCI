
data <- read.table('GeneExpressionTable.Normalized.txt2', header = F, sep='\t')

Sum <- apply(data[,5:12],1,sum) 

DATA <- cbind(data, Sum)


datachr6 <- DATA[which(DATA$V1 == 'chr6' & DATA$Sum > 8),]

datachrAuto <- DATA[which( DATA$V1 != 'chr6' & DATA$V1 != 'chrX' & DATA$V1 != 'chrY' & DATA$Sum > 10), ]

boxplot(datachr6$V5, datachr6$V6, datachr6$V7, datachr6$V8, datachr6$V9, datachr6$V10, datachr6$V11, datachr6$V12, outline=F, col='cyan', ylim=c(0,150))

boxplot(datachrAuto$V5, datachrAuto$V6, datachrAuto$V7, datachrAuto$V8, datachrAuto$V9, datachrAuto$V10, datachrAuto$V11, datachrAuto$V12, 
        outline=F, col='pink', ylim=c(0,150))


boxplot(datachr6$V5, datachrAuto$V5, datachr6$V6, datachrAuto$V6,datachr6$V7, datachrAuto$V7, datachr6$V8, datachrAuto$V8,  datachr6$V9, datachrAuto$V9, datachr6$V10, datachrAuto$V10, datachr6$V11, datachrAuto$V11, datachr6$V12, datachrAuto$V12, outline=F, col=rep(c('cyan', 'pink'), 8) , ylim=c(0,160) ) 
