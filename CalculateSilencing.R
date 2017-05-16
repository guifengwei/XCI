
DATA <- read.table('RawCounts_NormedCounts_NormedAllelicExpression.txt2', header = F, sep='\t')
Sum <- apply(DATA[,11:16],1,sum)  
DATA <- cbind(DATA, Sum)  

colnames(DATA) <- c('Chromosome', 'TSSStart', 'TSSStop', 'GeneName', 'RawCounts_NoDox', 'RawCounts_DoxA', 'RawCounts_DoxB', 'NormedCounts_NoDox', 'NormedCounts_DoxA', 'NormedCounts_DoxB', 'Genome1_NoDox', 'Genome1_DoxA', 'Genome1_DoxB', 'Genome2_NoDox', 'Genome2_DoxA', 'Genome2_DoxB', 'Sum')

x <- c(5,6,8,9) + 6 
y <- c(5,7,8,10) + 6

Cmp_data_x <- DATA[,x] 
Cmp_data_y <- DATA[,y]

source('CalculateSilencingScoreFunction.R')

Cmp_data_ss_x <- Calculate_Silencing_Score(Cmp_data_x)
Cmp_data_ss_y <- Calculate_Silencing_Score(Cmp_data_y)

transcriptome_ss <- cbind(DATA, Cmp_data_ss_x, Cmp_data_ss_y)

for(chrs in c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')){
	print(chrs)
	print(median(transcriptome_ss[transcriptome_ss$Chromosome==chrs,]$Cmp_data_ss_x))
	print(median(transcriptome_ss[transcriptome_ss$Chromosome==chrs,]$Cmp_data_ss_y))
}

D <- transcriptome_ss[transcriptome_ss$Chromosome =="chr7", ]

D$site = ceiling(D$TSSStart/1000000)
library(ggplot2)
p1<-ggplot(D, aes(x=factor(site), y= Cmp_data_ss_x) ) + geom_boxplot(outlier.shape = 3) + geom_point(position=position_jitter(width=0.25,height=0.001),shape=16, colour="purple", alpha=.55, size=2)

p2<-ggplot(D, aes(x=factor(site), y= Cmp_data_ss_y) ) + geom_boxplot(outlier.shape = 3) + geom_point(position=position_jitter(width=0.25,height=0.001),shape=16, colour="purple", alpha=.55, size=2)

library(gridExtra)
grid.arrange(p1,p2,p1,nrow=3, ncol=1)

for(i in 1:max(D$site)){
	print(c(i, median(D[D$site==i,]$Cmp_data_ss_x),  median(D[D$site==i,]$Cmp_data_ss_y)), quote=F)
	}

# sum(ifelse(D[which(D$site>120 & D$site<131),]$Cmp_data_ss_x>0.75,1,0))
# sum(ifelse(D[which(D$site>120 & D$site<131),]$Cmp_data_ss_x>=0.75,1,0))
# sum(ifelse(D[which(D$site>120 & D$site<131),]$Cmp_data_ss_x>=0.5,1,0)) 
# sum(ifelse(D[which(D$site>120 & D$site<131),]$Cmp_data_ss_x>=0.25,1,0))
