
Calculate_Normed_Split_Exp<- function(raw_count, normed_count){
	g1=sum(rbinom(raw_count, 1, 0.5))
	g2=raw_count-g1
	normed_g1 = normed_count * g1/(g1+g2)
	normed_g2 = normed_count * g2/(g1+g2)
	#### return normed_g1 and normed_g2
	c(normed_g1, normed_g2)
  }


FishersMethod = function(x){
	pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
}

Calculate_RepressionScore <- function(d1,d2,d3,d4){
	#data[i,5] = (data[i,4] - data[i,3])/data[i,3] - (data[i,2] - data[i,1])/data[i,1]
	A <- c(d1,d2,d3,d4) + 0.0001
	RS = (A[4] - A[3])/A[3] - (A[2] - A[1])/A[1]
	RS
	}

# RawCounts_NoDox, RawCounts_Dox1, RawCounts_Dox2, NormedCounts_NoDox, NormedCounts_Dox1, NormedCounts_Dox2, RS1, RS2
Permutation<- function(r1,r2,r3, e1,e2,e3, rs1, rs2){
	Permutation_N = 10000
	rs_p1 <- rep(0, Permutation_N)
	rs_p2 <- rep(0, Permutation_N)
    ######
	for(i in 1:Permutation_N){
		no_dox <- Calculate_Normed_Split_Exp(r1,e1)
		dox1   <- Calculate_Normed_Split_Exp(r2,e2)
		dox2   <- Calculate_Normed_Split_Exp(r3,e3)

		rs_p1[i] <- Calculate_RepressionScore(no_dox[1], dox1[1], no_dox[2], dox1[2])
		rs_p2[i] <- Calculate_RepressionScore(no_dox[1], dox2[1], no_dox[2], dox2[2])

		}
	p1 = sum(ifelse(rs_p1>=rs1,1,0))/Permutation_N
	p2 = sum(ifelse(rs_p2>=rs2,1,0))/Permutation_N
	c(p1,p2, FishersMethod(c(p1,p2)))
}


# for(i in 1:659){T[i,13:15]<-Permutation(T[i,5],T[i,6],T[i,7],T[i,8],T[i,9],T[i,10],T[i,11],T[i,12])} 
# T$qValue <- p.adjust(T$FisherP, method = 'BH', n = length(T$FisherP))

Data <- read.table("RawCounts_NormedCounts_NormedAllelicExpression_Chr2.txt2", header=T, sep='\t')
Data$p1 = 0
Data$p2 = 0
Data$FisherP = 0
n = dim(Data)[1]
for(i in 1:n){
	print(paste(round(i/n,4)*100, "%", sep=''))
	Data[i,21:23]<-Permutation(Data[i,5],Data[i,6], Data[i,7], Data[i,8], Data[i,9], Data[i,10], Data[i,18], Data[i,19])
	}

Data$qValue <- p.adjust(Data$FisherP, method = 'BH', n = length(Data$FisherP))

Data$Status <- ifelse(Data$qValue<0.05, 'Silenced', 'Escaped')

write.table(Data, 'RawCounts_NormedCounts_NormedAllelicExpression_Chr2.Pvalue.txt2', sep='\t', quote=F, col.names=T, row.names=F)



