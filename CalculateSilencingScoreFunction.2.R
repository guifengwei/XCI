
Calculate_Silencing_Score <- function(data, min_exp = 4)
{
	nrow = dim(data)[1]
	data = data + 0.0001
	data$Silencing = 0
	
	for(i in 1:nrow)
	{
#		
#		else if(data[i,2] + data[i,1] < min_exp & data[i,3] + data[i,4] < min_exp){
#			data[i,5] = 0 }
#		
#		else if(data[i,2] - data[i,1] > 0 & data[i,4] - data[i,3] > 0){
#			if( (data[i,2] - data[i,1])/data[i,1] >= (data[i,4] - data[i,3])/data[i,3]){
#				data[i,5] = 0}
#			else{
#				data[i,5] = 1 - ((data[i,4] - data[i,3])/data[i,3] - (data[i,2] - data[i,1])/data[i,1]) 
#				data[i,5] <- ifelse( data[i,5] <0 , 0, data[i,5])
#				}
#			}
#		
#		else if(data[i,4] - data[i,3] > 0 & data[i,2] - data[i,1] <= 0){
#			data[i,5] = 1
#			}
#		
#		else if(data[i,4] - data[i,3] <=0 & data[i,2] - data[i,1] > 0){
#			data[i,5] = 0}
#		
#		else if(data[i,4] - data[i,3] <0 & data[i,2] - data[i,1] < 0){
#			if(abs( data[i,2] - data[i,1])/data[i,1] > abs( data[i,4] - data[i,3])/data[i,3]){
#				data[i,5] = abs( data[i,2] - data[i,1])/data[i,1] - abs( data[i,4] - data[i,3])/data[i,3] }
#			else{
#				data[i,5] = 0}
#			}
#		}
        
################# Modify the Methods for calculating the silencing score
			d_score_NoDox = data[i,1]/(data[i,1] + data[i,3]) - 0.5
			d_score_Dox   = data[i,2]/(data[i,2] + data[i,4]) - 0.5
			
			data[i,5] = d_score_Dox - d_score_NoDox
			
			##if(data[i,5]>1){data[i,5] = 1}
			#else if(data[i,5]<0){data[i,5] = 0}
	}
	data$Silencing
}

