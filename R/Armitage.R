Armitage <-
function(genotypeStat){
     if(ncol(genotypeStat)==3 && nrow(genotypeStat)==2){
	     sum12=sum(genotypeStat)
         sum1=sum(genotypeStat[1,]) 
         sum2=sum(genotypeStat[2,])
         temp1=2*genotypeStat[1,1]+genotypeStat[1,2]
         temp2=2*genotypeStat[1,1]+genotypeStat[1,2]+2*genotypeStat[2,1]+genotypeStat[2,2]  
         temp3=4*genotypeStat[1,1]+genotypeStat[1,2]+4*genotypeStat[2,1]+genotypeStat[2,2]
         upp=sum12*(sum12*temp1-sum1*temp2)*(sum12*temp1-sum1*temp2)
         low=sum1*sum2*(sum12*temp3-temp2*temp2)
         statistic=upp/low
         pValue=pchisq(statistic,df=1,lower.tail=FALSE)
         list(statistic=statistic,pValue=pValue)
	 }else{
	     list(statistic=NA,pValue=NA)
	 }
}
