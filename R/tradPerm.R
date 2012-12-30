tradPerm <-
function(genotypeLine,affectionLine,fromCol,naString,sep,testMethod="chisq",repeatNum=1000){
     if(abs(repeatNum-round(repeatNum)) > 1e-10){
	     stop("'repeatNum' must be a integer.")
	 }

	 ## observe value calculate
     obsStat=genotypeStat(genotypeLine=genotypeLine,affectionLine=affectionLine,fromCol=fromCol,naString=naString,sep=sep)
	 obsGenotypeStat=obsStat$genotypeStat
	 colNum=ncol(obsGenotypeStat)
	 colName=colnames(obsGenotypeStat)
	 if(any(colName==naString)){
	     obsGenotypeStat=obsGenotypeStat[,-colNum]
	 }
	 
     switch(testMethod,
         chisq={
		     chiValue=matrix(0,nrow=1,ncol=repeatNum+1)
			 chiP=matrix(0,nrow=1,ncol=repeatNum+1)
	         temp=chisq.test(obsGenotypeStat)
			 chiValue[1,1]=temp$statistic
			 chiP[1,1]=temp$p.value
	     },
         fisher={
		     fisherP=matrix(0,nrow=1,ncol=repeatNum+1)
	         fisherP[1,1]=fisher.test(obsGenotypeStat)$p.value
	     },
		 OR={
		     alleleStat=obsStat$alleleStat
		     ORvalue=matrix(0,nrow=1,ncol=repeatNum+1)
			# CIvalue=matrix(0,nrow=1,ncol=repeatNum+1)
			 temp=OR(alleleStat[1,1],alleleStat[1,2],alleleStat[2,1],alleleStat[2,2])
			 ORvalue[1,1]=temp$OR
			# CIvalue[1,1]=temp$CI[1,2]-temp$CI[1,1]
		 },
		 Armitage={
		     trendValue=matrix(0,nrow=1,ncol=repeatNum+1)
			 trendP=matrix(0,nrow=1,ncol=repeatNum+1)
		     temp=Armitage(obsGenotypeStat)
		     trendValue[1,1]=temp$statistic
			 trendP[1,1]=temp$pValue
		 }
     )
	 
	 ## null distribution calculate
	 for(i in 2:(repeatNum+1)){
	     randData=permuteData(dataLine=genotypeLine,fromCol=fromCol)
	     expStat=genotypeStat(genotypeLine=randData,affectionLine=affectionLine,fromCol=fromCol,naString=naString,sep=sep)
	     expGenotypeStat=expStat$genotypeStat
	     colName=colnames(expGenotypeStat)
	     colNum=ncol(expGenotypeStat)
	     if(any(colName==naString)){
	         expGenotypeStat=expGenotypeStat[,-colNum]
	     }
		 
		 switch(testMethod,
             chisq={
	             temp=chisq.test(expGenotypeStat)
			     chiValue[1,i]=temp$statistic
			     chiP[1,i]=temp$p.value
	         },
             fisher={
	             fisherP[1,i]=fisher.test(expGenotypeStat)$p.value
	         },
			 OR={
		         expAlleleStat=expStat$alleleStat
				 temp=OR(expAlleleStat[1,1],expAlleleStat[1,2],expAlleleStat[2,1],expAlleleStat[2,2])
				 ORvalue[1,i]=temp$OR
			    # CIvalue[1,i]=temp$CI[1,2]-temp$CI[1,1]
		     },
		     Armitage={
		         temp=Armitage(expGenotypeStat)
		         trendValue[1,i]=temp$statistic
			     trendP[1,i]=temp$pValue
		     }
         )	 
	 }

	 ## calculate p value and return value
	 switch(testMethod,
         chisq={
	         obs=chiValue[1,1]
			 maxExp=as.matrix(chiValue[chiValue>obs])
			 pValue=nrow(maxExp)/repeatNum
			 list(pValue=pValue,obsStatistic=obs,obsP=chiP[1,1],permStatistic=chiValue[,2:(repeatNum+1)],permP=chiP[,2:(repeatNum+1)])
	     },
         fisher={
	         obs=fisherP[1,1]
			 minExp=as.matrix(fisherP[fisherP<obs])
			 pValue=nrow(minExp)/repeatNum
			 list(pValue=pValue,obsP=obs,permP=fisherP[,2:(repeatNum+1)])
	     },
	     OR={
	         obsOR=ORvalue[1,1]
			 maxExp=as.matrix(ORvalue[ORvalue>obsOR])
			 ORp=nrow(maxExp)/repeatNum
			# obsCI=CIvalue[1,1]
			# minExp=as.matrix(CIvalue[CIvalue<obsCI])
			# CIp=nrow(minExp)/repeatNum
			# list(ORp=ORp,CIp=CIp,obsOR=obsOR,obsCI=obsCI,permOR=ORvalue[,2:(repeatNum+1)],permCI=CIvalue[,2:(repeatNum+1)])
			 list(ORp=ORp,obsOR=obsOR,permOR=ORvalue[,2:(repeatNum+1)])
	     },
		 Armitage={
		     obs=trendValue[1,1]
			 maxExp=as.matrix(trendValue[trendValue>obs])
		     pValue=nrow(maxExp)/repeatNum
			 list(pValue=pValue,obsStatistic=obs,obsP=trendP[1,1],permStatistic=trendValue[,2:(repeatNum+1)],permP=trendP[,2:(repeatNum+1)])

		 }
     )
}
