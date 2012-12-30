MCPerm <-
function(genotypeStat,testMethod="chisq",repeatNum=1000){
     if(abs(repeatNum-round(repeatNum)) > 1e-10){
	     stop("'repeatNum' must be a integer.")
	 }
	 
	 genotypeCount=genotypeStat
	 delIndex=c()
	 for(n in 1:3){
	     if(all(genotypeCount[,n]==0)){
		     delIndex=c(delIndex,n)
		 }
	 }
	 if(!is.null(delIndex)){
	     genotypeCount=genotypeCount[,-delIndex,drop=FALSE]
	 }
	 
	 ## calculate true statistic value
	 switch(testMethod,
         chisq={
		     chiValue=matrix(0,nrow=1,ncol=repeatNum+1)
			 chiP=matrix(0,nrow=1,ncol=repeatNum+1)
	         temp=chisq.test(genotypeCount)
			 chiValue[1,1]=temp$statistic
			 chiP[1,1]=temp$p.value
	     },
         fisher={
		     fisherP=matrix(0,nrow=1,ncol=repeatNum+1)
	         fisherP[1,1]=fisher.test(genotypeCount)$p.value
	     },
	     OR={
		     caseAllele1=2*genotypeStat[1,1]+genotypeStat[1,2]
			 caseAllele2=2*genotypeStat[1,3]+genotypeStat[1,2]
			 controlAllele1=2*genotypeStat[2,1]+genotypeStat[2,2]
			 controlAllele2=2*genotypeStat[2,3]+genotypeStat[2,2]
			 ORvalue=matrix(0,nrow=1,ncol=repeatNum+1)
			# CIvalue=matrix(0,nrow=1,ncol=repeatNum+1)
			 temp=OR(caseAllele1,caseAllele2,controlAllele1,controlAllele2)
			 ORvalue[1,1]=temp$OR
			# CIvalue[1,1]=temp$CI[1,2]-temp$CI[1,1]
		 },
		 Armitage={
		     trendValue=matrix(0,nrow=1,ncol=repeatNum+1)
			 trendP=matrix(0,nrow=1,ncol=repeatNum+1)
		     temp=Armitage(genotypeStat)
		     trendValue[1,1]=temp$statistic
			 trendP[1,1]=temp$pValue
		 }
     )
	 
	 ## null distribution calculate
	 for(i in 2:(repeatNum+1)){
	     randMatrix=permuteMatrix(genotypeStat)
		 randCount=randMatrix
	     delIndex=c()
	     for(m in 1:3){
	         if(all(randCount[,m]==0)){
		         delIndex=c(delIndex,m)
		     }
	     }
		 if(!is.null(delIndex)){
		     randCount=randCount[,-delIndex,drop=FALSE]
		 }
	     
	     switch(testMethod,
             chisq={
	             temp=chisq.test(randCount)
			     chiValue[1,i]=temp$statistic
			     chiP[1,i]=temp$p.value
	         },
             fisher={
	             fisherP[1,i]=fisher.test(randCount)$p.value
	         },
	         OR={
		         caseAllele1=2*randMatrix[1,1]+randMatrix[1,2]
			     caseAllele2=2*randMatrix[1,3]+randMatrix[1,2]
			     controlAllele1=2*randMatrix[2,1]+randMatrix[2,2]
			     controlAllele2=2*randMatrix[2,3]+randMatrix[2,2]		   
			     temp=OR(caseAllele1,caseAllele2,controlAllele1,controlAllele2)
			     ORvalue[1,i]=temp$OR
			    # CIvalue[1,i]=temp$CI[1,2]-temp$CI[1,1]
		     },
		     Armitage={
		         temp=Armitage(randMatrix)
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
			 list(pValue=pValue,obsP=fisherP[1,1],permP=fisherP[,2:(repeatNum+1)])
	     },
	     OR={
	         obsOR=ORvalue[1,1]
			 maxExp=as.matrix(ORvalue[ORvalue>obsOR])
			 ORp=nrow(maxExp)/repeatNum
			# obsCI=CIvalue[1,1]
			# minExp=as.matrix(CIvalue[CIvalue<obsCI])
			# CIp=nrow(minExp)/repeatNum
			# list(ORp=ORp,CIp=CIp,obsOR=obsOR,obsCI=obsCI,permOR=ORvalue[,2:(repeatNum+1)],permCI=CIvalue[,2:(repeatNum+1)])
			list(pValue=ORp,obsOR=obsOR,permOR=ORvalue[,2:(repeatNum+1)])
	     },
		 Armitage={
		     obs=trendValue[1,1]
			 maxExp=as.matrix(trendValue[trendValue>obs])
		     pValue=nrow(maxExp)/repeatNum
			 list(pValue=pValue,obsStatistic=obs,obsP=trendP[1,1],permStatistic=trendValue[,2:(repeatNum+1)],permP=trendP[,2:(repeatNum+1)])

		 }
     )
}
