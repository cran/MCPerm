OR <-
function(caseAllele1,caseAllele2,controlAllele1,controlAllele2){
     allele1=caseAllele1+controlAllele1
	 allele2=caseAllele2+controlAllele2
	 if(allele1>allele2){
	     OR=(caseAllele2*controlAllele1)/(caseAllele1*controlAllele2)
	 }else{
	     OR=(caseAllele1*controlAllele2)/(caseAllele2*controlAllele1)
	 }
     LnOR=log(OR)
	 SELnOR=sqrt(1/caseAllele1+1/caseAllele2+1/controlAllele1+1/controlAllele2)
	 upLnOR=LnOR+1.96*SELnOR
	 lowLnOR=LnOR-1.96*SELnOR
	 lowerOR=exp(lowLnOR)
	 upperOR=exp(upLnOR)
	 CI=matrix(c(lowerOR,upperOR),nrow=1)
	 list(OR=OR,CI=CI)
}
