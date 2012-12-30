genotypeStat <-
function(genotypeLine,affectionLine,fromCol,naString,sep){
     if(is.null(genotypeLine)){
	     stop("'genotypeLine' cannot be empty.")
	 }
	 if(!is.data.frame(genotypeLine) && !is.matrix(genotypeLine)){
	     stop("'genotypeLine' must be a data.frame or matrix.")
	 }
	 if(nrow(genotypeLine)!=1){
	     stop("'genotypeLine' must be a data.frame or matrix with nrow=1.")
	 }
	 
	 if(is.null(affectionLine)){
	     stop("'affectionLine' cannot be empty.")
	 }
	 if(!is.data.frame(affectionLine) && !is.matrix(affectionLine)){
	     stop("'affectionLine' must be a data.frame or matrix.")
	 }
	 if(nrow(affectionLine)!=1){
	     stop("'affectionLine' must be a data.frame or matrix with nrow=1.")
	 }
	 
	 if(ncol(genotypeLine)!=ncol(affectionLine)){
	     stop("'genotypeLine' and 'affectionLine' must have the same length.")
	 }
	 
	 if(is.null(fromCol)){
	     stop("'fromCol' cannot be empty.")
	 }
	 if(abs(fromCol-round(fromCol)) > 1e-10){
	     stop("'fromCol' is the index of the first genotype data in the 'genotypeLine',so must be a integer.")
	 }
	 
	 if(is.null(naString)){
	     stop("'naString' cannot be empty.")
	 }
	 if(!is.character(naString)){
	     stop("'naString' must be a string.")
	 }
	 
	 if(is.null(sep)){
	     stop("'sep' cannot be empty.")
	 }
	 if(!is.character(sep)){
	     stop("'sep' must be a string.")
	 }
	 
	 if(is.matrix(genotypeLine)){
	     genotypeLine=as.data.frame(genotypeLine)
	 }
     if(is.matrix(affectionLine)){
	     affectionLine=as.data.frame(affectionLine)
	 }
     
	 # stat
	 colNum=ncol(genotypeLine)
	 genotypeCount=as.matrix(genotypeLine[,fromCol:colNum])
	 affection=as.matrix(affectionLine[,fromCol:colNum])
	 tempResult=as.matrix(table(affection,genotypeCount))
	 colNum=ncol(tempResult)
	 colName=colnames(tempResult)
	 
	 # deal with NA
	 existNA=0
	 if(any(colName==naString)){
	     existNA=1
		 if(colNum==1){
		     return(tempResult)
		 }
		 naIndex=which(colName==naString)
		 naResult=tempResult[,naIndex,drop=FALSE]
		 tempResult=tempResult[,-naIndex]
		 colName=colName[-naIndex]
		 colNum=ncol(tempResult)
	 }
	 
	 # deal with AG==GA
	 delIndex=c()
	 for(i in 1:colNum){
	     colNameSplit=unlist(strsplit(colName[i],sep))
		 if(length(colNameSplit)!=2){
		     stop("'sep' or 'naString' or 'genotypeLine' is wrong.Please see the putting data,then correct the parameter.")
		 }
	     if(colNameSplit[1]!=colNameSplit[2]){
		     colNameSplit=sort(colNameSplit)
			 colName[i]=paste(colNameSplit[1],colNameSplit[2],sep=sep)
			 colNameIndex=which(colName==colName[i])
			 if(length(colNameIndex)==2){
			     tempResult[,colNameIndex[1]]=tempResult[,colNameIndex[1]]+tempResult[,colNameIndex[2]]
				 delIndex=c(delIndex,colNameIndex[2])
			 }
		 }
	 }
	 if(!is.null(delIndex)){
	     tempResult=tempResult[,-delIndex]
	     colName=colName[-delIndex]
	 }
	 colnames(tempResult)=colName
     colName=sort(colName)
	 tempResult=tempResult[,colName]
	 colNum=ncol(tempResult)
	 
	 ### deal with allele
	 alleleName=c()
	 rowNum=nrow(tempResult)
	 alleleCount=matrix(nrow=rowNum,ncol=colNum*2)
	 for(i in 1:colNum){
	     alleleName=c(alleleName,unlist(strsplit(colName[i],sep)))
	 }
	 j=1
	 for(i in 1:colNum){
	     alleleCount[,j]=tempResult[,i]
		 j=j+1
		 alleleCount[,j]=tempResult[,i]
		 j=j+1
	 }
	 colnames(alleleCount)=alleleName
	 
	 allele=unique(alleleName)
	 alleleNum=length(allele)
	 alleleResult=matrix(nrow=rowNum,ncol=alleleNum)
	 for(i in 1:alleleNum){
	     index=which(colnames(alleleCount)==allele[i])
		 alleleResult[,i]=as.matrix(rowSums(alleleCount[,index,drop=FALSE]))
	 }
	 colnames(alleleResult)=allele
	 rownames(alleleResult)=rownames(tempResult)
	 
	 # return value
	 if(existNA==1){
	     genotypeResult=cbind(tempResult,naResult)
		 colnames(genotypeResult)=c(colName,naString)
	 }else{
	     genotypeResult=tempResult
	 }
	 
	 list(alleleStat=alleleResult,genotypeStat=genotypeResult)
}
