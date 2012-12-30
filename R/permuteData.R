permuteData <-
function(dataLine,fromCol){
     if(is.null(dataLine)){
	     stop("'dataLine' cannot be empty.")
	 }
	 if(!is.data.frame(dataLine) && !is.matrix(dataLine)){
	     stop("'dataLine' must be a data.frame or matrix.")
	 }
	 if(nrow(dataLine)!=1){
	     stop("'dataLine' must be a data.frame or matrix with nrow=1.")
	 }
	 
	 if(is.null(fromCol)){
	     stop("'fromCol' cannot be empty.")
	 }
	 if(abs(fromCol-round(fromCol)) > 1e-10){
	     stop("'fromCol' must be a integer.")
	 }
	 if(is.matrix(dataLine)){
	     dataLine=as.data.frame(dataLine)
	 }

	 colNum=ncol(dataLine)
	 randIndex=sample(colNum-fromCol+1)
	 tempData=dataLine[,fromCol:colNum]
	 tempData=tempData[,randIndex]
	 newDataLine=data.frame(dataLine[,1:(fromCol-1)],tempData)
	 # colnames(dataLine)=paste("x",1:colNum,sep="")
	 colnames(newDataLine)=paste("x",1:colNum,sep="")
	 # list(rawData=dataLine,newData=newDataLine)
         return(newDataLine)
}
