permuteMatrix <-
function(rawMatrix){
     if(!is.matrix(rawMatrix)){
	     stop("'rawMatrix' must be a matrix.")
	 }
     if(nrow(rawMatrix)!=2){
	     stop("'rawMatrix' must be a matrix with nrow=2.")
	 }
	 
	 colNum=ncol(rawMatrix)
	 matrixColSum=matrix(0,nrow=1,ncol=colNum)
	 for(i in 1:colNum){
	     matrixColSum[1,i]=rawMatrix[1,i]+rawMatrix[2,i]
	 }
	 matrixRowSum=rowSums(rawMatrix)
	 
	 newMatrix=matrix(0,nrow=2,ncol=colNum)
	 tempM=matrixRowSum[1][1]
	 tempN=matrixRowSum[2][1]
	 for(j in 1:colNum){
	     newMatrix[1,j]=rhyper(1,m=tempM,n=tempN,k=matrixColSum[1,j])
	     newMatrix[2,j]=matrixColSum[1,j]-newMatrix[1,j]
		 tempM=tempM-newMatrix[1,j]
		 tempN=tempN-newMatrix[2,j]
	 }
     return(newMatrix)
}
