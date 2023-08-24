### get proportion aligned using correlation coefficient
GetPropAligned <- function(input,reference,L=K){
  colnames(input)=colnames(reference)=seq(1,dim(input)[2],by=1)
  input_sum = colSums(input)
  idx = which(input_sum==0)
  if (length(idx)>0){
    for (i in 1:(length(idx))) {
      input[i, idx[[i]]] = 0.0001
    }
  }
  corMat = cor(input,reference,use="pairwise.complete.obs")
  prop_cor = rep(0,L)
  tmpmat=corMat
  tmpmat[is.na(tmpmat)] = rnorm(L,-1,0.01)
  if(L>2){
    for(i in 1:L){
      maxind = which(tmpmat == max(tmpmat), arr.ind = TRUE)[1,]
      prop_cor[maxind[1]] = colnames(corMat)[maxind[2]]
      tmpmat[maxind[1],]=tmpmat[,maxind[2]]=rep(-10,L)
      #print(prop_cor)
    }
  }else if(L==2){
    if(tmpmat[1,1]>0){
      prop_cor = c("1","2")
    }else{
      prop_cor = c("2","1")
    }
  }
  colnames(input) = prop_cor
  trans_input = input[,colnames(reference)]
  return(trans_input)
}