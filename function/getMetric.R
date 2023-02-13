library(philentropy)
getMetric <- function(trueProp,estProp){
  
  ##################################################
  ## measure each spot’s deconvolution accuracy
  ##################################################
  Nct = ncol(trueProp)
  Nst = nrow(trueProp)
  
  ## JSD
  jsd_vec = unlist(lapply(1:Nst,function(ist){
    jsd = suppressMessages(JSD(x = rbind(trueProp[ist,],estProp[ist,]), unit = "log2",est.prob = "empirical"))
  }))
  
  ## rmse
  rmse_vec = unlist(lapply(1:Nst,function(ist){
    rmse = sqrt(sum((trueProp[ist,]-estProp[ist,])^2)/Nct)
  }))
  
  ## RMSD
  mAD_vec = unlist(lapply(1:Nst,function(ist){
    mAD = mean(abs(trueProp[ist,]-estProp[ist,]))
  }))
  
  ## correlation
  cor_vec = unlist(lapply(1:Nst,function(ist){
    corr = cor(trueProp[ist,],estProp[ist,])
  }))
  
  ######################################################################
  ## measure the performance in distinguishing different cell types
  ######################################################################
  corr_ct = diag(cor(trueProp,estProp))
  MAE_ct = colSums(abs(trueProp - estProp)) / nrow(trueProp)
  
  return(list(jsd_vec = jsd_vec,rmse_vec = rmse_vec,mAD_vec = mAD_vec,cor_vec = cor_vec,corr_ct = corr_ct,MAE_ct = MAE_ct))
}