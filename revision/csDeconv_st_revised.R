normalize_col = function(Y){
  Y.norm = Y/Matrix::colSums(Y)
  return(Y.norm)
}
rowsum_0 = function(x){
  a = sum(x!=0)
  return(a)
}
source('/Users/fuyinghao/Documents/STsisal/Revision/function/deconfounding_revised.R')
csDeconv_st_revised = function (Y_raw,S.start,
                                K,
                                nMarker = 1000,
                                InitMarker = NULL,
                                TotalIter = 30,
                                bound_negative = FALSE)
{
  # if (is(Y_raw, "SummarizedExperiment")) {
  #   se <- Y_raw
  #   Y_raw <- assays(se)$counts
  # } else if (!is(Y_raw, "matrix")) {
  #   stop("Y_raw should be a matrix or a SummarizedExperiment object!")
  # }
  # if (is.null(rownames(Y_raw))) {
  #   row.names(Y_raw) <- seq(nrow(Y_raw))
  # }
  # if (is.null(InitMarker)) {
  #   if (nrow(Y_raw) < 2000) {
  #     InitMarker <- findRefinx(Y_raw, nmarker = nMarker)
  #   }
  #   else {
  #     tmp <- findRefinx(Y_raw, nmarker = nMarker * 2)
  #     # InitMarker <- tmp[nMarker + 1:nMarker]
  #     InitMarker <- tmp[1:nMarker]
  #   }
  # } else {
  #   if (sum(!(InitMarker %in% rownames(Y_raw))) > 0) {
  #     stop("Discrepancy between\n                    InitMarker and the row names of Y_raw!")
  #   }
  # }
  allProp <- list()
  updatedInx_list <- list()
  allRMSE <- rep(0, TotalIter + 1)
  Y <- Y_raw
  S.start = S.start[match(rownames(Y),rownames(S.start)),]
  Prop0 <- deconfounding(Y, K)
  allProp[[1]] <- Prop0$C$Matrix
  Y = clean_count(Y)
  out_all <- mycsfit(t(allProp[[1]]), t(Y))
  prof <- t(out_all$ghat)
  tmpmat <- prof %*% allProp[[1]]
  allRMSE[1] <- sqrt(mean((t(Y) - t(tmpmat)) ^ 2))
  message("+========================================+")
  message("+======= Total iterations = ", TotalIter, " ==========+")
  for (i in seq_len(TotalIter)) {
    message("Current iter = ", i)
    # if(qr(t(allProp[[i]])%*%allProp[[i]])$rank!=K){
    #   change_index <- which(apply(t(allProp[[i]]), 2,sum)==min(apply(t(allProp[[i]]),2,sum)))
    #   allProp[[i]][change_index,] <- runif(ncol(allProp[[i]])*length(change_index),0.01,0.02)
    #   allProp[[i]] <- normalize_col(allProp[[i]])
    # }
    updatedInx_list[[i]] <- DEVarSelect(Y, t(allProp[[i]]), nMarker,bound_negative = TRUE)
    Y <- Y[updatedInx_list[[i]],]
    Prop0 <- deconfounding(Y, K)
    allProp[[i + 1]] <- Prop0$C$Matrix
    outall = mycsfit(t(allProp[[i + 1]]), t(Y))
    prof <- t(outall$ghat)
    tmpmat <- prof %*% allProp[[i+1]]
    allRMSE[i + 1] <- sqrt(mean((t(Y) - t(tmpmat)) ^ 2))
  }
  min_idx <- which.min(allRMSE)
  Prop0 <- allProp[[min_idx]]
  if(min_idx == 1){
    updatedInx = 1:length(InitMarker)
  }else{ updatedInx = updatedInx_list[[min_idx-1]]
  }
  return(list(
    InitMarker = InitMarker,
    allRMSE = allRMSE,
    allProp = allProp,
    estProp = Prop0,
    updatedInx = updatedInx
  ))
}
