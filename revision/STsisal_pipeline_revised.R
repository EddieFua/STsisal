source('/Users/fuyinghao/Documents/STsisal/function/clean_count.R')
source('/Users/fuyinghao/Documents/STsisal/function/getCellNumber_st.R')
source('/Users/fuyinghao/Documents/STsisal/function/csDeconv_st.R')
source('/Users/fuyinghao/Documents/STsisal/function/mysisal.R')
source('/Users/fuyinghao/Documents/STsisal/function/simplenormalize.R')
source('/Users/fuyinghao/Documents/STsisal/function/sisal.R')
source('/Users/fuyinghao/Documents/STsisal/function/vca.R')
source('/Users/fuyinghao/Documents/STsisal/function/findRefinx.R')
source('/Users/fuyinghao/Documents/STsisal/Revision/function/deconfounding_revised.R')
library(deconf)
library(TOAST)
CornerToEstProp = function (corner)
{
  N_sample = dim(corner)[1]
  tmp <- nnls(corner, rep(1, N_sample))
  estProp <-
    diag(as.numeric(tmp$x), length(as.numeric(tmp$x))) %*% t(corner)
  estProp[estProp < 0] = 0
  estProp[estProp > 1] = 1
  for (i in 1:ncol(estProp)) {
    estProp[, i] = estProp[, i] / sum(estProp[, i])
  }
  return(t(estProp))
}
CornerToMarker = function (distances, topN)
{
  markerList <- apply(distances, 2, function(xx) {
    pure <- rownames(distances)[order(xx)[1:topN]]
    return(pure)
  })
  markerList <-
    split(markerList, rep(1:ncol(markerList), each = nrow(markerList)))
  return(markerList)
}

Rsumlog = function (a) 
{
  s <- a[1]
  for (i in 2:length(a)) s <- Raddlog(s, a[i])
  s
}

Raddlog = function (a, b) 
{
  result <- rep(0, length(a))
  idx1 <- a > b + 200
  result[idx1] <- a[idx1]
  idx2 <- b > a + 200
  result[idx2] <- b[idx2]
  idx0 <- !(idx1 | idx2)
  result[idx0] <- a[idx0] + log1p(exp(b[idx0] - a[idx0]))
  result
}

MyNaiveBayes = function (selProf, knowRef) 
{
  nbres <- matrix(0, ncol(knowRef), ncol(selProf))
  for (ii in 1:ncol(selProf)) {
    initrec <- rep(0, ncol(knowRef))
    for (jj in 1:ncol(knowRef)) {
      initrec[jj] <- -sum((knowRef[, jj] - selProf[, ii])^2)
    }
    nbres[, ii] <- exp(initrec)/exp(Rsumlog(initrec))
  }
  return(nbres)
}

GetCorRes = function (selProf, knowRefAll) 
{
  isc = intersect(rownames(selProf), rownames(knowRefAll))
  selProf = selProf[isc, ]
  knowRef <- knowRefAll[isc, ]
  cormat <- cor(knowRef, selProf)
  nbmat <- MyNaiveBayes(selProf, knowRef)
  summat <- (cormat + nbmat)/2
  summat_org <- summat
  assignLabel <- rep("Unassigned", ncol(selProf))
  for (i in 1:ncol(knowRefAll)) {
    thisidx <- which(summat == max(summat), arr.ind = TRUE)
    assignLabel[thisidx[2]] <- rownames(summat)[thisidx[1]]
    summat[thisidx[1], ] <- -1
    summat[, thisidx[2]] <- -1
  }
  return(list(assignLabel = assignLabel, probMat = summat_org))
}
source('/Users/fuyinghao/Documents/STsisal/Revision/function/csDeconv_st_revised.R')
STsisal_pipeline_revised = function(Y_raw,S.start, C.start,K = NULL, knowRef = NULL, possibleCellNumber = 3:15,n_marker = 1000,TotalIter = 5){
  Y_raw = clean_count(Y_raw)
  if (is.null(K)) {
    Kres = getCellNumber_st(Y_raw, possibleCellNumber)
    K = Kres$bestK
  }
  OutRF <- csDeconv_st_revised(Y_raw,S.start, C.start, K = K, TotalIter = TotalIter, bound_negative = TRUE,nMarker = n_marker)
  reduceYmat <- Y_raw[OutRF$InitMarker, ]
  for (i in 1:length(OutRF$updatedInx)) {
    reduceYmat <- reduceYmat[OutRF$updatedInx[[i]], ]
  }
  tsisal_res <- mysisal(reduceYmat, K = K, topN = 50)
  estProp <- tsisal_res$estProp
  selMarker <- tsisal_res$selMarker
  if (!is.null(knowRef)) {
    out_all <- mycsfit(estProp, t(Y_raw))
    prof <- t(out_all$ghat)
    rownames(prof) <- rownames(Y_raw)
    selProf <- prof[unlist(selMarker), ]
    labres <- GetCorRes(selProf, knowRefAll = knowRef)
    colnames(estProp) <- labres$assignLabel
  }
  return(list(estProp = estProp, selMarker = selMarker, K = K))
}







