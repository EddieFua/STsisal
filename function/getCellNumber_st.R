source('/Users/fuyinghao/Documents/STsisal/function/compute_aic.R')
source('/Users/fuyinghao/Documents/STsisal/function/mycsfit.R')
source('/Users/fuyinghao/Documents/STsisal/function/clean_count.R')
getCellNumber_st = function (Y.raw, possibleCellNumber = 2:15) 
{
  allAIC = c()
  for (K in possibleCellNumber) {
    Y.raw <- clean_count(Y.raw)
    out <- csDeconv_st(Y.raw, K = K, TotalIter = 5, bound_negative = TRUE)
    estProp = mysisal(Y.raw[out$InitMarker, ][out$updatedInx, ], K = K, topN = 50)$estProp
    aic = compute_aic(estProp, Y.raw[out$InitMarker, ][out$updatedInx, ])
    allAIC = append(allAIC, aic)
  }
  bestK = possibleCellNumber[which.min(allAIC)]
  return(list(bestK = bestK, allAIC = allAIC))
}