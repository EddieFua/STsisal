getCellNumber_st = function (Y.raw, possibleCellNumber = 2:15) 
{
  allAIC = c()
  for (K in possibleCellNumber) {
    out <- csDeconv_st(Y.raw, K = K, TotalIter = 5, bound_negative = TRUE)
    estProp = mysisal(Y.raw[Out$InitMarker, ][Out$updatedInx, ], K = K, topN = 50)$estProp
    aic = compute_aic(estProp, Y.raw[Out$InitMarker, ][Out$updatedInx, ])
    allAIC = append(allAIC, aic)
  }
  bestK = possibleCellNumber[which.min(allAIC)]
  return(list(bestK = bestK, allAIC = allAIC))
}