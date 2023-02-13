library(matrixStats)
library(MCMCpack)

###############################
## get cell type proportion of each layer
###############################
getOneLayerProp <- function(pattern_gp_label,ipt,ct.select,noise_p,ntotal = 10){
  
  nobs = sum(pattern_gp_label == ipt) ## total number of spot in this layer
  numCT = sample(1:length(ct.select),1) ## total number of cell types in each layer: main cell type + colocalized cell types
  
  if(ipt == 1){
    main_ct = "Neurons" ### defined one dominant cell type in layer 1
    concen = rep(1,numCT) ### set the concentration parameters, currently fixed to be 1 for all cell types
  }else if (ipt == 2){ 
    main_ct = "Astrocytes" ### defined one dominant cell type in layer 2
    concen = rep(1,numCT) 
  }else if(ipt == 3){
    main_ct = "Oligos" ### defined one dominant cell type in layer 3
    concen = rep(1,numCT) 
  }
  
  ## proportion of all spots in this layer
  prop_this_layer = rdirichlet(nobs,concen)
  rownames(prop_this_layer) = names(pattern_gp_label)[pattern_gp_label == ipt]
  ct_this_layer = c(main_ct,sample(ct.select[ct.select != main_ct],numCT-1))
  
  mix_spot = sample(rownames(prop_this_layer),round(nobs * noise_p)) ## noisy spots
  fix_spot = setdiff(rownames(prop_this_layer),mix_spot)
  if(length(mix_spot) > 0){
    prop_this_layer[mix_spot,] = rdirichlet(length(mix_spot),rep(1,numCT))
  }
  
  if(length(fix_spot) > 0){
    prop_this_layer[fix_spot,] = t(sapply(fix_spot,function(i){
      prop_this_layer[i,][order(prop_this_layer[i,],decreasing = T)] ### for non-noisy positions, order the proportion to assign the largest proportion to the dominant cell type
    }))
  }
  colnames(prop_this_layer) = ct_this_layer 
  
  ## We fixed the total number of cells on each spot to be 10
  sample = round(ntotal * prop_this_layer,digits = 0)
  return(sample)
}

## get cell type proportion matrix of all layers
getProportion <- function(pattern_gp_label,ct.select,noise_p = 0.2,ntotal = 10){
  
  ## proportion matrix
  sample.matrix = matrix(0,nrow = length(pattern_gp_label),ncol = length(ct.select))
  rownames(sample.matrix) = names(pattern_gp_label)
  colnames(sample.matrix) = ct.select

  pattern1 = getOneLayerProp(pattern_gp_label,1,ct.select,noise_p = noise_p,ntotal=10)
  pattern2 = getOneLayerProp(pattern_gp_label,2,ct.select,noise_p = noise_p,ntotal=10)
  pattern3 = getOneLayerProp(pattern_gp_label,3,ct.select,noise_p = noise_p,ntotal=10)
  
  idx1 = match(colnames(sample.matrix),colnames(pattern1))
  idx2 = match(colnames(sample.matrix),colnames(pattern2))
  idx3 = match(colnames(sample.matrix),colnames(pattern3))
  sample.matrix[pattern_gp_label == 1,] = pattern1[,idx1]
  sample.matrix[pattern_gp_label == 2,] = pattern2[,idx2]
  sample.matrix[pattern_gp_label == 3,] = pattern3[,idx3]
  
  sample.matrix[is.na(sample.matrix)] = 0
  prop.matrix = sweep(sample.matrix,1,rowSums(sample.matrix),"/")
  return(res = list(sample.matrix = sample.matrix,prop.matrix = prop.matrix))
}

###############################
## get gene expression profiles of all spots
##############################
getSpotMix <- function(eset.sub.split1,pattern_gp_label,ct.select,noise_p = 0.2,ntotal = 10){
  
  phenoData <- eset.sub.split1@phenoData@data ## using split 1 to sample the single cell RNAseq data
  K <- length(unique(ct.select)) ## number of cell types
  message(paste('Using',K,'cell types to generate pseudo spatial dataset'))
  
  ## get proportion matrix
  out = getProportion(pattern_gp_label,ct.select,noise_p = noise_p,ntotal = 10)
  sample.matrix = out$sample.matrix
  true.prop = out$prop.matrix
  
  ## generate ST data
  nGene = nrow(exprs(eset.sub.split1))
  nSpot = nrow(true.prop)
  obs.Y = matrix(0,nrow = nGene,ncol = nSpot)
  rownames(obs.Y) = rownames(exprs(eset.sub.split1))
  colnames(obs.Y) = rownames(true.prop)
  
  for(i in 1:nrow(sample.matrix)){
    one_spot_expr = c()
    for(ct in ct.select){
      nCell = sample.matrix[i,ct]
      if(nCell > 0 & nCell <= sum(phenoData$cellType == ct)){
        index = sample(which(phenoData$cellType == ct),nCell,replace = FALSE)
      }else if(nCell > 0 & nCell > sum(phenoData$cellType == ct)){
        index = sample(which(phenoData$cellType == ct),nCell,replace = TRUE)
      }else{
        index = 0
      }
      one_spot_expr = cbind(one_spot_expr,exprs(eset.sub.split1)[,index])
    }
    obs.Y[,i] = rowSums(one_spot_expr)
  }
  return(res = list(obs.Y = obs.Y,true.prop = true.prop,sample.matrix = sample.matrix))
}



