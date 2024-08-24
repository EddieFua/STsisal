rm(list = ls())
library(ape)
source('/Users/fuyinghao/Documents/STsisal/function/sim_functions.R')
## mOB single cell and spatial data for simulation
load('/Users/fuyinghao/Documents/STsisal/sim_data/sim_raw_data/pattern_gp_label.RData')
##position
load('/Users/fuyinghao/Documents/STsisal/sim_data/sim_raw_data/sim_MOB_location.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/sim_raw_data/split.scRNAseq.forSimu.RData')
readRDS('/Users/fuyinghao/Documents/STsisal/sim_data/sim_raw_data/GSE109447.sceset.scenario5.RDS')
ct_select = levels(split$eset.sub.split1@phenoData@data$cellType)[-c(2,4)]
c_noise = c(0,0.2,0.4)
library(Matrix)
library(TOAST)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CARD)
library(spacexr)
library(SPOTlight)
library(Giotto)
library(MuSiC)
library(STdeconvolve)

source('/Users/fuyinghao/Documents/STsisal/pipeline/STsisal_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/pipeline/Deconvolution_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/function/getMetric.R')
source('/Users/fuyinghao/Documents/STsisal/function/GetPropAligned.R')
i_STsisal_rmse = list()
i_CARD_rmse = list()
i_RCTD_rmse = list()
i_STdeconvolve_rmse = list()

i_STsisal_corr = list()
i_CARD_corr = list()
i_RCTD_corr = list()
i_STdeconvolve_corr = list()

i_STsisal_cor_vec=list()
i_CARD_cor_vec = list()
i_RCTD_cor_vec = list()
i_STdeconvolve_cor_vec = list()

i_STsisal_jsd = list()
i_CARD_jsd = list()
i_RCTD_jsd = list()
i_STdeconvolve_jsd = list()

i_STsisal_MAE = list()
i_CARD_MAE = list()
i_RCTD_MAE = list()
i_STdeconvolve_MAE = list()
for (i in 1:3) {
  for (j in 1:20) {
    simdata = getSpotMix(split$eset.sub.split1,pattern_gp_label,ct_select,noise_p = c_noise[i])
    obs.Y = simdata$obs.Y
    true.prop = simdata$true.prop
    write.table(obs.Y,file = paste0('/Users/fuyinghao/Documents/STsisal/sim_data/sim_data_2/',c_noise[i],'/',j,'.txt'))
    write.table(true.prop,file = paste0('/Users/fuyinghao/Documents/STsisal/sim_data/sim_data_2/',c_noise[i],'/Y_',j,'.txt'))
  }
}

load('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res.RData')

for(i in 1:3){
  j_STsisal = list()
  j_CARD = list()
  j_RCTD = list()
  j_STdeconvolve = list()
  for (j in 1:20) {
    sc_meta = split$eset.sub.split1@phenoData@data
    sc_count = split$eset.sub.split1@assayData$exprs
    colnames(sc_meta) = c("cellname","cellType","sampleInfo")
    idx = which(sc_meta$cellType%in%ct_select)
    sc_meta = sc_meta[idx,]
    sc_count = sc_count[,idx]
    sc_meta$cellType = as.character(sc_meta$cellType)
    sc_meta$cellType = as.factor(sc_meta$cellType)
    spatial_location = location
    ##generate sim data
    obs.Y = read.table(paste0('/Users/fuyinghao/Documents/STsisal/sim_data/sim_data_2/',c_noise[i],'/',j,'.txt'))
    true.prop = read.table(paste0('/Users/fuyinghao/Documents/STsisal/sim_data/sim_data_2/',c_noise[i],'/Y_',j,'.txt'))
    simdata = list(obs.Y = as.matrix(obs.Y),true.prop = as.matrix(true.prop))
    rownames(spatial_location) = colnames(simdata$obs.Y)
    simdata$true.prop[1, which(colSums(simdata$true.prop) == 0)] = 0.00001
    simdata$true.prop = simdata$true.prop/rowSums(simdata$true.prop)
    ##STsisal deconvolve
    STsisal_res = STsisal_pipeline(simdata$obs.Y,6)
    STsisal_res$estProp = GetPropAligned(STsisal_res$estProp,simdata$true.prop,6)
    j_STsisal[[j]] = getMetric(simdata$true.prop,STsisal_res$estProp)  
    ##CARD deconvolve
    CARD_res = CARD_pipeline(sc_count,sc_meta,simdata$obs.Y,spatial_location)
    CARD_res = CARD_res[,colnames(simdata$true.prop)]
    j_CARD[[j]] = getMetric(simdata$true.prop,CARD_res)
    ##RCTD deconvolve
    RCTD_res = RCTD_pipeline(sc_count,sc_meta,simdata$obs.Y,spatial_location)
    j_RCTD[[j]] = getMetric(simdata$true.prop,RCTD_res)
    # ##Spotlight deconvolve
    # SPOTlight_res = SPOTlight_pipeline(sc_count,sc_meta,simdata$obs.Y,spatial_location)
    # j_SPOTlight[[j]] = getMetric(simdata$true.prop,SPOTlight_res)
    ##STdeconvolve deconvolve
    STdeconvolve_res = STdeconvolve_pipeline(simdata$obs.Y,6)
    STdeconvolve_res = GetPropAligned(STdeconvolve_res,simdata$true.prop,6)
    j_STdeconvolve[[j]] = getMetric(simdata$true.prop,STdeconvolve_res)
  }
  #rmse
  i_STsisal_rmse[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$rmse_vec))
  i_CARD_rmse[[i]] = sapply(1:20, function(i)mean(j_CARD[[i]]$rmse_vec))
  i_RCTD_rmse[[i]] = sapply(1:20, function(i)mean(j_RCTD[[i]]$rmse_vec))
  i_STdeconvolve_rmse[[i]] = sapply(1:20, function(i)mean(j_STdeconvolve[[i]]$rmse_vec))
  #cor_cell type
  i_STsisal_corr[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$corr_ct))
  i_CARD_corr[[i]] = sapply(1:20, function(i)mean(j_CARD[[i]]$corr_ct))
  i_RCTD_corr[[i]] = sapply(1:20, function(i)mean(j_RCTD[[i]]$corr_ct))
  i_STdeconvolve_corr[[i]] = sapply(1:20, function(i)mean(j_STdeconvolve[[i]]$corr_ct))
  #cor_sample
  i_STsisal_cor_vec[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$cor_vec))
  i_CARD_cor_vec[[i]] = sapply(1:20, function(i)mean(j_CARD[[i]]$cor_vec))
  i_RCTD_cor_vec[[i]] = sapply(1:20, function(i)mean(j_RCTD[[i]]$cor_vec))
  i_STdeconvolve_cor_vec[[i]] = sapply(1:20, function(i)mean(j_STdeconvolve[[i]]$cor_vec))
  #jsd
  i_STsisal_jsd[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$jsd_vec))
  i_CARD_jsd[[i]] = sapply(1:20, function(i)mean(j_CARD[[i]]$jsd_vec))
  i_RCTD_jsd[[i]] = sapply(1:20, function(i)mean(j_RCTD[[i]]$jsd_vec))
  i_STdeconvolve_jsd[[i]] = sapply(1:20, function(i)mean(j_STdeconvolve[[i]]$jsd_vec))
  #MAE
  i_STsisal_MAE[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$MAE_ct))
  i_CARD_MAE[[i]] = sapply(1:20, function(i)mean(j_CARD[[i]]$MAE_ct))
  i_RCTD_MAE[[i]] = sapply(1:20, function(i)mean(j_RCTD[[i]]$MAE_ct))
  i_STdeconvolve_MAE[[i]] = sapply(1:20, function(i)mean(j_STdeconvolve[[i]]$MAE_ct))
}

save(i_STsisal_rmse,
     i_CARD_rmse,
     i_RCTD_rmse,
     i_STdeconvolve_rmse,
     i_STsisal_corr,
     i_CARD_corr,
     i_RCTD_corr,
     i_STdeconvolve_corr,
     i_STsisal_cor_vec,
     i_CARD_cor_vec,
     i_RCTD_cor_vec,
     i_STdeconvolve_cor_vec,
     i_STsisal_jsd ,
     i_CARD_jsd,
     i_RCTD_jsd,
     i_STdeconvolve_jsd,
     i_STsisal_MAE,
     i_CARD_MAE,
     i_RCTD_MAE,
     i_STdeconvolve_MAE,
     file = '/Users/fuyinghao/Documents/STsisal/exp_res/sim_res.RData')

load('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res.RData')

rownames(spatial_location) = rownames(simdata$true.prop)
rownames(STsisal_res$estProp)  = rownames(spatial_location)
rownames(true.prop)  = rownames(spatial_location)
spatial_pieplot(true.prop,spatial_location,colors)
source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
colors = RColorBrewer::brewer.pal(12,'Set3')[1:6]
i_STsisal_corr[[1]]
spatial_pieplot(STsisal_res$estProp,spatial_location,colors)
spatial_pieplot(simdata$true.prop,spatial_location,colors)
mean(getMetric(simdata$true.prop,STsisal_res$estProp)$rmse_vec)
spatial_pieplot(STdeconvolve_res ,location,colors)
spatial_pieplot(RCTD_res,location,colors)
spatial_pieplot(CARD_res,location,colors)
spatial_pieplot(RCTD_res,location,colors)

max_type = function(prop){
  idx = apply(prop,1,function(i) which(i == max(i))[1],simplify = TRUE)
  for (i in 1:length(idx)) {
    tmp = rep(0,ncol(prop))
    tmp[idx[i]] = 1
    prop[i,] = tmp
  }
  return(prop)
}
b = max_type(STsisal_res$estProp)
a = max_type(simdata$true.prop)
spatial_pieplot(a,spatial_location,colors)
spatial_pieplot(b,spatial_location,colors)
