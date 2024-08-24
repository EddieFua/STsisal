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
# c_noise = c(0,0.2,0.4)
c_noise = c(0,0.3,0.6)
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

source('/Users/fuyinghao/Documents/STsisal/Revision/function/STsisal_pipeline_revised.R')
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
# for (i in 1:3) {
#   for (j in 1:20) {
#     simdata = getSpotMix(split$eset.sub.split1,pattern_gp_label,ct_select,noise_p = c_noise[i])
#     obs.Y = simdata$obs.Y
#     true.prop = simdata$true.prop
#     write.table(obs.Y,file = paste0('/Users/fuyinghao/Documents/STsisal/sim_data/sim_data_2/',c_noise[i],'/',j,'.txt'))
#     write.table(true.prop,file = paste0('/Users/fuyinghao/Documents/STsisal/sim_data/sim_data_2/',c_noise[i],'/Y_',j,'.txt'))
#   }
# }

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
    CARD_obj = createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = simdata$obs.Y,
      spatial_location = spatial_location,
      ct.varname = "cellType",
      ct.select = unique(sc_meta$cellType),
      sample.varname = 'sampleInfo',
      minCountGene = 90,
      minCountSpot = 5)
    CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
    n.cell.types = 6
    I = simdata$obs.Y
    S.start <- CARD_obj@algorithm_matrix$B
    st_data = as.matrix(simdata$obs.Y)[match(rownames(S.start),rownames(as.matrix(simdata$obs.Y))),]
    STsisal_res = STsisal_pipeline_revised(st_data,S.start,6)
    STsisal_res$estProp = GetPropAligned(STsisal_res$estProp,simdata$true.prop,6)
    j_STsisal[[j]] = getMetric(simdata$true.prop,STsisal_res$estProp)
  }
  #rmse
  i_STsisal_rmse[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$rmse_vec))

  i_STsisal_corr[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$corr_ct))
  #cor_sample
  i_STsisal_cor_vec[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$cor_vec))
  #jsd
  i_STsisal_jsd[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$jsd_vec))
  #MAE
  i_STsisal_MAE[[i]] = sapply(1:20, function(i)mean(j_STsisal[[i]]$MAE_ct))
}

i_STsisal_rmse_revised = i_STsisal_rmse

i_STsisal_corr_revised = i_STsisal_corr
#cor_sample
i_STsisal_cor_vec_revised = i_STsisal_cor_vec
#jsd
i_STsisal_jsd_revised = i_STsisal_jsd
#MAE
i_STsisal_MAE_revised = i_STsisal_MAE


load('/Users/fuyinghao/Documents/STsisal/exp_res/sim_res.RData')
save(
  i_STsisal_rmse,
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
  i_STdeconvolve_MAE,i_STsisal_rmse_revised,

  i_STsisal_corr_revised,
  #cor_sample
  i_STsisal_cor_vec_revised,
  #jsd
  i_STsisal_jsd_revised,
  i_STsisal_MAE_revised ,
  file = '/Users/fuyinghao/Documents/STsisal/Revision/res/sim_res.RData'
)


library(ggplot2)
load('/Users/fuyinghao/Documents/STsisal/Revision/res/sim_res.RData')
color = RColorBrewer::brewer.pal(5,'Set3')[1:3]
##rmse
rmse =  rbind(STsisal = cbind(unlist(i_STsisal_rmse),rep('STsisal',60)),
              STsisal_revised = cbind(unlist(i_STsisal_rmse_revised),rep('STsisal_revised',60)),
              STdeconvolve = cbind(unlist(i_STdeconvolve_rmse),rep('STdeconvolve',60)),
              CARD = cbind(unlist(i_CARD_rmse),rep('CARD',60)),
              RCTD = cbind(unlist(i_RCTD_rmse),rep('RCTD',60)))
rmse = cbind(rmse,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),5))
colnames(rmse) = c('RMSE','Method','Heterogenity rate')
rmse = as.data.frame(rmse)
rmse$Method = as.factor(rmse$Method)
rmse$`Heterogenity rate` = as.factor(rmse$`Heterogenity rate`)
rmse$RMSE = as.numeric(rmse$RMSE)
p_rmse = ggplot(rmse, aes(x=Method, y=RMSE, fill=`Heterogenity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "RMSE")+
  theme(legend.position = "none")

##corr
corr =  rbind(STsisal = cbind(unlist(i_STsisal_corr),rep('STsisal',60)),
              STsisal_revised = cbind(unlist(i_STsisal_corr_revised),rep('STsisal_revised',60)),
              STdeconvolve = cbind(unlist(i_STdeconvolve_corr),rep('STdeconvolve',60)),
              CARD = cbind(unlist(i_CARD_corr),rep('CARD',60)),
              RCTD = cbind(unlist(i_RCTD_corr),rep('RCTD',60)))
corr = cbind(corr,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),5))
colnames(corr) = c('Corr','Method','Heterogenity rate')
corr = as.data.frame(corr)
corr$Method = as.factor(corr$Method)
corr$`Heterogenity rate` = as.factor(corr$`Heterogenity rate`)
corr$Corr = as.numeric(corr$Corr)
p_corr = ggplot(corr, aes(x=Method, y=Corr, fill=`Heterogenity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "Corr")+
  theme(legend.position = "none")

##jsd
jsd =  rbind(STsisal = cbind(unlist(i_STsisal_jsd),rep('STsisal',60)),
             STsisal_revised = cbind(unlist(i_STsisal_jsd_revised),rep('STsisal_revised',60)),
             STdeconvolve = cbind(unlist(i_STdeconvolve_jsd),rep('STdeconvolve',60)),
             CARD = cbind(unlist(i_CARD_jsd),rep('CARD',60)),
             RCTD = cbind(unlist(i_RCTD_jsd),rep('RCTD',60)))
jsd = cbind(jsd,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),5))
colnames(jsd) = c('Jsd','Method','Heterogenity rate')
jsd = as.data.frame(jsd)
jsd$Method = as.factor(jsd$Method)
jsd$`Heterogenity rate` = as.factor(jsd$`Heterogenity rate`)
jsd$Jsd = as.numeric(jsd$Jsd)
p_jsd = ggplot(jsd, aes(x=Method, y=Jsd, fill=`Heterogenity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "JSD")+
  theme(legend.position = "none")

##MAE
MAE =  rbind(STsisal = cbind(unlist(i_STsisal_MAE),rep('STsisal',60)),
             STsisal_revised = cbind(unlist(i_STsisal_MAE_revised),rep('STsisal_revised',60)),
             STdeconvolve = cbind(unlist(i_STdeconvolve_MAE),rep('STdeconvolve',60)),
             CARD = cbind(unlist(i_CARD_MAE),rep('CARD',60)),
             RCTD = cbind(unlist(i_RCTD_MAE),rep('RCTD',60)))
MAE = cbind(MAE,noise = rep(c(rep(0,20),rep(0.3,20),rep(0.6,20)),5))
colnames(MAE) = c('MAE','Method','Heterogenity rate')
MAE = as.data.frame(MAE)
MAE$Method = as.factor(MAE$Method)
MAE$`Heterogenity rate` = as.factor(MAE$`Heterogenity rate`)
MAE$MAE = as.numeric(MAE$MAE)
p_MAE = ggplot(MAE, aes(x=Method, y=MAE, fill=`Heterogenity rate`)) + 
  geom_boxplot()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_brewer(palette="BuPu")+labs(title = "MAE")+
  theme(legend.position = "none")


library(cowplot)

# 设置图像大小和图例位置
plot_size <- 5  # 图像大小为 5 英寸
legend_position <- "top"  # 图例位置设为上方

adjust_theme <- function(p) {
  p + theme(
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),  # 加大加粗X轴标题
    axis.title.y = element_text(size = 14, face = "bold"),  # 加大加粗Y轴标题
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),  # 旋转X轴刻度文字
    axis.text.y = element_text(size = 12, face = "bold"),  # 加大加粗Y轴刻度文字
    # legend.position = "top",
    legend.title = element_text(size = 12),  # 可以根据需要调整图例标题大小
    legend.text = element_text(size = 10)    # 可以根据需要调整图例文本大小
  )
}

# 应用adjust_theme函数调整每个图形
p_rmse <- adjust_theme(p_rmse)
p_corr <- adjust_theme(p_corr)
p_jsd <- adjust_theme(p_jsd)
p_MAE <- adjust_theme(p_MAE)


# 组合图像
plot_grid(p_rmse, p_corr, p_jsd, p_MAE, nrow = 2, ncol = 2, align = "hv") +
  theme(legend.position = legend_position) +
  cowplot::get_legend(p_rmse)



# 保存图像
ggsave('/Users/fuyinghao/Documents/STsisal/Revision/figure/M_revised.pdf', width = plot_size * 2, height = plot_size*2)

