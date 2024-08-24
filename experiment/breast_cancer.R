rm(list = ls())
source('/Users/fuyinghao/Documents/STsisal/pipeline/STsisal_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/pipeline/Deconvolution_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/function/getMetric.R')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P2/breast cancer_Endecon/breast.sc.cell.label.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P2/breast cancer_Endecon/breast.sc.ref.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P2/breast cancer_Endecon/breast.spot.annotation.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P2/breast cancer_Endecon/breast.st.loc.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P2/breast cancer_Endecon/breast.st.RData')

st_count = breast.st
spatial_location = breast.st.loc
sc_count = breast.sc.ref
sc_meta = breast.sc.cell.label
sc_meta = as.data.frame(sc_meta)

colnames(sc_meta) = 'cellType'
spot_names <- paste0("spot", 1:nrow(spatial_location))
colnames(st_count) <- spot_names
rownames(spatial_location) <- spot_names
colnames(spatial_location) <- c("x", "y")

sampleInfo <- rep("sample1", ncol(sc_count))
names(sampleInfo) <- colnames(sc_count)
sc_meta <- data.frame(cellID = colnames(sc_count), cellType = sc_meta, sampleInfo = sampleInfo)
rownames(sc_meta) <- colnames(sc_count)

##STsisal deconvolve
STsisal_res = STsisal_pipeline(st_count,8)
##CARD deconvolve
CARD_res = CARD_pipeline(sc_count,sc_meta,st_count,spatial_location)
##RCTD deconvolve
RCTD_res = RCTD_pipeline(sc_count,sc_meta,st_count,spatial_location)
##STdeconvolve deconvolve
STdeconvolve_res = STdeconvolve_pipeline(st_count,8)





# source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
# location = spatial_location
# colors = RColorBrewer::brewer.pal(7,'Set1')[1:8]
# p1 = spatial_pieplot(STsisal_res$estProp,location,colors)
# p2 = spatial_pieplot(STdeconvolve_res ,location,colors)
# p3 = spatial_pieplot(RCTD_res,location,colors)
# p4 = spatial_pieplot(CARD_res,location,colors)
# cowplot::plot_grid(p1,p2,p3,p4,ncol = 2,nrow = 2,labels = c('STsisal', 'Stdeconvolve','RCTD','CARD'), align = "hv")
# ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/brease_cancer_comparison.pdf',width = 8*2,height = 8*2)
# 
# max_type = function(prop){
#   idx = apply(prop,1,function(i) which(i == max(i)))
#   for (i in 1:length(idx)) {
#     tmp = rep(0,ncol(prop))
#     tmp[idx[i]] = 1
#     prop[i,] = tmp
#   }
#   return(prop)
# }
# 
# STsisal_res$estProp_max = max_type(STsisal_res$estProp)
# CARD_res_max = max_type(CARD_res)
# RCTD_res_max = max_type(RCTD_res)
# STdeconvolve_res_max = max_type(STdeconvolve_res)
# source('/Users/fuyinghao/Downloads/TOAST/R/GetPropAligned.R')
# 
# p1_n = spatial_pieplot(STsisal_res$estProp_max,location,colors)
# p2_n = spatial_pieplot(STdeconvolve_res_max ,location,colors)
# p3_n = spatial_pieplot(RCTD_res_max,location,colors)
# p4_n = spatial_pieplot(CARD_res_max,location,colors)
# cowplot::plot_grid(p1_n,p2_n,p3_n,p4_n,ncol = 4,nrow = 1,labels = c('STsisal', 'STdeconvolve','RCTD','CARD'),align = "hv",label_size = 30,label_x = 0.4)
# ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/breast_cancer_comparison_maxtype.pdf',width = 8*4,height = 8*1)


library(gridExtra)
col_names_ori <- colnames(CARD_res)
# col_names_ori[8] = 'Epithelial'
legend_key_size <- unit(1.5, "cm")  # 定义较大的标尺大小
legend_text_size <- 18  # 定义较小的图例文本字体大小
plots_list <- list()
for (i in 1:ncol(CARD_res)) {
  proportion_range <- range(CARD_res[, i])
  new_dat <- data.frame(value = CARD_res[, i], loc_x = spatial_location[, 1], loc_y = spatial_location[, 2])
  colnames(new_dat) <- c('Proportion', 'loc_x', 'loc_y')
  p = ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Proportion)) +
          geom_point(shape = 21, size = 12, color = 'white') +
          scale_fill_continuous(type = "viridis", limits = proportion_range) +  # 设置图例范围
          ggtitle(col_names_ori[i]) +  # Add column name as the title
          theme(plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"), 
                panel.background = element_blank(), plot.background = element_blank(), 
                panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
                axis.title = element_blank(), legend.title = element_text(size = 20, face = "bold"), 
                legend.text = element_text(size = legend_text_size), 
                legend.key = element_rect(colour = "transparent", fill = "white"), 
                legend.key.size = legend_key_size, 
                legend.spacing.x = unit(0.5, "cm"), 
                legend.spacing.y = unit(0.5, "cm"), 
                strip.text = element_text(size = 16, face = "bold"), 
                legend.position = "bottom") +
          theme(plot.title = element_text(size = 30, vjust = 0.3, hjust = 0.5, face = "bold"))
  plots_list[[i]] <- p
}

# 保存到PDF
pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/breast_CARD.pdf'), width = 9*4, height = 9*2)
grid.arrange(grobs = plots_list, ncol = 4, nrow = 2)
dev.off()


col_names <- colnames(RCTD_res)
RCTD_res = RCTD_res[,match(col_names_ori,col_names)]
col_names <- colnames(RCTD_res)
# col_names[8] = 'Epithelial'
legend_key_size <- unit(1.5, "cm")  # 定义较大的标尺大小
legend_text_size <- 18  # 定义较小的图例文本字体大小
plots_list <- list()
for (i in 1:ncol(RCTD_res)) {
  proportion_range <- range(RCTD_res[, i])
  new_dat <- data.frame(value = RCTD_res[, i], loc_x = spatial_location[, 1], loc_y = spatial_location[, 2])
  colnames(new_dat) <- c('Proportion', 'loc_x', 'loc_y')
  p = ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Proportion)) +
    geom_point(shape = 21, size = 12, color = 'white') +
    scale_fill_continuous(type = "viridis", limits = proportion_range) +  # 设置图例范围
    ggtitle(col_names[i]) +  # Add column name as the title
    theme(plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"), 
          panel.background = element_blank(), plot.background = element_blank(), 
          panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.title = element_blank(), legend.title = element_text(size = 20, face = "bold"), 
          legend.text = element_text(size = legend_text_size), 
          legend.key = element_rect(colour = "transparent", fill = "white"), 
          legend.key.size = legend_key_size, 
          legend.spacing.x = unit(0.5, "cm"), 
          legend.spacing.y = unit(0.5, "cm"), 
          strip.text = element_text(size = 16, face = "bold"), 
          legend.position = "bottom") +
    theme(plot.title = element_text(size = 30, vjust = 0.3, hjust = 0.5, face = "bold"))
  plots_list[[i]] <- p
}

# 保存到PDF
pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/breast_STsisal.pdf'), width = 9*4, height = 9*2)
grid.arrange(grobs = plots_list, ncol = 4, nrow = 2)
dev.off()
col_names = colnames(STsisal_res$estProp)
plots_list <- list()
for (i in 1:ncol(STsisal_res$estProp)) {
  proportion_range <- range(STsisal_res$estProp[, i])
  new_dat <- data.frame(value = STsisal_res$estProp[, i], loc_x = spatial_location[, 1], loc_y = spatial_location[, 2])
  colnames(new_dat) <- c('Proportion', 'loc_x', 'loc_y')
  p = ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Proportion)) +
    geom_point(shape = 21, size = 12, color = 'white') +
    scale_fill_continuous(type = "viridis", limits = proportion_range) +  # 设置图例范围
    ggtitle(col_names[i]) +  # Add column name as the title
    theme(plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"), 
          panel.background = element_blank(), plot.background = element_blank(), 
          panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.title = element_blank(), legend.title = element_text(size = 20, face = "bold"), 
          legend.text = element_text(size = legend_text_size), 
          legend.key = element_rect(colour = "transparent", fill = "white"), 
          legend.key.size = legend_key_size, 
          legend.spacing.x = unit(0.5, "cm"), 
          legend.spacing.y = unit(0.5, "cm"), 
          strip.text = element_text(size = 16, face = "bold"), 
          legend.position = "bottom") +
    theme(plot.title = element_text(size = 30, vjust = 0.3, hjust = 0.5, face = "bold"))
  plots_list[[i]] <- p
}

# 保存到PDF
pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/breast_STsisal.pdf'), width = 9*4, height = 9*2)
grid.arrange(grobs = plots_list, ncol = 4, nrow = 2)
dev.off()


load('/Users/fuyinghao/Documents/STsisal/exp_res/STsisal_breast.RData')

