rm(list = ls())
library(ape)
source('/Users/fuyinghao/Documents/STsisal/pipeline/STsisal_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/pipeline/Deconvolution_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/function/getMetric.R')
source('/Users/fuyinghao/Documents/STsisal/function/sim_functions.R')
source('/Users/fuyinghao/Downloads/TOAST/R/getCellNumber.R')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/MOB/Rep12_MOB_count_matrix-1.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/MOB/MOB.dge.sceset.RData')
library(SingleCellExperiment)
library(pbmcapply)
ct.varname = "cellType"
sample.varname = "sampleInfo"
sc_count = assays(sce)$counts
sc_meta = colData(sce)
spatial_count = MOB_raw
location <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(MOB_raw),split="x"),"[",1)),y=as.numeric(sapply(strsplit(colnames(MOB_raw),split="x"),"[",2)))
rownames(location) = colnames(MOB_raw)
spatial_location = location
table(sc_meta$cellType)
# idx = sc_meta$cellType!='EPL-IN'
sc_count = sc_count[,match(rownames(sc_meta),colnames(sc_count))]
# sc_meta = sc_meta[idx,]
# sc_count = sc_count[,idx]
sc_meta$cellType = as.factor(sc_meta$cellType)
levels(sc_meta$cellType)[3] = "M-TC"

library(STdeconvolve)
library(TOAST)
knowRefAll = sc_count
colnames(knowRefAll) = sc_meta$cellType[match(colnames(sc_count),rownames(sc_meta))]
##CARD deconvolve
CARD_res = CARD_pipeline(sc_count,sc_meta,spatial_count,spatial_location)


colors = RColorBrewer::brewer.pal(7,'Set1')[1:4]
# load('/Users/fuyinghao/Documents/STsisal/exp_res/mob_res.RData')
# STdeconvolve deconvolve
STdeconvolve_res = STdeconvolve_pipeline(as.matrix(spatial_count),5)
load('/Users/fuyinghao/Documents/STsisal/exp_res/mob_res.RData')
##STsisal deconvolve
STsisal_res = STsisal_pipeline(as.matrix(spatial_count),5,n_marker = 500)

##RCTD deconvolve
RCTD_res = RCTD_pipeline(sc_count,sc_meta,spatial_count,spatial_location)
##SPOTlight devonvolve
# SPOTlight_res = SPOTlight_pipeline(sc_count,sc_meta,spatial_count,spatial_location)
# save(STsisal_res,STdeconvolve_res,CARD_res,RCTD_res,file= '/Users/fuyinghao/Documents/STsisal/exp_res/mob_res.RData')
colors = RColorBrewer::brewer.pal(7,'Set1')[1:4]

source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
colors = RColorBrewer::brewer.pal(7,'Set1')[1:5]
p1 = spatial_pieplot(STsisal_res$estProp,size= 0.45,location,colors)
p2 = spatial_pieplot(STdeconvolve_res ,location,size = 0.45,colors)
p3 = spatial_pieplot(RCTD_res,location,size = 0.45,colors)
p4 = spatial_pieplot(CARD_res,location,size = 0.45,colors)
#p4 = spatial_pieplot(SPOTlight_res,size= 0.45,location,colors)
cowplot::plot_grid(p1,p2,p3,p4,ncol = 4,nrow = 1)
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/MOB_comparison.pdf',width = 8*4,height = 8*1)

p = spatial_pieplot(STsisal_res$estProp,size= 0.45,location,colors)
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/MOB_STsisal.pdf',width = 10*1,height = 10*1)

max_type = function(prop){
  idx = apply(prop,1,function(i) which(i == max(i)))
  for (i in 1:length(idx)) {
    tmp = rep(0,ncol(prop))
    tmp[idx[i]] = 1
    prop[i,] = tmp
  }
  return(prop)
}

STsisal_res$estProp_max = max_type(STsisal_res$estProp)
CARD_res_max = max_type(CARD_res)
RCTD_res_max = max_type(RCTD_res)
STdeconvolve_res_max = max_type(STdeconvolve_res)
# SPOTlight_res_max = max_type(SPOTlight_res)

p_n = spatial_pieplot(STsisal_res$estProp_max,size= 0.45,location,colors)
p2_n = spatial_pieplot(STdeconvolve_res_max ,location,size = 0.45,colors)
p3_n = spatial_pieplot(RCTD_res_max,location,size = 0.45,colors[c(2,3,4,5)])
p4_n = spatial_pieplot(CARD_res_max,location,size = 0.45,colors[c(2,3,4,5)])
# p5_n = spatial_pieplot(SPOTlight_res_max,location,size = 0.45,colors[c(1,4)])
cowplot::plot_grid(p_n,p2_n,p3_n,p4_n,ncol = 4,nrow = 1)
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/MOB_comparison_maxtype.pdf',width = 8*4,height = 8*1)





load('/Users/fuyinghao/Documents/STsisal/exp_res/mob_res.RData')
STsisal_res$selMarker
library(viridis)
legend_key_size <- unit(1.5, "cm")  # 定义较大的标尺大小
legend_text_size <- 18  # 定义较小的图例文本字体大小
new_palette <- inferno(n = 100)
st_gene = STsisal_res$selMarker$`3`
min_max = function(dat){
  dat = (dat-min(dat))/(max(dat)-min(dat))
  return(dat)
}
st_count = spatial_count
for (i in st_gene) {
  value = min_max(st_count[which(rownames(st_count) == i), ])
  new_dat <- data.frame(value = value, loc_x = location[, 1], loc_y = location[, 2])
  colnames(new_dat) <- c('Count', 'loc_x', 'loc_y')
  pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/MOB_gene/', i, '.pdf'), width = 10, height = 10)
  print(ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Count)) +
          geom_point(shape = 21, size = 12, color = 'white') +
          scale_fill_viridis_c(option = "inferno") +  # 使用新的调色板
          theme(plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"),
                panel.background = element_blank(), plot.background = element_blank(),
                panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                axis.title = element_blank(), legend.title = element_text(size = 20, face = "bold"),
                legend.text = element_text(size = legend_text_size),
                legend.key = element_rect(colour = "transparent", fill = "white"),
                legend.key.size = legend_key_size,
                legend.spacing.x = unit(0.5, "cm"),  # 增加图例项之间的水平间距
                legend.spacing.y = unit(0.5, "cm"),  # 增加图例项之间的垂直间距
                strip.text = element_text(size = 16, face = "bold"),
                legend.position = "bottom") +
          theme(plot.title = element_text(size = 30, vjust = 0.3, hjust = 0.5, face = "bold")))
  dev.off()
}

# 获取Proportion值的范围
proportion_range <- range(STsisal_res$estProp)

for (i in 1:ncol(STsisal_res$estProp)) {
  new_dat <- data.frame(value = STsisal_res$estProp[, i], loc_x = location[, 1], loc_y = location[, 2])
  colnames(new_dat) <- c('Proportion', 'loc_x', 'loc_y')
  pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/MOB/', i, '.pdf'), width = 10, height = 10)
  print(ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Proportion)) +
          geom_point(shape = 21, size = 12, color = 'white') +
          scale_fill_continuous(type = "viridis", limits = proportion_range) +  # 设置图例范围
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
          theme(plot.title = element_text(size = 30, vjust = 0.3, hjust = 0.5, face = "bold")))
  dev.off()
}





# p = heatmap(t(STsisal_res$estProp), col = heat.colors(256), Rowv = FALSE, Colv = FALSE)
df <- data.frame(STsisal_res$estProp)
data = cor(df)
colnames(data) = c('V1','V2','V3','V4','V5')
rownames(data) = c('V1','V2','V3','V4','V5')
data1 <- melt(data)# 创建热图

ggplot(data = data1, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "BuPu", direction = 1, limits = c(-1, 1), na.value = NA) + # Set limits for correlation coefficients
  labs(x = "Cell Type", y = "Cell Type", fill = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18)
  )
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/MOB/heatmap.pdf',width = 10,height = 8)


unique(sc_meta$cellType)
df = as.data.frame(table(sce@colData$cellType))
colnames(df) = c('Celltype','Number')
df$Number = as.numeric(df$Number)
df$Celltype = as.factor(df$Celltype)
df = df[order(df$Number,decreasing = T),]
library(ggplot2)
library(extrafont)
color = RColorBrewer::brewer.pal(9,'BuPu')[3:7]
myLabel = as.vector(df$Celltype)   
myLabel = paste(myLabel, "(", round(df$Number / sum(df$Number) * 100, 2), "%)        ", sep = "")
ggplot(df, aes(x = "", y = Number,fill = Celltype)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  labs(x = "", y = "", title = "") + 
  coord_polar("y", start=0)+
  theme_void() +
  labs(title = "Cell type number of MOB", fill = "Category")+
  # geom_text(aes(label = myLabel), position = position_stack(vjust = 0.5), size = 4, color = "white") +
  scale_fill_manual(values=color,labels = myLabel)+
  theme(legend.position = "right", legend.title = element_text(size = 15,face = "bold"), 
        legend.text = element_text(size = 15,face = "bold"),
        plot.title = element_text(size = 20,vjust = 0.3, hjust = 0.5,face = "bold"))
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/ct_number_mob.pdf',width = 8,height = 8, dpi = 300)





