rm(list = ls())
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/Figure4A_layer_annote.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/PDAC_A_metadata.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/PDAC_A_sc.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/PDAC_A_st_count.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/PDAC_A_st_location.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/PDAC_B_metadata.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/PDAC_B_sc.RData')
sc_meta_A = data.frame(cellID = colnames(sc_data_A),cellType = dat_A_meta,sampleInfo = rep('sampleInfo',length(dat_A_meta)))
sc_meta_B = data.frame(cellID = colnames(sc_data_B),cellType = dat_B_meta,sampleInfo = rep('sampleInfo',length(dat_B_meta)))
rownames(sc_meta_A) = colnames(sc_data_A)
rownames(sc_meta_B) = colnames(sc_data_B)
sc_data_A = as.matrix(sc_data_A)
sc_data_B = as.matrix(sc_data_B)
rownames(spatial_location) = colnames(st_count)
spatial_location = as.data.frame(spatial_location)
colnames(spatial_location) = c('x','y')
st_count = as.matrix(st_count)
####proper reference
source('/Users/fuyinghao/Documents/STsisal/pipeline/STsisal_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/pipeline/Deconvolution_pipeline.R')
##STsisal deconvolve
STsisal_res = STsisal_pipeline(st_count,20,n_marker = 500,TotalIter = 10)
##CARD deconvolve
CARD_res = CARD_pipeline(sc_data_A,sc_meta_A,st_count,spatial_location)
##RCTD deconvolve
RCTD_res = RCTD_pipeline(sc_data_A,sc_meta_A,st_count,spatial_location)
##Spotlight
# SPOTlight_res = SPOTlight_pipeline(sc_data_A,sc_meta_A,st_count,spatial_location)
##STdeconvolve
STdeconvolve_res = STdeconvolve_pipeline(st_count,20)
save(STsisal_res,CARD_res,RCTD_res,SPOTlight_res,STdeconvolve_res,file = '/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/A_res.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/A_res.RData')
colnames(RCTD_res) = c('Acinar cells','Cancer clone A cells','Cancer clone B cells','Ductal high hypoxic cells','Ductal centroacinar cells','Ductal antigen presenting cells','Ductal terminal cells','Endocrine cells','Endothelial cells',
                       'Fibroblasts','Microphage A cells','Microphage B cells','Mast cells',
                             'Myeloid dendritic A cells','Myeloid dendritic B cells','Monocytes','Plasmacytoid dendritic cells','Red blood cells',' T and natural killer (NK) cells','Tuft cells')

source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
location = spatial_location
colors_20 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#393b79", "#637939", "#8c6d31", "#bd9e39", "#ad494a", "#e7ba52", "#946cbb", "#1f77b4", "#aec7e8", "#ffbb78")
p1 = spatial_pieplot(STsisal_res$estProp,location,size = 0.5,colors_20)
p2 = spatial_pieplot(STdeconvolve_res,location,size = 0.5,colors_20)
p3 = spatial_pieplot(RCTD_res,location,size = 0.5,colors_20)
p4 = spatial_pieplot(CARD_res,location,size = 0.5,colors_20)
p5 = spatial_pieplot(SPOTlight_res,location,size = 0.5,colors_20)
cowplot::plot_grid(p1,p2,p3,p4,p5,ncol = 5,nrow = 1, hgap = 2, vgap = 2)
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC_pro_comparison.pdf',width = 10*5,height = 10*1,limitsize = FALSE)


# delete_zero = function(input){
#   id = which(colSums(input)==0)
#   input[,id] = 0.01
#   input = input/rowSums(input)
#   return(input)
# }
# plot_heatmap = function(input,reference,K){
#   library(ggcorrplot)
#   input = GetPropAligned(input,reference,K)
#   corMat = cor(input,reference,use="pairwise.complete.obs")
#   prop_cor = rep(0,K)
#   tmpmat=corMat
#   tmpmat[is.na(tmpmat)] = rnorm(K,-1,0.01)
#   p = ggcorrplot(tmpmat,lab=T)
#                  #colors = scale_fill_brewer(palette="BuPu"))
#   return(p)
# }
# 
# STsisal_res$estProp = delete_zero(STsisal_res$estProp)
# h_pro = plot_heatmap(CARD_res,STsisal_res$estProp,20)
# pdf('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC_pro_corr.pdf',width = 10,height = 10)
# print(h_pro)
# dev.off()


####improper reference
##STsisal deconvolve
STsisal_impro_res = STsisal_pipeline(st_count,13,n_marker = 500,TotalIter = 15)
##STdeconvolve
STdeconvolve_impro_res = STdeconvolve_pipeline(st_count,13)
colnames(STdeconvolve_impro_res) = colnames(STsisal_res$estProp)[1:13]
##CARD deconvolve
CARD_impro_res = CARD_pipeline(sc_data_B,sc_meta_B,st_count,spatial_location)
##RCTD deconvolve
RCTD_impro_res = RCTD_pipeline(sc_data_B,sc_meta_B,st_count,spatial_location)
##Spotlight
SPOTlight_impro_res = SPOTlight_pipeline(sc_data_B,sc_meta_B,st_count,spatial_location)
save(STsisal_impro_res,STdeconvolve_impro_res,CARD_impro_res,RCTD_impro_res,SPOTlight_impro_res,file = '/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/B_res.RData')

load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/B_res.RData')
colors =  c(RColorBrewer::brewer.pal(12,'Set3')[1:12],RColorBrewer::brewer.pal(3,'Set1')[1])
colors_13 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#393b79", "#637939", "#8c6d31")
source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
location = spatial_location
colnames(RCTD_impro_res) = c('Acinar cells','Cancer clone A cells','Ductal centroacinar cells','Ductal antigen presenting cells','Ductal terminal cells','Endocrine cells','Endothelial cells','Microphage cells','Mast cells',
                             'Myeloid dendritic cells','Monocytes','Red blood cells','Tuft cells')
p1 = spatial_pieplot(STsisal_impro_res$estProp,location,size = 0.5,colors_13)
p2 = spatial_pieplot(STdeconvolve_impro_res,location,size = 0.5,colors_13)
p3 = spatial_pieplot(RCTD_impro_res,location,size = 0.5,colors_13)
p4 = spatial_pieplot(CARD_impro_res,location,size = 0.5,colors_13)
p5 = spatial_pieplot(SPOTlight_impro_res,location,size = 0.5,colors_13)
cowplot::plot_grid(p3,ncol = 1,nrow = 1, hgap = 5, vgap = 5)
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/part_impro_comparison.pdf',width = 10*5,height = 10*1,limitsize = FALSE)


STsisal_impro_res$estProp = delete_zero(STsisal_impro_res$estProp)
h_impro = plot_heatmap(CARD_impro_res,STsisal_impro_res$estProp,13)
pdf('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC_impro_corr.pdf',width = 10,height = 10)
print(h_impro)
dev.off()

input = GetPropAligned(CARD_impro_res,STsisal_impro_res$estProp,13)
corMat_impro = cor(input,STsisal_impro_res$estProp,use="pairwise.complete.obs")

input = GetPropAligned(CARD_res,STsisal_res$estProp,20)
corMat_pro = cor(input,STsisal_res$estProp,use="pairwise.complete.obs")
mean(diag(corMat_pro))
mean(diag(corMat_impro))

load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/A_res.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/B_res.RData')
####Comparison
get_reference = function(sc_count,sc_meta){
  reference = matrix(NA,ncol = length(unique(sc_meta$cellType)),nrow = nrow(sc_count))
  ct = unique(sc_meta$cellType)
  rownames(reference) = rownames(sc_count)
  colnames(reference) = ct
  for (i in 1:ncol(reference)){
    reference[,i] = rowMeans(sc_count[,sc_meta[which(sc_meta$cellType==ct[i]),'cellID']])
  }
  return(reference)
}
reference_A = get_reference(sc_data_A,sc_meta_A)
reference_B = get_reference(sc_data_B,sc_meta_B)

sel_marker_A = unique(unlist(STsisal_res$selMarker))
sel_marker_B = unique(unlist(STsisal_impro_res$selMarker))
reduce_A = reference_A[sel_marker_A,]
reduce_B = reference_B[sel_marker_B,]


reduce_A = as.data.frame(t(rbind(reduce_A,group = factor(1:ncol(reduce_A)))))
reduce_B = as.data.frame(t(rbind(reduce_B,group = factor(1:ncol(reduce_B)))))

p_A = apply(reduce_A[,1:(ncol(reduce_A)-1)], 2, function(i){summary(aov(i~group,data = reduce_A))[[1]]$Pr[[1]]})
p_B = apply(reduce_B[,1:(ncol(reduce_B)-1)], 2, function(i){summary(aov(i~group,data = reduce_B))[[1]]$Pr[[1]]})

g = c(rep('Matched reference', length(p_A)),rep('Mismatched reference',length(p_B)))
p = c(p_A,p_B)
p = as.data.frame(cbind(P_value = p, group = factor(g)))
p$group[p$group==1] = 'Matched reference'
p$group[p$group==2] = 'Mismatched reference'
p$group = as.factor(p$group)
colnames(p) = c('P_value','Reference_class')
colors <- brewer.pal(3, "Pastel1")[1:2]
ggplot2::ggplot(p, aes(x=Reference_class, y=P_value, fill=Reference_class)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  labs(x = "Reference class", y = "P value") +
  scale_x_discrete(labels = c('Matched reference', 'Mismatched reference')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 16),axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16),legend.position = "none")

ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC_p_value.pdf',height = 8,width = 8)


load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/PDAC/Figure4A_layer_annote.RData')
layer_manual_PDAC
colors = c("#9467bd", "#8c564b", "#e377c2", "#17becf")

ggplot(layer_manual_PDAC,aes(x = x, y = y,
                   fill = factor(Region)))+
  geom_point(shape = 21,size = 10,color = 'white')+
  scale_fill_manual(values = colors)+
  theme(plot.margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm"), 
        panel.background = element_blank(), plot.background = element_blank(), 
        panel.border = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title = element_blank(), legend.title = element_text(size = 16, 
                                                                  face = "bold"), legend.text = element_text(size = 15), 
        legend.key = element_rect(colour = "transparent", 
                                  fill = "white"), legend.key.size = unit(0.45, 
                                                                          "cm"), strip.text = element_text(size = 16, face = "bold"), 
        legend.position = "bottom" )+ guides(fill = guide_legend(title = "Cell Type", nrow = 2, ncol = 2))+theme(plot.title = element_text(size = 30,vjust = 0.3, hjust = 0.5,face = "bold"))
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC_ground_truth.pdf',width = 10,height = 10)



for (i in 1:length(sel_marker_A)){
  id = which(rownames(st_count)==sel_marker_A[i])
  new_dat = data.frame(value = st_count[id,], loc_x = location[, 1], loc_y = location[, 2])
  pdf(paste0('/Users/fuyinghao/Documents/kernel/res_notcombn/',i,'.pdf'),width=10,height=10)
  print(ggplot(new_dat,aes(x = loc_x, y = loc_y,
                           fill = value))+
          geom_point(shape = 21,size = 10)+
          scale_fill_continuous(type = "viridis")  +
          theme_bw())
  dev.off()
}

legend_key_size <- unit(1.5, "cm")  # 定义较大的标尺大小
legend_text_size <- 18  # 定义较小的图例文本字体大小

for (i in 1:ncol(STsisal_res$estProp)) {
  new_dat <- data.frame(value = STsisal_res$estProp[, i], loc_x = location[, 1], loc_y = location[, 2])
  colnames(new_dat) <- c('Proportion', 'loc_x', 'loc_y')
  pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC/', i, '.pdf'), width = 10, height = 10)
  print(ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Proportion)) +
          geom_point(shape = 21, size = 10.2, color = 'white') +
          scale_fill_continuous(type = "viridis") +
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
library(viridis)  # 加载viridis库，用于更改调色板

# 设置新的调色板
new_palette <- inferno(n = 100)
st_gene = c("TM4SF1"   ,"APOL1"    ,"TFF1"     ,"TFF3"    , "AQP5"     ,"SLC38A10")
min_max = function(dat){
  dat = (dat-min(dat))/(max(dat)-min(dat))
  return(dat)
}
for (i in st_gene) {
  new_dat <- data.frame(value = min_max(st_count[which(rownames(st_count) == i), ]), loc_x = location[, 1], loc_y = location[, 2])
  colnames(new_dat) <- c('Count', 'loc_x', 'loc_y')
  pdf(paste0('/Users/fuyinghao/Documents/STsisal/exp_res/PDAC_gene/', i, '.pdf'), width = 10, height = 10)
  print(ggplot(new_dat, aes(x = loc_x, y = loc_y, fill = Count)) +
          geom_point(shape = 21, size = 10.2, color = 'white') +
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




#####cell type heat map
pro_marker = unlist(STsisal_res$selMarker)
impro_marker = unlist(STsisal_impro_res$selMarker)
length(pro_marker)
length(impro_marker)
intersect(pro_marker,impro_marker)
colors <- brewer.pal(3, "Pastel1")[1:2]



####get optimal estimated cell type number 
# load('/Users/fuyinghao/Documents/STsisal/exp_res/AIC_PDAC.RData')
# Kres = getCellNumber_st(st_count, 10:20)
st_count = clean_count(st_count)
Kres = getCellNumber_st(st_count,10:25)
K = Kres$bestK
Kres$allAIC[1:11]

st_count = clean_count(st_count)
Kres_add = getCellNumber_st(st_count,15)
K = Kres$bestK
Kres$allAIC
Kres_add$allAIC
allAIC = c(Kres$allAIC,Kres_add$allAIC[-1])
save(allAIC,file = '/Users/fuyinghao/Documents/STsisal/Revision/res/AIC.RData')
dat = data.frame(cell_type_number = 10:30,AIC = allAIC)
p = ggplot(dat, aes(x = cell_type_number, y = AIC)) +
  geom_line(color = colors[1], size = 1.5) +
  geom_point(shape = 19, size = 5, color = colors[2]) +
  labs(x = "Cell Type Number", y = "AIC", size = 16) +  # 设置字体大小为 12
  theme(axis.text = element_text(size = 16, face = "bold"),axis.title.y = element_text(size = 18, face = "bold"),axis.title.x = element_text(size = 18, face = "bold"))+
  scale_x_continuous(breaks = seq(min(dat$cell_type_number), max(dat$cell_type_number), 1), 
                     minor_breaks = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank())
# save(Kres,file = '/Users/fuyinghao/Documents/STsisal/exp_res/AIC_PDAC.RData')
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/AIC_PDAC.pdf', plot = p, width = 10, height = 6, dpi = 1000)


p = ggplot(dat, aes(x = cell_type_number, y = AIC)) +
  geom_line(color = colors[1], size = 1.5) +  # Line for AIC trends
  geom_point(shape = 19, size = 5, color = colors[2]) +  # Points with a color distinction
  scale_color_manual(values = colors[2], guide = FALSE) +  # Customize point colors, no legend for colors
  labs(
    x = "Cell Type Number",
    y = "AIC"
    # title = "AIC Across Different Cell Type Numbers"  # Adding a title for clarity
  ) +
  theme_bw(base_size = 16) +  # Use theme_bw as the base theme with adjusted base size
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold")  # Center the plot title
  ) +
  scale_x_continuous(
    breaks = seq(min(dat$cell_type_number), max(dat$cell_type_number), 1), 
    minor_breaks = NULL
  )



delete_zero = function(input){
  id = which(colSums(input)==0)
  input[,id] = 0.01
  input = input/rowSums(input)
  return(input)
}
STsisal_res$estProp = delete_zero(STsisal_res$estProp)
out_all <- mycsfit(STsisal_res$estProp, t(st_count))
prof <- t(out_all$ghat)
rownames(prof) <- rownames(Y_raw)
selProf <- prof[unlist(selMarker), ]
labres <- GetCorRes(selProf, knowRefAll = knowRef)
colnames(estProp) <- labres$assignLabel






