rm(list = ls())
# load('/Users/fuyinghao/Documents/STsisal/exp_res/mob_res.RData')
source('/Users/fuyinghao/Documents/STsisal/pipeline/STsisal_pipeline.R')
# source('/Users/fuyinghao/Documents/STsisal/Revision/function/STsisal_pipeline_revised.R')
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
sc_count = sc_count[,match(rownames(sc_meta),colnames(sc_count))]
sc_meta$cellType = as.factor(sc_meta$cellType)
levels(sc_meta$cellType)[3] = "M-TC"
idx = which(sc_meta$cellType=="EPL-IN")
library(Seurat)
mob <- CreateSeuratObject(counts = sc_count, project = "pbmc3k", min.cells = 3, min.features = 200)
mob <- FindVariableFeatures(mob, selection.method = "vst", nfeatures = 2000)
sc_count = as.matrix(GetAssayData(mob,slot = 'counts'))
sc_count = sc_count[,-which(colnames(sc_count)%in%rownames(sc_meta)[which(sc_meta$cellType=="EPL-IN")])]
sc_meta = sc_meta[-idx,]
sc_meta$cellType = droplevels(sc_meta$cellType)
knowRef =  clean_count(as.matrix(sc_count))
colnames(knowRef) = sc_meta$cellType[match(colnames(sc_count),rownames(sc_meta))]
Y_raw = clean_count(as.matrix(spatial_count))
STsisal_res = STsisal_pipeline(as.matrix(spatial_count),5,n_marker = 500)
estProp = STsisal_res$estProp
selMarker = STsisal_res$selMarker
out_all <- mycsfit(estProp, t(Y_raw))
prof <- t(out_all$ghat)
rownames(prof) <- rownames(Y_raw)
selProf <- prof[unlist(selMarker), ]
knowRefAll = knowRef
isc = intersect(rownames(selProf), rownames(knowRefAll))
selProf = selProf[isc, ]
knowRef <- knowRefAll[isc, ]
cormat <- cor(knowRef, selProf)
# nbmat <- MyNaiveBayes(selProf, knowRef)
summat <- cormat
summat_org <- summat
assignLabel <- rep("Unassigned", ncol(selProf))
for (i in 1:ncol(knowRefAll)) {
  thisidx <- which(summat == max(summat), arr.ind = TRUE)
  assignLabel[thisidx[2]] <- rownames(summat)[thisidx[1]]
  summat[thisidx[1], ] <- -1
  summat[, thisidx[2]] <- -1
}

colnames(estProp) <- assignLabel

estProp[,2] = estProp[,2]+estProp[,4]
estProp = estProp[,-c(4)]


source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
colors = RColorBrewer::brewer.pal(7,'Set1')[1:5]
p1 = spatial_pieplot(estProp,size= 0.45,location,colors)
max_type = function(prop){
  idx = apply(prop,1,function(i) which(i == max(i)))
  for (i in 1:length(idx)) {
    tmp = rep(0,ncol(prop))
    tmp[idx[i]] = 1
    prop[i,] = tmp
  }
  return(prop)
}

estProp_max = max_type(estProp)
p2 = spatial_pieplot(estProp_max,size= 0.45,location,colors[1:5])
cowplot::plot_grid(p1,p2,ncol = 2,nrow = 1)
ggsave('/Users/fuyinghao/Documents/STsisal/Revision/figure/MOB5_not_aligned.pdf', width = 20, height = 8)

label_matrix = as.matrix(apply(estProp_max,1,function(i) which(i==1)),ncol = 1)
label_matrix[,2] = as.numeric(label_matrix[,2])
spatial_location


##X#generate label matrix
# label_matrix=read.csv('/Users/fuyinghao/Documents/STsisal/Revision/res/MOB_true_label.csv')
# label_matrix = rbind(colnames(label_matrix),label_matrix)
# label_matrix$X1[1] = 1
# label_matrix = cbind(spatial_location,label = label_matrix$X1)
# # 如果你的label_matrix是一个matrix，你需要将它转换为data frame
# label_df <- as.data.frame(label_matrix)
# 
# # 绘制散点图，使用label作为颜色分组
# ggplot(label_df, aes(x = x, y = y, color = factor(label))) +
#   geom_point(size= 6) +  # 添加点
#   scale_color_brewer(palette = "Set1") +  # 使用预设的颜色方案
#   theme_minimal() +  # 使用简洁主题
#   labs(color = "Label")  # 添加图例标题
# save(label_matrix,file = '/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/MOB/True_label.RData')

library(mclust)
load('/Users/fuyinghao/Documents/STsisal/Revision/res/mob_STsisal_4.RData')
STsisal_label = apply(estProp_max,1,function(i) which(i==1))
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/MOB/True_label.RData')
ari_STsisal <- adjustedRandIndex(STsisal_label, label_matrix$label)
source('/Users/fuyinghao/Documents/STsisal/pipeline/Deconvolution_pipeline.R')
STdeconvolve_res = STdeconvolve_pipeline(as.matrix(spatial_count),4)
STdeconvolve_label = apply(max_type(STdeconvolve_res),1,function(i) which(i==1))
ari_STdeconvolve <- adjustedRandIndex(STdeconvolve_label, label_matrix$label)

CARD_res = CARD_pipeline(sc_count,sc_meta,spatial_count,spatial_location)
CARD_label = apply(max_type(CARD_res),1,function(i) which(i==1))
ari_CARD <- adjustedRandIndex(CARD_label, label_matrix$label)
RCTD_res = RCTD_pipeline(sc_count,sc_meta,spatial_count,spatial_location)
RCTD_label = apply(max_type(RCTD_res),1,function(i) which(i==1))
ari_RCTD <- adjustedRandIndex(RCTD_label, label_matrix$label)

# save(STsisal_res,file = '/Users/fuyinghao/Documents/STsisal/Revision/res/mob_STsisal_4.RData')

ari_values <- c(RCTD = ari_RCTD, CARD = ari_CARD, STdeconvolve = ari_STdeconvolve, STsisal = ari_STsisal)

# 将ARI值转换为数据框
ari_data <- data.frame(Method = names(ari_values), ARI = ari_values)

calculatePurity <- function(cluster_labels, true_labels) {
  # 确保聚类标签和真实标签长度相同
  if(length(cluster_labels) != length(true_labels)) {
    stop("Length of cluster labels and true labels must be the same")
  }
  
  # 组合聚类标签和真实标签
  combined_labels <- data.frame(cluster_labels, true_labels)
  
  # 计算每个聚类的纯度
  purity_sum <- 0
  clusters <- unique(combined_labels$cluster_labels)
  for(cluster in clusters) {
    # 在每个聚类中找到最常见的真实标签
    subset <- subset(combined_labels, cluster_labels == cluster)
    most_common_label <- names(sort(table(subset$true_labels), decreasing = TRUE))[1]
    # 统计与最常见真实标签一致的数据点数量
    correct_labels_count <- sum(subset$true_labels == most_common_label)
    # 累加到总和中
    purity_sum <- purity_sum + correct_labels_count
  }
  
  # 计算总纯度
  purity <- purity_sum / length(cluster_labels)
  
  return(purity)
}
purity_STsisal = calculatePurity(STsisal_label,label_matrix$label)
purity_STdeconvolve = calculatePurity(STdeconvolve_label,label_matrix$label)
purity_CARD = calculatePurity(CARD_label,label_matrix$label)
purity_RCTD = calculatePurity(RCTD_label,label_matrix$label)
purity_values <- c(RCTD = purity_RCTD, CARD = purity_CARD, STdeconvolve = purity_STdeconvolve, STsisal = purity_STsisal)

# 将ARI值转换为数据框
purity_data <- data.frame(Method = names(purity_values), Purity = purity_values)

ari_data$Metric <- 'ARI'
purity_data$Metric <- 'Purity'
colnames(ari_data)[2] = 'Value'
colnames(purity_data)[2] = 'Value'
combined_data <- rbind(ari_data, purity_data)
library(RColorBrewer)
method_colors <- brewer.pal(4, "Dark2")
method_colors <- c('#9FBA95', '#EDB8B0', '#B696B6', '#80C1C4')
method_colors <- setNames(method_colors, unique(combined_data$Method))

# Create the combined plot
combined_plot <- ggplot(combined_data, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) + 
  geom_text(aes(label = round(Value, 2)), position = position_stack(vjust = 1.1), color = "black",size = 8) +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey80", color = NA),
    text = element_text(color = "black",size = 16, face = "bold"),
    strip.text = element_text(size = 30, face = "bold"), # Bold title box text
    axis.text.x = element_text(size = 20, face = "bold"), # Slanted x-axis text
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.spacing = unit(3, "lines"), # Adjust spacing between panels
    panel.border = element_rect(linetype = "solid", color = "black",fill= NA,size = 1),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    legend.position = "none")
ggsave("/Users/fuyinghao/Documents/STsisal/Revision/figure/combined_ari_purity_plot.png", combined_plot, width = 20, height = 6)





  