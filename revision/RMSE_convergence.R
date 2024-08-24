rm(list = ls())
source('/Users/fuyinghao/Documents/STsisal/function/clean_count.R')
source('/Users/fuyinghao/Documents/STsisal/function/getCellNumber_st.R')
source('/Users/fuyinghao/Documents/STsisal/function/csDeconv_st.R')
source('/Users/fuyinghao/Documents/STsisal/function/mysisal.R')
source('/Users/fuyinghao/Documents/STsisal/function/simplenormalize.R')
source('/Users/fuyinghao/Documents/STsisal/function/sisal.R')
source('/Users/fuyinghao/Documents/STsisal/function/vca.R')
source('/Users/fuyinghao/Documents/STsisal/function/findRefinx.R')
library(deconf)
library(TOAST)
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/MOB/Rep12_MOB_count_matrix-1.RData')
load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/MOB/MOB.dge.sceset.RData')
library(SingleCellExperiment)
library(pbmcapply)
ct.varname = "cellType"
sample.varname = "sampleInfo"
sc_count = assays(sce)$counts
sc_meta = colData(sce)
st_count = MOB_raw
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


library(deconf)
library(TOAST)

Y_raw = as.matrix(st_count)
# Y_raw = clean_count(Y_raw)
K = 5
n_marker = 500
TotalIter = 20
if (is.null(K)) {
  Kres = getCellNumber_st(Y_raw, possibleCellNumber)
  K = Kres$bestK
}
library(CARD)
CARD_obj = createCARDObject(
        sc_count = sc_count,
        sc_meta = sc_meta,
        spatial_count = st_count,
        spatial_location = spatial_location,
        ct.varname = "cellType",
        ct.select = unique(sc_meta$cellType),
        sample.varname = 'sampleInfo',
        minCountGene = 90,
        minCountSpot = 5)
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
n.cell.types = 5
I = st_count
S.start <- CARD_obj@algorithm_matrix$B
Y_raw = Y_raw[match(rownames(S.start),rownames(Y_raw)),]
source('/Users/fuyinghao/Documents/STsisal/Revision/function/csDeconv_st_revised.R')
rmse_matrix = matrix(NA,nrow = 20,ncol = 21)
for(i in 1:20){
  OutRF <- csDeconv_st_revised(Y_raw,S.start , K = K, TotalIter = TotalIter, bound_negative = TRUE,nMarker = n_marker)
  rmse_matrix[i,] = OutRF$allRMSE
}
rmse_matrix = reshape2::melt(rmse_matrix)
rmse_matrix = cbind(rmse_matrix,n_marker = rep(1000,nrow(rmse_matrix)))
rmse_matrix = cbind(rmse_matrix,group = rep('Reference-based',nrow(rmse_matrix)))
# save(rmse_data_long,file = '/Users/fuyinghao/Documents/STsisal/Revision/res/RMSE_convergence.RData')
load('/Users/fuyinghao/Documents/STsisal/Revision/res/RMSE_convergence.RData')

# Your existing ggplot code
new_data = rmse_data_long[rmse_data_long$n_marker==1000,]
new_data = cbind(new_data,group = rep('Randomly initialized',nrow(new_data)))

new_data = rbind(new_data,rmse_matrix)
ggplot(new_data, aes(x = Var2, y = value, color = group)) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "line", size = 1) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "point", size = 3) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, position = position_dodge(0.2)) +
  # scale_color_manual(values = nature_colors, name = "Number of Markers") +  # Set legend name
  labs(x = "Iteration",
       y = "RMSE") +
  theme_minimal() +
  theme(legend.position = "bottom",  # Change the legend position if needed
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))+
  theme(text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
ggsave('/Users/fuyinghao/Documents/STsisal/Revision/figure/RMSE_convergence.pdf', width = 8, height = 8)





