rm(list = ls())
sc = anndata::read_h5ad('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/stereoseq_zebrafish/zf24_scRNA.h5ad')
sc_count = t(as.matrix(sc$X))
sc_meta = sc$obs
location = read.csv('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/stereoseq_zebrafish/ST_5_loc.csv')
st_count = read.csv('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P1/stereoseq_zebrafish/ST_5.csv')  

source('/Users/fuyinghao/Documents/STsisal/pipeline/STsisal_pipeline.R')
source('/Users/fuyinghao/Documents/STsisal/pipeline/Deconvolution_pipeline.R')

rownames(st_count) = st_count$X
st_count = st_count[,-1]
rownames(location) = location$X.1
location = location[,-1]
sc_meta = data.frame(cellID = rownames(sc_meta),cellType = sc_meta$celltype_new,sampleInfo = rep('sampleInfo',length(sc_meta)))
rownames(sc_meta) = sc_meta$cellID
length(unique(sc_meta$celltype_new))
colnames(location) = c('x','y')
##STsisal deconvolve
STsisal_res = STsisal_pipeline(t(as.matrix(st_count)),23,n_marker = 1200)
##STdeconvolve
STdeconvolve_impro_res = STdeconvolve_pipeline(t(as.matrix(st_count)),23)
##CARD deconvolve
CARD_res = CARD_pipeline(sc_count,sc_meta,t(as.matrix(st_count)),location)
##RCTD deconvolve
RCTD_res = RCTD_pipeline(sc_count,sc_meta,t(as.matrix(st_count)),location)
##Spotlight
SPOTlight_impro_res = SPOTlight_pipeline(sc_data_B,sc_meta_B,st_count,spatial_location)
save(STsisal_res,STdeconvolve_impro_res,CARD_impro_res,RCTD_impro_res,SPOTlight_impro_res,file = '/Users/fuyinghao/Documents/STsisal/sim_data/final data/P4/B_res.RData')


source('/Users/fuyinghao/Documents/STsisal/function/spatial_pieplot.R')
p1 = spatial_pieplot(STsisal_res$estProp,location)
p3 = spatial_pieplot(RCTD_res,location)
p4 = spatial_pieplot(CARD_res,location)
cowplot::plot_grid(p1,p3,p4,ncol = 3,nrow = 1,labels = c('STsisal','RCTD','CARD'), align = "hv")
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/stereo_mob_comparison.pdf',width = 10*3,height = 10*1,dpi = 3000)


p1 = vizAllTopics(theta = STsisal_res$estProp,
                  pos = location,
                  r = 60,#45
                  lwd = 0,
                  showLegend = TRUE,
                  plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
  
  ## outer border
  ggplot2::geom_rect(data = location,
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +
  
  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +
  
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")


p2 = vizAllTopics(theta = CARD_res,
                  pos = location,
                  r = 60,
                  lwd = 0,
                  showLegend = TRUE,
                  plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
  
  ## outer border
  ggplot2::geom_rect(data = location,
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +
  
  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +
  
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")
# cowplot::plot_grid(p1,p2,ncol = 2,nrow = 1,labels = c('STsisal', 'Stdeconvolve'), align = "hv")
p3 = vizAllTopics(theta = RCTD_res,
                  pos = location,
                  r = 60,
                  lwd = 0,
                  showLegend = TRUE,
                  plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
  
  ## outer border
  ggplot2::geom_rect(data = location,
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +
  
  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +
  
  ## remove the pixel "groups", which is the color aesthetic for the pixel borders
  ggplot2::guides(colour = "none")

cowplot::plot_grid(p1,p2,p3,ncol = 3,nrow = 1,labels = c('STsisal','RCTD','CARD'), align = "hv")
ggsave('/Users/fuyinghao/Documents/STsisal/exp_res/stereo_mob_comparison.pdf',width = 20*3,height = 20*1,dpi = 3000)

dat = read.table('/Users/fuyinghao/Downloads/sc_zf10_meta_data.txt')
head(dat)
