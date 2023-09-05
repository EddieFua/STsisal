# STsisal
##STsisal: A Reference-Free Deconvolution Pipeline for Spatial Transcriptomics 
![Pipeline](pipeline.jpg)

We propose a  novel reference-free deconvolution method explicitly designed for ST data. STsisal employs a combination of marker gene selection, mixing ratio decomposition, and cell type characteristic matrix analysis to identify cell types accurately. STsisal is implemented as an open-acess pipeline.

Installation
------------
STsisal is mainly depends on the previous research TOAST and deconf, and deconf has been deprecated. You can download from this repository.

Usage
------------
After download TOAST and deconf, you can source STsisal pipeline directly.
``` r
source(~/pipeline/STsisal_pipeline.R)
source(~/pipeline/STsisal_pipeline.R)
load('~/data/P4/PDAC/Figure4A_layer_annote.RData')
load('~/data/P4/PDAC/PDAC_A_metadata.RData')
load('~/data/P4/PDAC/PDAC_A_sc.RData')
load('~/data/P4/PDAC/PDAC_A_st_count.RData')
load('~/data/P4/PDAC/PDAC_A_st_location.RData')
load('~/data/P4/PDAC/PDAC_B_metadata.RData')
load('~/data/P4/PDAC/PDAC_B_sc.RData')

##data preparation
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

##STsisal deconvolve
STsisal_res = STsisal_pipeline(st_count,20,n_marker = 500,TotalIter = 10)

``` 
