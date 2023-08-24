library(Matrix)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(CARD)
library(spacexr)
library(SPOTlight)
library(Giotto)
library(STdeconvolve)
library(SPOTlight)
library(SeuratObject)
##*****************************************************
## sc_count: scRNA-seq count data must be in the format of matrix or sparseMatrix, while each row represents a gene and each column represents a cell.
## sc_meta: scRNAseq meta data must be in the format of data frame while each row represents a cell. sc_meta data must contain the column indicating the cell type assignment for each cell (e.g., “cellType” column in the example sc_meta data). 
#           Sample/subject information should be provided, if there is only one sample, we can add a column by sc_meta$sampleInfo = "sample1"
## spatial_count：spatial transcriptomics count data, along with spatial location information.
## spatial_location: spatial location information
##*******************************************************

CARD_pipeline <- function(sc_count,sc_meta,spatial_count,spatial_location){
  CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = 'sampleInfo',
    minCountGene = 90,
    minCountSpot = 5) 
  
  ## Deconvolution using CARD 
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  estProp = CARD_obj@Proportion_CARD
  return(estProp = estProp)
}


RCTD_pipeline <- function(sc_count,sc_meta,spatial_count,spatial_location){
  spatial_count <- as(spatial_count,"dgCMatrix")
  puck <- SpatialRNA(spatial_location, spatial_count, colSums(spatial_count),require_int = F)
  sc_count = as(sc_count,"dgCMatrix")
  # nUMI = sc_meta$nUMI
  cell_types = as.factor(sc_meta$cellType)
  names(cell_types) = rownames(sc_meta)
  reference <- Reference(sc_count, cell_types,require_int = FALSE)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8, test_mode = FALSE,CELL_MIN_INSTANCE = 3,UMI_min = 90) # here puck is the SpatialRNA object, and reference is the Reference object.
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  weights = sweep(results$weights, 1, rowSums(results$weights), '/')
  estProp = matrix(weights@x,ncol = ncol(weights))
  colnames(estProp) = colnames(weights)
  rownames(estProp) = rownames(weights)
  return(estProp)
}

STdeconvolve_pipeline <- function(spatial_count,K){
  
  ## remove pixels with too few genes
  counts <- cleanCounts(spatial_count, min.lib.size = 90)
  ## feature select for genes
  corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
  ## choose optimal number of cell-types
  ldas <- fitLDA(t(as.matrix(corpus)), Ks = K)
  ## get best model results
  optLDA <- optimalModel(models = ldas, opt = "min")
  ## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
  results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
  deconProp <- results$theta
  deconGexp <- results$beta
  return(estProp = deconProp)
}

##download code
###
# install.packages("devtools")
# library(devtools)
# install_github("https://github.com/MarcElosua/SPOTlight")
###
SPOTlight_pipeline <- function(sc_count,sc_meta,st_count,spatial_location){
  pbmc <- CreateSeuratObject(counts = sc_count)
  Seurat::Idents(object = pbmc) <- sc_meta$cellType
  pbmc = Seurat::SCTransform(pbmc, verbose = FALSE)
  cluster_markers_all <- Seurat::FindAllMarkers(object = pbmc, 
                                                assay = "SCT",
                                                slot = "data",
                                                verbose = TRUE,
                                                only.pos = TRUE,
                                                logfc.threshold = 0.1)
  pbmc$subclass = sc_meta$cellType
  st_data <- CreateSeuratObject(counts = as.matrix(st_count))
  st_data = Seurat::SCTransform(st_data, verbose = FALSE)
  spotlight_ls_pbmc <- SPOTlight(
    x = pbmc,
    y = st_data@assays$RNA@counts,
    groups = sc_meta$cellType,
    mgs = cluster_markers_all,
    hvg = NULL,
    weight_id = "avg_log2FC",
    group_id = "cluster",
    gene_id = "gene")
  estProp = spotlight_ls_pbmc[["mat"]]
  return(estProp)
}




# ###
# library(remotes)  # if not installed:install.packages('remotes')
# remotes::install_github("RubD/Giotto") 
# ###
SpatialDWLS <- function(sc_count,sc_meta,spatial_count,spatial_location){
  
  ## st data
  spe = matrix(spatial_count,ncol = ncol(spatial_count))
  rownames(spe) = rownames(spatial_count)
  colnames(spe) = colnames(spatial_count)
  instrs = createGiottoInstructions(python_path = "/Users/fuyinghao/miniforge3/bin/python")
  st_data <- createGiottoObject(raw_exprs = spe,instructions = instrs)
  st_data <- normalizeGiotto(gobject = st_data)
  st_data <- calculateHVG(gobject = st_data)
  gene_metadata = fDataDT(st_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
  st_data <- runUMAP(st_data, dimensions_to_use = 1:10)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000) ##******** error
  
  ## sc data
  sce = matrix(sc_count,ncol = ncol(sc_count))
  rownames(sce) = rownames(sc_count)
  colnames(sce) = colnames(sc_count)
  sc_data <- createGiottoObject(raw_exprs = sc_count,instructions = instrs)
  sc_data <- normalizeGiotto(gobject = sc_data)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata = fDataDT(sc_data)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  sc_data <- runPCA(gobject = sc_data, genes_to_use = featgenes, scale_unit = F)
  signPCA(sc_data, genes_to_use = featgenes, scale_unit = F)
  sc_data@cell_metadata$leiden_clus <- as.character(sc_meta[,"cellType"])
  scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                     method = 'scran',
                                                     expression_values = 'normalized',
                                                     cluster_column = 'leiden_clus')
  Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])
  norm_exp<-2^(sc_data@norm_expr)-1
  id<-sc_data@cell_metadata$leiden_clus
  ExprSubset<-norm_exp[Sig_scran,]
  Sig_exp<-NULL
  for (i in unique(id)){
    Sig_exp<-cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp)<-unique(id)
  st_data <- runDWLSDeconv(st_data,sign_matrix = Sig_exp, n_cell = 20) ## wrong
  return(st_data@spatial_enrichment$DWLS)
}


### get proportion aligned using correlation coefficient
GetPropAligned <- function(input,reference,L=K){
  colnames(input)=colnames(reference)=seq(1,dim(input)[2],by=1)
  corMat = cor(input,reference,use="pairwise.complete.obs")
  prop_cor = rep(0,L)
  tmpmat=corMat
  tmpmat[is.na(tmpmat)] = rnorm(L,-1,0.01)
  if(L>2){
    for(i in 1:L){
      maxind = which(tmpmat == max(tmpmat), arr.ind = TRUE)
      prop_cor[maxind[1]] = colnames(corMat)[maxind[2]]
      tmpmat[maxind[1],]=tmpmat[,maxind[2]]=rep(-1,L)
    }
  }else if(L==2){
    if(tmpmat[1,1]>0){
      prop_cor = c("1","2")
    }else{
      prop_cor = c("2","1")
    }
  }
  colnames(input) = prop_cor
  trans_input = input[,colnames(reference)]
  return(trans_input)
}

