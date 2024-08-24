rm(list = ls())
load('/Users/fuyinghao/Documents/STsisal/exp_res/STsisal_breast.RData')
library("org.Hs.eg.db")
library(GOstats)
library(Category)
library(qvalue)
genes = STsisal_res$selMarker$`5`
info <- AnnotationDbi::select(org.Hs.eg.db, keys=genes,
                              columns=c("ENSEMBL"),
                              keytype="SYMBOL")
ens <- info$ENSEMBL

load('/Users/fuyinghao/Documents/STsisal/sim_data/final data/P2/breast cancer_Endecon/breast.st.RData')
all_genes = rownames(breast.st)
names = AnnotationDbi::select(org.Hs.eg.db, keys=all_genes,
                              columns=c("ENSEMBL"),
                              keytype="SYMBOL")
df<-data.frame(name = names$ENSEMBL)

analyze<- function(ens,df){
  sel.entrez<-  unique(na.omit(ens)) # selected gene names (15 character)
  all.entrez<-  unique(na.omit(df$name)) # all gene names in the data matrix under study
  
  # if there is more than one matched, then all are taken to the next step
  
  sel.entrez<-unlist(mget(substr(sel.entrez,1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))
  sel.entrez<-unique(sel.entrez[!is.na(sel.entrez)])
  all.entrez<-unlist(mget(substr(all.entrez,1,15), org.Hs.egENSEMBL2EG,ifnotfound=NA))
  all.entrez<-unique(all.entrez[!is.na(all.entrez)])
  
  # this is hypergeometric testing, by GOstats
  
  params <- new("GOHyperGParams", geneIds=sel.entrez[!is.na(sel.entrez)], universeGeneIds=all.entrez[!is.na(all.entrez)], ontology="BP", pvalueCutoff=0.05,conditional=T, testDirection="over", annotation="org.Hs.eg.db")
  
  over = hyperGTest(params)
  ov<-summary(over)
  ov<-ov[ov[,6]<=500 & ov[,6]>=12,]
  
  for(i in 2:4) ov[,i]<-signif(ov[,i], 3)
  
  # this is the fill in which genes caused the pathway to be significant
  ov<-cbind(ov, idx = rep('',nrow(ov)), name = rep('', nrow(ov)))
  ov[,8]<-as.vector(ov[,8])
  ov[,9]<-as.vector(ov[,9])
  
  for(i in 1:nrow(ov))
  {
    this.genes<-mget(ov[i,1],org.Hs.egGO2ALLEGS)[[1]]
    this.genes<-unique(this.genes[this.genes %in% sel.entrez])
    that<-this.genes[1]
    if(length(this.genes)>1)
    {
      for(j in 2:length(this.genes)) that<-paste(that, ", ", this.genes[j], sep="")
    }
    ov[i,8]<-that
    
    this.genes<-unlist(mget(this.genes, org.Hs.egSYMBOL))
    that<-this.genes[1]
    if(length(this.genes)>1)
    {
      for(j in 2:length(this.genes)) that<-paste(that, ", ", this.genes[j], sep="")
    }
    ov[i,9]<-that
    
  }
  ov <- cbind(ov, p_adjusted = rep(0,nrow(ov)))
  ov$p_adjusted[which(is.na(ov$Pvalue)==FALSE)] <- p.adjust(ov$Pvalue,"BH")
  return(ov)
  
}
res <- analyze(ens,df)
ov = res
write.csv(ov,file = '/Users/fuyinghao/Documents/STsisal/exp_res/EA.csv')
