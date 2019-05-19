source("~/Rstudio course/functions.R")
library(data.table)
library(limma)
library(org.Hs.eg.db)
library(Biobase)

# GSE83645 experiment was about psoriasis and uninvolved skin of people with psoriasis

file1 <- read.table('./GSE83645_rawdatacount.txt')
genes <- row.names(file1)
condition <- GEOquery::getGEO("GSE83645", AnnotGPL = TRUE)[[1]]
str(pData(condition))
file1 <- as.matrix(file1)
es <- ExpressionSet(file1)
# 2 types of tissue: psoriasis and uninvolved
pData(es)$tissue <- condition$`tissue type:ch1`
# 5 patiants
pData(es)$patient <- condition$`patient id:ch1`
# 
pData(es)$body_part <- condition$`body part:ch1`

#----------------------PCA plot based on limma
library(limma)
library(ggplot2)
library(ggrepel)
es.qnorm <- es
fData(es)$symbol <- genes
exprs(es.qnorm) <- normalizeBetweenArrays(log2(exprs(es.qnorm) + 1), method="quantile")

es.qnorm.top12K <- es.qnorm
fData(es.qnorm.top12K)$mean <- apply(exprs(es.qnorm.top12K), 1, mean)
es.qnorm.top12K <- es.qnorm.top12K[order(fData(es.qnorm.top12K)$mean, decreasing = TRUE), ]
head(exprs(es.qnorm.top12K))
es.qnorm.top12K <- es.qnorm.top12K[1:12000,]
p <- pcaPlot(es.qnorm.top12K, 1, 2) + 
  aes(color=patient) + 
  geom_text_repel(aes(label=tissue)) 
print(p)
#-------------------------------------------PCA plot based on DESeq + dif.ex. 
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=exprs(es),
                              colData=pData(es),
                              design=~tissue)
dds
dds <- DESeq(dds)
dds
plotDispEsts(dds)
vst <- varianceStabilizingTransformation(dds)
plotPCA(vst, intgroup = 'tissue')


unique(dds$tissue)
de <- results(dds, contrast = c("tissue", "psoriasis", "uninvolved"), cooksCutoff = F)
head(de)
de <- data.table(ID=rownames(de), as.data.table(de))
head(de)


#-----------------------------------------Pathway analysis with msigdbr library 
# fgseaMultilevel from githab devtools::install_github("ctlab/fgsea")
stats <- de[, setNames(stat, entrez)]
library(msigdbr)
# GO BP pathways from MSigDB
m_df <- msigdbr(species = "Homo sapiens")
m_df
pathways <- split(m_df$human_gene_symbol, m_df$gs_name)

library(fgsea)
fr <- fgsea(pathways, stats, nperm = 1000, nproc=4, minSize=15, maxSize=500)
fr[order(pval)]
frML <- fgseaMultilevel(pathways, stats,
                        sampleSize = 100,
                        nproc=4, minSize=15, maxSize=500)
frML[order(pval)]

frML[padj < 0.01]

# Work for very long time, bc I have too many pathways 
collapsedPathways <- collapsePathways(fr[order(pval)][padj < 0.01], pathways, stats)
str(collapsedPathways)


mainPathways <- frML[pathway %in% collapsedPathways$mainPathways][
  order(sign(ES)*log(pval)), pathway]

frMain <- frML[match(mainPathways, pathway)]
frMain[, leadingEdge := lapply(leadingEdge, mapIds, 
                               x=org.Hs.eg.db, keytype="ENTREZID", column="SYMBOL")]
dir.create("gsea")
fwrite(frMain, file="gsea/Control.vs.IFN-y-IL4-12h.tsv", sep="\t", sep2=c("", " ", ""))

pdf("gsea/Control.vs.IFN-y-IL4-12h.pdf", width=12, height=2 + length(mainPathways) * 0.25)
plotGseaTable(pathways = pathways[mainPathways], stats = stats, fgseaRes=frMain, gseaParam = 0.5)
dev.off()

fr[, leadingEdge := NULL]
fwrite(frML[order(sign(ES)*log(pval))], file="./gsea/Control.vs.IFN-y-IL4-12h.full.tsv", sep="\t") 

#---------------------------------Enrichment
library(ggplot2)
plotEnrichment(pathways[["SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6"]], stats) + 
  ggtitle("SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6")

plotEnrichment(pathways[["REACTOME_CELL_CYCLE"]], stats) + 
  ggtitle("REACTOME_CELL_CYCLE")

plotEnrichment(pathways[["GO_INNATE_IMMUNE_RESPONSE"]], stats) + 
  ggtitle("GO_INNATE_IMMUNE_RESPONSE")






