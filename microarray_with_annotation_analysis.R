library(GEOquery)
library(limma)
library(org.Hs.eg.db)

# for collapseBy
source("functions.R")
es <- getGEO("GSE36287", AnnotGPL = TRUE)[[1]]
str(experimentData(es))
str(pData(es))
head(fData(es))
es$`treatment:ch1`
# we should make all our conditions like 'X_X', no spaces and other symbols
es$condition <- gsub("\\-", "_", es$`treatment:ch1`)
es$condition

es <- collapseBy(es, fData(es)$`Gene ID`, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]
#---------------------------------------------------------
fData(es) <- data.frame(row.names = rownames(es))
# We should check our annotation package, org.Mm.eg.db for Mus Muscus, org.Hs.eg.db for Homo sapiense
fData(es)$symbol <- mapIds(org.Hs.eg.db, keys = rownames(es), column = "SYMBOL", keytype = "ENTREZID")

es.qnorm <- es

summary(exprs(es.qnorm))
exprs(es.qnorm) <- normalizeBetweenArrays(log2(exprs(es.qnorm)+1), method="quantile")
summary(exprs(es.qnorm))

es.qnorm.top12K <- es.qnorm
es.qnorm.top12K <- es.qnorm.top12K[head(order(apply(exprs(es.qnorm.top12K), 1, mean), 
                                              decreasing = TRUE), 12000), ]

es.qnorm.top12K <- es.qnorm.top12K[1:12000,]
write.gct(es.qnorm.top12K, file="./es.qnorm.top12k.gct")


p <- pcaPlot(es.qnorm.top12K, 1, 2) + 
  aes(color=condition) 
print(p)

es.design <- model.matrix(~0+condition, data=pData(es.qnorm.top12K))

fit <- lmFit(es.qnorm.top12K, es.design)
# We shoul choose out 2 conditions to compare, it is betere to place treatment before control
fit2 <- contrasts.fit(fit, makeContrasts2(c("condition", "IFNg_treated", "untreated control"),levels=es.design))
fit2 <- eBayes(fit2)
de <- topTable(fit2, adjust.method="BH", number=Inf)
head(de)
library(data.table)
de <- as.data.table(de, keep.rownames=TRUE)
de[symbol == "WARS"]

#----------------------------------------------------------
