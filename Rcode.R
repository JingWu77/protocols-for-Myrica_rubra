########################################################################################################################
# Documents and scripts were written by: Jing Wu
# For manuscript: Jing Wu, et al. (2026). "Genomic origin, domestication, and gene expression control of fruit traits in Chinese bayberry with bHLH-IBH1 mediated anthocyanin biosynthesis" (under review). 
# email: 12307018@zju.edu.cn
# Jun Chen Lab. 
########################################################################################################################



########################################################################################################################
# Part 1 GWAS analysis
########################################################################################################################
## for fruit color
library(GAPIT)
library(CMplot)
geno_raw <- read.table("my_geno_data.raw", header = TRUE)
myGD <- geno_raw[, -c(1, 3:6)]
colnames(myGD)[1] <- "Taxa"
bim <- read.table("myrica_135.snp.allfiltered.LD.bim", header = FALSE)
myGM <- bim[, c(2, 1, 4)]
colnames(myGM) <- c("SNP", "Chromosome", "Position")
pheno_all <- read.csv("135_pheno_new.csv", header = FALSE)
ids <- read.table("135_id.txt", header = FALSE)
myY <- data.frame(Taxa = ids[,1], Color = pheno_all[,1])
if(any(is.na(myGD_data))) {
  fill_mean <- function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  }
  myGD_data <- apply(myGD_data, 2, fill_mean)
}
col_vars <- apply(myGD_data, 2, var)
keep_logic <- col_vars > 0
myGD_data <- myGD_data[, keep_logic]
myGM <- myGM[keep_logic, ]
myGD <- data.frame(Taxa = myGD_taxa, myGD_data, stringsAsFactors = FALSE)
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM,
  PCA.total = 3,          
  model = c("GLM", "MLM", "MLMM", "FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  file.output = TRUE,
)
## for fruit size
library(GAPIT)
library(CMplot)
geno_raw <- read.table("my_geno_64.raw", header = TRUE)
myGD <- geno_raw[, -c(1, 3:6)] 
colnames(myGD)[1] <- "Taxa"
bim <- read.table("myrica_64_cleaned.bim", header = FALSE)
myGM <- bim[, c(2, 1, 4)]
colnames(myGM) <- c("SNP", "Chromosome", "Position")
pheno_data <- read.table("64_c.txt", header = FALSE)
ids <- read.table("64_id.txt", header = FALSE)
myY <- data.frame(Taxa = as.character(ids[,1]), 
                  Size = as.numeric(pheno_data[,1]))
myGD_taxa <- as.character(myGD[, 1])
myGD_numeric <- as.matrix(myGD[, -1])
fill_mean <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}
myGD_numeric <- apply(myGD_numeric, 2, fill_mean)
vars <- apply(myGD_numeric, 2, var)
keep_idx <- which(vars > 0)
myGD <- data.frame(Taxa = myGD_taxa, myGD_numeric[, keep_idx], stringsAsFactors = FALSE)
myGM <- myGM[keep_idx, ]
myGAPIT <- GAPIT(
  Y = myY,
  GD = myGD,
  GM = myGM,
  PCA.total = 3,           
  model = c("GLM", "MLM", "MLMM", "FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  file.output = TRUE
)

########################################################################################################################
# Part 2 Positive selection and visualization
########################################################################################################################
## xpclr-select top 5% gene
library("GenomicRanges")
library(plotrix)
library("data.table")
setwd("/data/wujing/xpclr/")
out <- read.table("./xpclr_top5.txt", header = T)
length(out$chrom)
gene <- read.table("./ref.gene.bed", header = T)
gene_GR <- GRanges(seqnames = gene$chrom, IRanges(start = gene$start, end = gene$end), strand = gene$strand, id = gene$id)
out_GR <- GRanges(seqnames = out$chrom, IRanges(start = out$start, end = out$end), snp = out$SNP, xpclr = out$XPCLR)
out_GR2 <- out_GR 
cw_union <- GenomicRanges::union (out_GR,out_GR2)
write.table(cw_union, file="cw.xpclr.non.overlapping.top1.txt", quote = FALSE, row.names = FALSE) 
cw_gene_op <- findOverlaps(cw_union, gene_GR) 
cw_op_gene <- cw_union[cw_gene_op@from,] 
cw_op_gene 
gene_op_cw<- gene_GR[cw_gene_op@to,] 
gene_op_cw 
write.table(gene_op_cw, file="gene_op_xpclr_cw.txt", quote = FALSE, row.names = FALSE)
## xpehh-slide window
library("devtools")
library(dplyr)
library(tidyr)
library(ggplot2)
library(vroom)
library("data.table")
cw_xpehh <- fread('./cw.final.xpehh.txt', sep = "\t", header = TRUE)
pos_win <- winScan(x = cw_xpehh,groups = "chr",position = "pos",values = c("normxpehh"), win_size = 20000,win_step = 2000,funs = c("mean"),cores = 40)
write.table(pos_win, file = "./cw.xpehh.norm.mean.txt", quote = FALSE, row.names = FALSE)
## xpclr-select top 5% gene
library("GenomicRanges")
library(plotrix)
library("data.table")
cw_xpehh <- read.table("./cw.xpehh.norm.mean.txt", header = T)
length(cw_xpehh$chr)
gene <- read.table("./ref.gene.bed", header = T)
gene_GR <- GRanges(seqnames = gene$chrom, IRanges(start = gene$start, end = gene$end), strand = gene$strand, id = gene$id)
sort_cw <- cw_xpehh[order(cw_xpehh$normxpehh_mean, decreasing = T),]
out <- sort_cw[1:7259,]
summary(out$normxpehh_mean)
out_GR <- GRanges(seqnames = out$chr, IRanges(start = out$win_start, end = out$win_end), normxpehh_n = out$normxpehh_n, normxpehh_mean = out$normxpehh_mean)
out_GR2 <- out_GR 
cw_union <- GenomicRanges::union (out_GR,out_GR2)
write.table(cw_union, file="cw.xpehh.non.overlapping.txt", quote = FALSE, row.names = FALSE) 
###
cw_gene_op <- findOverlaps(cw_union, gene_GR) 
cw_op_gene <- cw_union[cw_gene_op@from,] 
cw_op_gene 
gene_op_cw<- gene_GR[cw_gene_op@to,] 
gene_op_cw 
write.table(gene_op_cw, file="gene_op_xpehh_cw.txt", quote = FALSE, row.names = FALSE)
## visualize by cmplot
xpclr <- read.table("myrica_xpclr_final.txt", fill = TRUE, header = T)
xpclr <- na.omit(xpclr)
xpclr_quantile <- quantile(xpclr2$XPCLR, 0.95)
xpclr <- xpclr[,c(1,2,5,6,13)]
colnames(xpclr) <- c("SNP","chrom","start","end","XPCLR")
SNPs1 <- c("2_00138001_00158000",
           "2_33806001_33826000",
           "1_17408001_17428000",
           "1_17720001_17740000",
           "1_22032001_22052000",
           "3_31950001_31970000",
           "7_04152001_04172000",
           "7_04156001_04176000",
           "7_25940001_25960000",
           "6_06930001_06950000",
           "6_29544001_29564000",
           "6_31278001_31298000",
           "6_31280001_31300000")        
SNPs2 <- c("2_19488001_19508000",
           "5_07002001_07022000",
           "5_07004001_07024000",
           "3_07730001_07750000",
           "6_30656001_30676000")
SNPs3 <- c("2_26480001_26500000",
           "1_07336001_07356000",
           "1_34650001_34670000",
           "1_41682001_41702000",
           "1_41710001_41730000",
           "3_07136001_07156000",
           "3_14222001_14242000",
           "3_32318001_32338000",
           "6_29012001_29032000",
           "3_38182001_38202000",
           "7_04026001_04046000",
           "7_01352001_01372000",
           "6_31204001_31224000",
           "3_14012001_14032000",
           "8_16748001_16768000")
SNPs4 <- c("5_28360001_28380000",
           "5_28364001_28384000",
           "5_30664001_30684000",
           "2_33780001_33800000",
           "2_33740001_33760000",
           "2_33754001_33774000",
           "3_29468001_29488000",
           "7_04652001_04672000",
           "3_28446001_28466000",
           "3_28454001_28474000",
           "8_05524001_05544000",
           "1_41152001_41172000",
           "1_41172001_41192000",
           "1_41234001_41254000",
           "7_11942001_11962000",
           "7_11948001_11968000")
SNPs5 <- c("6_32712001_32732000",
           "2_20066001_20086000",
           "5_28648001_28668000",
           "4_02908001_02928000",
           "7_10364001_10384000",
           "1_42240001_42260000",
           "3_38184001_38204000",
           "8_08860001_08880000",
           "6_28788001_28808000",
           "4_24462001_24482000",
           "2_19530001_19550000")
SNPs6 <- c("5_28298001_28318000",
           "1_07978001_07998000",
           "6_31206001_31226000")
genes1 <- c("RFR1",
            "BAN",
            "RNS1",
            "RNS1",
            "EIN3",
            "AAE13",
            "KFB",
            "KFB",
            "ARM",
            "MYB3",
            "AHP4",
            "TT7",
            "TT7")
genes2 <- c("CYP98A3",
            "SHT",
            "SHT",
            "CYP82C2",
            "ALDH2C4")
genes3 <- c("WDL4",
            "OFP13",
            "CYP71",
            "PSY1R",
            "PSY1R",
            "DA1",
            "ER",
            "IQD26",
            "IQD26",
            "CYP714A1",
            "ARFA1F",
            "BGL1",
            "GA3OX1",
            "CKX7",
            "SPL12")
genes4 <- c("GLXI-LIKE;9",
            "LOI1",
            "GGPS1",
            "TPS03",
            "TPS14",
            "TPS14",
            "CYP76C4",
            "CYP71B36",
            "TKL1",
            "TKL1",
            "LOX1",
            "LOX2",
            "LOX2",
            "LOX2",
            "UGT85A2",
            "UGT85A2")
genes5 <- c("SPS1F",
            "AGLU1",
            "GULLO6",
            "THI1",
            "THI1",
            "THFS",
            "CYP714A1",
            "GME",
            "MIOX4",
            "UMAMIT41",
            "UMAMIT9")
genes6 <- c("PFL",
            "SE",
            "LUG")
colors=c("#B2B6C1","#F0EFED","#B2B6C1","#F0EFED","#B2B6C1","#F0EFED","#B2B6C1","#F0EFED")
highlight_SNPs <- c(SNPs1, SNPs2,SNPs3,SNPs4,SNPs5,SNPs6)
highlight_genes <- c(genes1, genes2,genes3,genes4,genes5,genes6)
highlight_colors <- c(rep("#44036f", 13), rep("#ff7400", 5),rep("#00782d", 15),rep("#0b61a4", 16),rep("#b52d43", 11),rep("#83a000", 3))
pdf("XPCLR_0.95_s12.pdf",width = 12,height = 12)
CMplot(xpclr2,
       plot.type="m",
       LOG10=FALSE,
       col=colors,
       points.alpha=10,
       threshold=xpclr_quantile,threshold.lty=2,threshold.lwd=1,threshold.col="black",
       amplify=FALSE,band=0.5,cex=0.5,ylim=c(0,20),main="",file="pdf",
       ylab="XPCLR",chr.border=FALSE,dpi=300,file.output=FALSE,
       verbose=TRUE,
       highlight = highlight_SNPs,
       highlight.col = highlight_colors, 
       highlight.cex = 1,
       highlight.pch = 16,
       highlight.text = highlight_genes,
       highlight.text.cex = 1.1,
       highlight.text.font = 4, 
       highlight.text.col = highlight_colors)
dev.off()

########################################################################################################################
# Part 3 Gene differential expression analysis using R package 'edgeR'
########################################################################################################################
## input and filter data
library(edgeR)
x<-read.csv("allfruit_counts.csv",sep=",",row.names="Geneid")
group<-c(rep("A",27),rep("B",27),rep("C",27),rep("D",27),rep("E",27),rep("F",27))
y <- DGEList(counts=x, group=group)
keep <- rowSums(cpm(y)>1) >= 1   
y <- y[keep,,keep.lib.sizes=FALSE]
design <- model.matrix(~ 0 + group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)   
fit <- glmQLFit(y,design,robust=TRUE) 
y$common.dispersion
## differential expression analysis of each two groups
#group A vs B
qlf_AvsB <- glmQLFTest(fit,contrast=c(-1,1,0,0,0,0))
FDR_AvsB <- p.adjust(qlf_AvsB$table$PValue, method="BH")  
Sig_AvsB <- qlf_AvsB$table 
Sig_AvsB$padj=FDR_AvsB   
Sig_AvsB <- qlf_AvsB$table[which(FDR_AvsB<0.05),]  
write.csv(Sig_AvsB,file="Sig_AvsB_0.05_edgeR.csv") 
df <- read.table("Sig_AvsB_0.05_edgeR.csv", sep=",",header = T)
gene_diff <- df[order(df$logFC, decreasing = TRUE),]
gene_diff[which(gene_diff$logFC >= 1),'sig'] <- 'up'
gene_diff[which(gene_diff$logFC <= -1),'sig'] <- 'down'
gene_diff[which(abs(gene_diff$logFC) <= 1),'sig'] <- 'none'
gene_diff_select <- subset(gene_diff, sig %in% c('up', 'down'))
write.csv(gene_diff_select, file = 'AvsB_005_1_DEG.csv',row.names = F, quote = F)
#group A vs C
qlf_AvsC <- glmQLFTest(fit,contrast=c(-1,1,0,0,0,0))
FDR_AvsC <- p.adjust(qlf_AvsC$table$PValue, method="BH")  
Sig_AvsC <- qlf_AvsC$table 
Sig_AvsC$padj=FDR_AvsC   
Sig_AvsC <- qlf_AvsC$table[which(FDR_AvsC<0.05),]  
write.csv(Sig_AvsC,file="Sig_AvsC_0.05_edgeR.csv") 
df <- read.table("Sig_AvsC_0.05_edgeR.csv", sep=",",header = T)
gene_diff <- df[order(df$logFC, decreasing = TRUE),]
gene_diff[which(gene_diff$logFC >= 1),'sig'] <- 'up'
gene_diff[which(gene_diff$logFC <= -1),'sig'] <- 'down'
gene_diff[which(abs(gene_diff$logFC) <= 1),'sig'] <- 'none'
gene_diff_select <- subset(gene_diff, sig %in% c('up', 'down'))
write.csv(gene_diff_select, file = 'AvsC_005_1_DEG.csv',row.names = F, quote = F)
#PS: Run this code repeatedly until all the combinations of the two groups have been compared.
## select the union of all differentially expressed genes dataset
library(tidyverse)
library(pacman)
df <- read.csv("allfruit_DEGs.csv", header = T) %>% as_tibble
str(df)
list <- as.list(df)
x1 <- list$G1G2
x2 <- list$G1G3
x3 <- list$G1G4
x4 <- list$G1G5
x5 <- list$G1G6
x6 <- list$G2G3
x7 <- list$G2G4
x8 <- list$G2G5
x9 <- list$G2G6
x10 <- list$G3G4
x11 <- list$G3G5
x12 <- list$G3G6
x13 <- list$G4G5
x14 <- list$G4G6
x15 <- list$G5G6
a = union(union(x1,x2),x3)
b = union(union(a,x4),x5)
c = union(union(b,x6),x7)
d = union(union(c,x8),x9)
e = union(union(d,x10),x11)
f = union(union(e,x12),x13)
g = union(union(f,x14),x15)
write.csv(g,file = "allfruit_DEGs_union.csv")

########################################################################################################################
# Part 4 WGCNA analysis
########################################################################################################################
##step 1 input gene and traits data
library(WGCNA)
library(FactoMineR)
library(factoextra)
library(tidyverse) 
library(data.table) 
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores())
dir.create(tempdir())
tpm0 <- read.table("G6_DEG_tpm.csv", sep=",",header = T)
gene <- tpm0$Geneid
rownames(tpm0) = tpm0[,1]
tpm0 <- tpm0[,-1]
data <- log2(tpm0+1)
datTraits <- read.table("alltraits.csv", sep=",",header = T)
rownames(datTraits2) <- datTraits[,1]
datTraits2 = datTraits[,-1]
# find missing
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
# pca and sample clustering
group_list <- datTraits$Group
dat.pca <- PCA(datExpr0, graph = F)
pca <- fviz_pca_ind(dat.pca,title = "Principal Component Analysis",addEllipses = FALSE,legend.title = "Groups", geom.ind = "point",pointsize = 5,labelsize = 5,fill.ind = group_list,col.ind = group_list,repel = TRUE, mean.point=F)+theme_void()+coord_fixed(ratio = 1)+scale_shape_manual(values=c(16,16,16,16,16,16))
pca <- pca + theme(panel.border = element_rect(color = "black", size = 1, fill = NA),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),axis.title.y =element_text(size=14,angle = 90))
sampleTree <- hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, cex.axis = 1, cex.main = 1,cex.lab=1)
sample_colors <- numbers2colors(as.numeric(factor(datTraits$Group)), colors = rainbow(length(table(datTraits$Group))), signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
pdf("Sample dendrogram and trait.pdf",width = 8,height = 6)
plotDendroAndColors(sampleTree, sample_colors, groupLabels = "trait", marAll = c(1, 4, 3, 1),main = "Sample dendrogram and trait" )
dev.off()
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")
##step2 select the best power value
R.sq_cutoff = 0.85
powers <- c(seq(1,20,by = 1), seq(22,30,by = 2))
sft <- pickSoftThreshold(datExpr,networkType = "unsigned",powerVector = powers,verbose = 5)
sft$powerEstimate
par(mfrow = c(1,2))
cex1 = 0.9
pdf("sft_power.pdf",width = 16,height = 12)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
dev.off()
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
power=9
save(sft, power, file = "step2_power_value.Rdata")
##step3 construct a weighted co-expression network
net = blockwiseModules(datExpr, power = power,TOMType = "unsigned", maxBlockSize = ncol(datExpr), minModuleSize = 80,reassignThreshold = 0,mergeCutHeight = 0.25,numericLabels = TRUE,saveTOMs = F,saveTOMFileBase = "DKTOM",verbose = 3)
table(net$colors)
moduleColors <- labels2colors(net$colors)
table(moduleColors)
pdf("20m_genes_modules_ClusterDendrogram.pdf",width = 16,height = 12)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
save(net, moduleColors, file = "20m_step3_genes_modules.Rdata")
##step4  association of gene module and phenotype 
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits2, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
names(MEs)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")  
dim(textMatrix) = dim(moduleTraitCor)
pdf("20m_Module-trait-relationship_heatmap.pdf",width = 2*length(colnames(datTraits2)), height = 0.6*length(names(MEs)) )
par(mar=c(10, 9, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = colnames(datTraits2),yLabels = names(MEs),ySymbols = names(MEs),colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1,1),main = "Module-trait relationships")
dev.off()
save(datTraits,datTraits2, file = "20m_step4_design.Rdata")
gene_module <- data.frame(gene=colnames(datExpr),module=moduleColors)
write.csv(gene_module,file = "20m_gene_moduleColors.csv",row.names = F, quote = F)

########################################################################################################################
# Part 5 MFUZZ analysis
########################################################################################################################
## for DK
library("Mfuzz")
gene <- read.table("dk_allfruit_8colormodule_count_meantpm.csv",header = T,row.names=1,sep=",")
gene_counts <- data.matrix(gene)
eset <- new("ExpressionSet",exprs = gene_counts)
gene.r <- filter.NA(eset, thres=0.25)                   
gene.f <- fill.NA(gene.r,mode="mean")               
counts <- filter.std(gene.f,min.std=0)                  
gene.s <- standardise(counts)                           
c <- 6                                                                   
m <- mestimate(gene.s)                                      
cDK <- mfuzz(gene.s, c = c, m = m)
cDK$size                                                                  
cDK$cluster[cDK$cluster == 1]                                  
write.table(cDK$cluster,"output_dk_6_new.txt",quote=F,row.names=T,col.names=F,sep="\t")       
pdf("dk_mfuzz_6_new.pdf",width = 12, height = 6)
mfuzz.plot2(gene.s, cl=cDK,mfrow=c(2,3),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()
## for SJ
library("Mfuzz")
gene <- read.table("sj_allfruit_8colormodule_count_meantpm.csv",header = T,row.names=1,sep=",")
gene_counts <- data.matrix(gene)
eset <- new("ExpressionSet",exprs = gene_counts)
gene.r <- filter.NA(eset, thres=0.25)                   
gene.f <- fill.NA(gene.r,mode="mean")               
counts <- filter.std(gene.f,min.std=0)                  
gene.s <- standardise(counts)                           
c <- 6                                                                   
m <- mestimate(gene.s)                                      
cSJ <- mfuzz(gene.s, c = c, m = m)
cSJ$size                                                                  
cSJ$cluster[cSJ$cluster == 1]                                  
write.table(cSJ$cluster,"output_SJ_6_new.txt",quote=F,row.names=T,col.names=F,sep="\t")       
pdf("SJmfuzz_6_new.pdf",width = 12, height = 6)
mfuzz.plot2(gene.s, cl=cSJ,mfrow=c(2,3),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()
## for BQ
library("Mfuzz")
gene <- read.table("bq_allfruit_8colormodule_count_meantpm.csv",header = T,row.names=1,sep=",")
gene_counts <- data.matrix(gene)
eset <- new("ExpressionSet",exprs = gene_counts)
gene.r <- filter.NA(eset, thres=0.25)                   
gene.f <- fill.NA(gene.r,mode="mean")               
counts <- filter.std(gene.f,min.std=0)                  
gene.s <- standardise(counts)                           
c <- 6                                                                   
m <- mestimate(gene.s)                                      
cBQ <- mfuzz(gene.s, c = c, m = m)
cBQ$size                                                                  
cBQ$cluster[cBQ$cluster == 1]                                  
write.table(cBQ$cluster,"output_BQ_6_new.txt",quote=F,row.names=T,col.names=F,sep="\t")       
pdf("BQmfuzz_6_new.pdf",width = 12, height = 6)
mfuzz.plot2(gene.s, cl=cBQ,mfrow=c(2,3),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()




########################################################################################################################
# Part 6 Spatial transcriptome analysis
########################################################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(glmGamPoi)
library(devtools)
## for WJ16 sample
setwd("D:/CJlab/myrica/spatial/WJ16")
data_dir <-"E:/spatial/WJ16"
file_name <-"filtered_feature_bc_matrix.h5"
myrica <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice ="fruit1")
myrica@project.name <-"fruit1"
Idents(myrica) <-"fruit1"
myrica$orig.ident <-"fruit1"
p1 <- VlnPlot(myrica, features ="nCount_Spatial",pt.size = 0,cols ="tomato") + NoLegend()
p2 <- SpatialFeaturePlot(myrica, features ="nCount_Spatial",pt.size.factor = 4,stroke = 0.0) +theme(legend.position ="right")
p1 | p2
p3 <- VlnPlot(myrica, features ="nFeature_Spatial",pt.size = 0,cols ="tomato") + NoLegend()
p4 <- SpatialFeaturePlot(myrica, features ="nFeature_Spatial",pt.size.factor = 4,stroke = 0.0) +theme(legend.position ="right")
p3 | p4
myrica=myrica[,unname(which(colSums(GetAssayData(myrica))!=0))]
myrica <- SCTransform(myrica, assay = "Spatial", verbose = FALSE)
top10000 <- FindVariableFeatures(myrica, nfeatures = 10000)
top10000 <- ScaleData(top10000)
top10000 <- RunPCA(top10000, features = VariableFeatures(object = top10000), verbose = FALSE, npcs= 100)
top10000 <- FindNeighbors(top10000, reduction = "pca", dims = 1:30)
top10000 <- FindClusters(top10000, resolution = 0.5, verbose = FALSE)
top10000 <- RunUMAP(top10000, reduction = "pca", n.neighbors = 50, min.dist = 0.01, dims = 1:20)
p9 <- DimPlot(top10000, reduction = "umap", label = TRUE)+ scale_color_manual(values = cluster_colors)
p10 <- SpatialDimPlot(top10000, label = FALSE,pt.size.factor = 3,stroke = 0.0)+scale_fill_manual(values = cluster_colors)
top10000 <- RunTSNE(top10000, reduction = "pca", n.neighbors = 50, min_dist = 0.01, dims = 1:20,check_duplicates = FALSE)
p11 <- DimPlot(top10000, reduction = "tsne", label = TRUE)
p12 <- SpatialDimPlot(top10000, label = FALSE,pt.size.factor = 2.8,stroke = 0.0)
markers <- FindAllMarkers(myrica, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25,test.use = "wilcox")
write.csv(markers, "WJ16_cluster_markers_0.1.csv", row.names = FALSE)
## for WJ26 sample
setwd("D:/CJlab/myrica/spatial/WJ26")
data_dir <-"E:/spatial/E20240400-02-01/processFile_RNA/WJ01-26"
file_name <-"filtered_feature_bc_matrix.h5"
myrica <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice ="fruit1")
myrica@project.name <-"fruit2"
Idents(myrica) <-"fruit2"
myrica$orig.ident <-"fruit2"
myrica=myrica[,unname(which(colSums(GetAssayData(myrica))!=0))]
myrica <- SCTransform(myrica, assay = "Spatial", verbose = FALSE)
top10000 <- FindVariableFeatures(myrica, nfeatures = 10000)
top10000 <- ScaleData(top10000)
top10000 <- RunPCA(top10000, features = VariableFeatures(object = top10000), verbose = FALSE, npcs= 100)
top10000 <- FindNeighbors(top10000, reduction = "pca", dims = 1:30)
top10000 <- FindClusters(top10000, resolution = 0.5, verbose = FALSE)
top10000 <- RunUMAP(top10000, reduction = "pca", n.neighbors = 50, min.dist = 0.01, dims = 1:30)
p9 <- DimPlot(top10000, reduction = "umap", label = TRUE)+ scale_color_manual(values = cluster_colors)
p10 <- SpatialDimPlot(top10000, label = FALSE,pt.size.factor = 3,stroke = 0.0)+scale_fill_manual(values = cluster_colors)
top10000 <- RunTSNE(top10000, reduction = "pca", n.neighbors = 50, min_dist = 0.01, dims = 1:30,check_duplicates = FALSE)
p11 <- DimPlot(top10000, reduction = "tsne", label = TRUE)
p12 <- SpatialDimPlot(top10000, label = FALSE,pt.size.factor = 2.8,stroke = 0.0)
markers <- FindAllMarkers(myrica, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25,test.use = "wilcox")
write.csv(markers, "WJ26_cluster_markers_0.1.csv", row.names = FALSE)
## for WJ43 sample
setwd("D:/CJlab/myrica/spatial/WJ43")
data_dir <-"E:/spatial/WJ43"
file_name <-"filtered_feature_bc_matrix.h5"
myrica <- Load10X_Spatial(data.dir = data_dir, filename = file_name,slice ="fruit1")
myrica@project.name <-"fruit3"
Idents(myrica) <-"fruit3"
myrica$orig.ident <-"fruit3"
myrica=myrica[,unname(which(colSums(GetAssayData(myrica))!=0))]
myrica <- SCTransform(myrica, assay = "Spatial", verbose = FALSE)
top10000 <- FindVariableFeatures(myrica, nfeatures = 10000)
top10000 <- ScaleData(top10000)
top10000 <- RunPCA(top10000, features = VariableFeatures(object = top10000), verbose = FALSE, npcs= 100)
top10000 <- FindNeighbors(top10000, reduction = "pca", dims = 1:30)
top10000 <- FindClusters(top10000, resolution = 0.5, verbose = FALSE)
top10000 <- RunUMAP(top10000, reduction = "pca", n.neighbors = 50, min.dist = 0.01, dims = 1:50)
p9 <- DimPlot(top10000, reduction = "umap", label = TRUE)+ scale_color_manual(values = cluster_colors)
p10 <- SpatialDimPlot(top10000, label = FALSE,pt.size.factor = 3,stroke = 0.0)+scale_fill_manual(values = cluster_colors)
top10000 <- RunTSNE(top10000, reduction = "pca", n.neighbors = 50, min_dist = 0.01, dims = 1:50,check_duplicates = FALSE)
p11 <- DimPlot(top10000, reduction = "tsne", label = TRUE)
p12 <- SpatialDimPlot(top10000, label = FALSE,pt.size.factor = 2.8,stroke = 0.0)
markers <- FindAllMarkers(myrica, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25,test.use = "wilcox")
write.csv(markers, "WJ43_cluster_markers_0.1.csv", row.names = FALSE)
