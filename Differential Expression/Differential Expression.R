#############################

##Author Name: Conor Cremin
##DE Analysis
##Dataset:COVOD-19

#############################

#Packages
if( !require("gplots")){
  BiocManager::install("gplots")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if( !require("DESeq2")){
  BiocManager::install("DESeq2")
}
if( !require("RColorBrewer")){
  BiocManager::install("RColorBrewer")
}
if( !require("limma")){
  BiocManager::install("limma")
}
if( !require("pheatmap")){
  BiocManager::install("pheatmap")
}
if( !require("gridExtra")){
  BiocManager::install("gridExtra")
}
if( !require("grid")){
  BiocManager::install("grid")
}
if( !require("ggplotify")){
  BiocManager::install("ggplotify")
}
if( !require("EnhancedVolcano")){
  BiocManager::install("EnhancedVolcano")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("circlize")){
  BiocManager::install("circlize")
}
if( !require("viridis")){
  BiocManager::install("viridis")
}
if( !require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}

#Load files for Analysis:

countdata = read.csv("raw.counts.csv", check.names = F, row.names = 1)
metadata = read.csv("metadata.csv")
countdata <- countdata[,as.character(metadata$SampleID)] 
keep <- rowSums(countdata) > 10
countdata <- countdata[keep,]
hamster_human = read.csv("hamester_to_human_orthologues.csv", stringsAsFactors = F)
hamster_human = hamster_human[!duplicated(hamster_human$Gene.stable.ID),]
hamster_human = hamster_human[!duplicated(hamster_human$Human.gene.stable.ID),]
hamster_human = hamster_human[hamster_human$Human.gene.stable.ID != "",]
covid_names = read.delim("gene.names.txt", header = F, stringsAsFactors = F)

#Quality Countrol Check (Matching sample names):

all(metadata$SampleID == colnames(countdata))

#DE- Using DESeq2

dds <- DESeqDataSetFromMatrix(countData = countdata, 
                              colData =metadata, design = ~ Condition)
dds <-DESeq(dds)

vsd <- vst(dds, blind = TRUE) 

#sPCA Analysis of DE conditions using log2FC differentially expressed genes:

g = unique(metadata$Condition)[c(3,1,2,4,5)]
myplots = list()
genesets.up = matrix(0,nrow = length(rownames(countdata)), ncol = length(g))
genesets.up = data.frame(genesets.up)
rownames(genesets.up) = rownames(countdata)
for (i in 1:length(g)) {
  group1 = g[i]
  group2 = "Control"
  res <- results(dds, contrast = c("Condition", group1, group2), pAdjustMethod = "BH")
  res <- na.omit(res, cols = c("log2FoldChange", "padj"))
  res <- res[res$log2FoldChange >= 1 & res$padj <= 0.05,]
  m =match(rownames(res), rownames(genesets.up))
  f.a = !is.na(m)
  f.t =m[f.a]
  genesets.up[f.t,i] = res[f.a,2]
  colnames(genesets.up)[i] = as.character(g[i])
}

genes.pca <- prcomp(t(genesets.up),center = T, scale. = F)
df = as.data.frame(genes.pca$x)
col = c('red3','green3','DarkRed','LightSalmon3','LightSlateBlue')
names(col) = rownames(df)

ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3,aes(colour = names(col)) ) +
  labs(x = paste("PC1", paste0("(",round(as.numeric(genes.pca$sdev[1])*1, digits = 1), "%)")),
       y = paste("PC2", paste0("(",round(as.numeric(genes.pca$sdev[2])*1, digits = 1), "%)"))
  ) +
  scale_fill_manual(values=col,
                    breaks=names(col),
                    labels=names(col)) +
  theme_bw()+
  geom_vline(xintercept =0, size = 1) +
  geom_hline(yintercept = 0, size =1) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title=element_text(size=14,face="bold"),
        axis.line = element_line(size =.5),
        axis.ticks = element_line(size = 1),
        axis.text.x = element_text(color = "black", face = "bold", size = 15), 
        axis.text.y = element_text(color = "black", face = "bold", size = 15), 
        legend.title = element_text(color = "black",face = "bold", size = 18,),
        legend.text = element_text(color = "black", face = "bold", size = 12),
        legend.key.size = unit(.5, "cm"),
        legend.key.width = unit(.5,"cm") ,
        legend.position = "bottom")
dev.off()

#Volcano Plots:

group1 = "" #Specify
group2 = "" #Specify
res <- results(dds, contrast = c("Condition", group1, group2), pAdjustMethod = "BH")
res <- na.omit(res, cols = c("log2FoldChange", "padj"))
id = paste(group1, "V", group2)
EnhancedVolcano(as.data.frame(res),
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           lab = rownames(res),
                           selectLab = covid_names$V2,
                           labSize = 3,
                           pointSize = 1,
                           drawConnectors = T,
                           endsConnectors = "last",
                           subtitle = NULL,
                           pCutoff = 0.05,
                           FCcutoff = 1,
                           cutoffLineType = 2,
                           title = paste(id),
                           captionLabSize = 12,
                           col=c('grey30', 'green3', 'royalblue', 'red3'),
                           colAlpha = 1,
                           border = 'partial',
                           #borderWidth = 1,
                           legendPosition = NULL,
                           legendLabSize = 0,
                           legendLabels = c('NS', expression(Log[2]~FC),
                                            'p-value', expression(p-value~and~log[2]~FC))) +
    theme_classic()+
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
#Make data.frame for GSEA Analysis of DEGs from differential Expression. 
#(Genes identified as differentially expressed are given a score of 1 and those that are not are given a 0)

g = unique(metadata$Condition)[c(3,1,2,4,5)]
genesets.up = matrix(0,nrow = length(hamster_human$Golden.Hamster.gene.stable.ID), ncol = length(g))
genesets.up = data.frame(genesets.up)
rownames(genesets.up) = hamster_human$Golden.Hamster.gene.stable.ID
for (i in 1:length(g)) {
  group1 = g[i]         #Change if needed
  group2 = "Control"     #Change if needed
  res <- results(dds, contrast = c("Condition", group1, group2), pAdjustMethod = "BH")
  res <- na.omit(res, cols = c("log2FoldChange", "padj"))
  res <- res[res$log2FoldChange >= 1 & res$padj <= 0.05,]   #Alter this depending on your decided DE Thresholds!!!!!!!!!!!
  m =match(rownames(res), rownames(genesets.up))
  f.a = !is.na(m)
  f.t =m[f.a]
  genesets.up[f.t,i] = 1
  colnames(genesets.up)[i] = as.character(g[i])
}

#genesets.up can be used for GSEA!!!!