library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library("ggplot2")
library("ChIPpeakAnno")
library(gprofiler2)
library(diffloop)
library(GenomicRanges)
library(Repitools)
library(DOSE)
library("EnsDb.Hsapiens.v86")
library("org.Hs.eg.db")
library("enrichplot")
library(msigdbr)

file1 = "H3K27ac.bed"
file1_name = basename(file1)
  
df.file1 = read.table(file1, header = F, sep="\t", stringsAsFactors = F)
colnames(df.file1)[1] =  "chr"
colnames(df.file1)[2] =  "start"
colnames(df.file1)[3] =  "end"
colnames(df.file1)[4] =  "score"

# plot(density(df.file$enrich))
head(df.file1, 2)
file1_name


# Reading the file with the peaks :
peak2 = toGRanges(df.file1)

# peak = readPeakFile(FILE, as="GRanges")
# head(peak)
# peak3 = makeGRangesFromDataFrame(df.file1, keep.extra.columns=TRUE) 



# Annotations based on TxDb database

peak_ucsc = addchr(peak2)

peakAnno_ucsc <- annotatePeak(peak_ucsc, tssRegion=c(-3000, 3000),
                              TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                              annoDb="org.Hs.eg.db",  
                              addFlankGeneInfo = TRUE,
  # flankDistance = 1000000,
  sameStrand = FALSE,
  ignoreOverlap = FALSE,
  ignoreUpstream = FALSE,
  ignoreDownstream = FALSE,
  # overlap = "all",
  verbose = TRUE)

peakAnno_ucsc_GR <- as.GRanges(peakAnno_ucsc)
peakAnno_ucsc_DF <- as.data.frame(peakAnno_ucsc)

dim(df.file)
length(peakAnno_ucsc_GR)
dim(peakAnno_ucsc_DF)[1]

write.table(peakAnno_ucsc_DF, file=paste(FILE_NAME, "genome.ucsc.annotations.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)  
length(unique(peakAnno_ucsc_DF$SYMBOL))

# Annotation Pie
plotAnnoPie(peakAnno_ucsc)

png(paste(file1_name, "genome.ucsc.annotations.pie.png", sep="."), width = 400, height = 400, units = "px")
plotAnnoPie(peakAnno_ucsc)
dev.off()

png(paste(file1_name, "genome.ucsc.annotations.upset.png", sep="."), width = 400, height = 400, units = "px")
upsetplot(peakAnno_ucsc)
dev.off()

upsetplot(peakAnno_ucsc)



# Functional enrichment analysis
# Selecting the genes that we will work with

# Functional enrichment analysis
# Selecting the genes that we will work with

length(as.data.frame(peakAnno_ucsc)$geneId)
GENES = unique(as.data.frame(peakAnno_ucsc)$geneId)
length(GENES)

result = try({
  
egoo <- enrichGO(gene    = GENES,
                OrgDb     = org.Hs.eg.db,
                # keyType  = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1, readable = TRUE)

dotplot(egoo)
ggsave(paste(FILE_NAME, "genome.ucsc.enrich.GO.BP.png", sep="."), 
               width = 2800,
               height = 2800,
               units = "px",
               dpi = 300 )
    
write.table(as.data.frame(egoo), file=paste(FILE_NAME, "genome.ucsc.enrich.GO.BP.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

}, silent = TRUE)

result = try({
  
egoo2 <- enrichGO(gene    = GENES,
                OrgDb     = org.Hs.eg.db,
                # keyType  = 'ENSEMBL',
                ont           = "MF",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1, readable = TRUE)
    

dotplot(egoo2)
    
ggsave(paste(FILE_NAME, "genome.ucsc.enrich.GO.MF.png", sep="."), 
   width = 2800,
   height = 2800,
   units = "px",
   dpi = 300 ) 
    
write.table(as.data.frame(egoo2), file=paste(FILE_NAME, "genome.ucsc.enrich.GO.MF.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

}, silent = TRUE)



# KEGG pathways
# over-representation analysis 
# gene set enrichment analysis 
# for pathway visualization : library("pathview")

result = try({
hsa <- search_kegg_organism('Homo sapiens', by='scientific_name')
})

result = try({
  
mkk1 <- enrichKEGG(gene = GENES,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)
  

dotplot(mkk1)

ggsave( paste(FILE_NAME, "genome.ucsc.kegg1.png", sep="."), 
   width = 2800,
   height = 2800,
   units = "px",
   dpi = 300 ) 
    
write.table(as.data.frame(mkk1), 
            file = paste(FILE_NAME, "genome.ucsc.enrich.KEGG1.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    
}, silent = TRUE)



result = try({

mkk2 <- gseKEGG(gene = GENES,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)

dotplot(mkk2)

ggsave( paste(FILE_NAME, "genome.ucsc.kegg2.png", sep="."), 
        width = 2800,
        height = 2800,
        units = "px",
        dpi = 300 ) 
write.table(as.data.frame(mkk2), file=paste(FILE_NAME, "genome.ucsc.gse.KEGG2.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
   
}, silent = TRUE)




result = try({
    
mkk3 <- enrichMKEGG(gene = GENES,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)
  
dotplot(mkk3)
ggsave( paste(FILE_NAME, "genome.ucsc.kegg3.png", sep="."), 
   width = 2800,
   height = 2800,
   units = "px",
   dpi = 300 ) 
    
write.table(as.data.frame(mkk3), file=paste(FILE_NAME, "genome.ucsc.module.KEGG3.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
     
}, silent = TRUE)



# WikiPathways analysis
# over-representation analysis 
# gene set enrichment analysis

result = try({
    
wp1 = enrichWP(GENES, organism = "Homo sapiens") 

dotplot(wp1)
ggsave( paste(FILE_NAME, "genome.ucsc.wiki1.png", sep="."), 
   width = 2800,
   height = 2800,
   units = "px",
   dpi = 300 ) 
    
write.table(as.data.frame(wp1), file=paste(FILE_NAME, "genome.ucsc.enrich.Wiki1.txt", sep="."),
   append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    
}, silent = TRUE)



result = try({
    
wp2 = gseWP(GENES, organism = "Homo sapiens")
dotplot(wp2)

ggsave( paste(FILE_NAME, "genome.ucsc.wiki2.png", sep="."), 
   width = 2800,
   height = 2800,
   units = "px",
   dpi = 300 ) 
    
write.table(as.data.frame(wp2), file=paste(FILE_NAME, "genome.ucsc.gse.Wiki2.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    
}, silent = TRUE)



# Over-representation analysis for the network of cancer genes (NCG)

result = try({
  
ncg <- enrichNCG(GENES) 
head(ncg)
if (dim(as.data.frame(ncg))[1] != 0)
{
    
dotplot(ncg)  
    
ggsave( paste(FILE_NAME, "genome.ucsc.enrich.cancer.genes.png", sep="."), 
    width = 2800,
    height = 2800,
    units = "px",
    dpi = 300 ) 
    
write.table(as.data.frame(ncg), file=paste(FILE_NAME, "genome.ucsc.enrich.NCG.txt", sep="."),
   append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}
    
}, silent = TRUE)



# Disease over-representation analysis
result = try({
  
disease <- enrichDO(gene = GENES,
              ont           = "DO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              universe      = GENES,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.1,
              readable      = FALSE)
    
if (dim(as.data.frame(disease))[1] != 0)
{
  
dotplot(disease) 
ggsave( paste(FILE_NAME, "genome.ucsc.enrich.DISEASE.genes.png", sep="."),
  width = 2800,
  height = 2800,
  units = "px",
  dpi = 300 ) 
  
write.table(as.data.frame(disease), file=paste(FILE_NAME, "genome.ucsc.enrich.DISEASE.txt", sep="."),
            append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}
    
}, silent = TRUE)



result = try({

m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
    dplyr::select(gs_name, entrez_gene)
    
em <- enricher(GENES, TERM2GENE=m_t2g)
dotplot(em)  
ggsave(paste(FILE_NAME, "genome.ucsc.enrich.MSIGDBR.C6.png", sep="."), 
    width = 2000,
   height = 2000,
   units = "px",
  dpi = 300 ) 
   
write.table(as.data.frame(em), file=paste(FILE_NAME, "genome.ucsc.enrich.MSIGDBR.txt", sep="."),
     append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    
em2 <- pairwise_termsim(em)
    
treeplot(em2) 
ggsave(paste(FILE_NAME, "genome.ucsc.enrich.MSIGDBR.C6.tree.png", sep="."), 
    width = 2000,
   height = 4000,
  units = "px",
  dpi = 300 ) 
   
}, silent = TRUE)



result = try({

m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
                 dplyr::select(gs_name, entrez_gene)
   
em <- enricher(GENES, TERM2GENE=m_t2g)

dotplot(em)

ggsave(paste(FILE_NAME, "genome.ucsc.enrich.MSIGDBR.H.png", sep="."), 
   width = 2000,
   height = 4000,
    units = "px",
   dpi = 300 ) 
   
write.table(as.data.frame(em), file=paste(FILE_NAME, "genome.ucsc.enrich.MSIGDBR.H.txt", sep="."),
    append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	

em2 <- pairwise_termsim(em)

treeplot(em2) 

ggsave(paste(FILE_NAME, "genome.ucsc.enrich.MSIGDBR.H.tree.png", sep="."), 
   width = 2000,
    height = 4000,
    units = "px",
   dpi = 300 ) 
    
}, silent = TRUE)