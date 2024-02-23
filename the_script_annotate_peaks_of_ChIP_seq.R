suppressPackageStartupMessages({
  
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library("ggplot2")
library("ChIPpeakAnno")
library(gprofiler2)
library(GenomicRanges)
library(Repitools)
library(DOSE)
library("EnsDb.Hsapiens.v86")
library("org.Hs.eg.db")
library("enrichplot")
library(msigdbr)
library(diffloop)
  
})

##########################################################
##########################################################

# DIRECTORY="//"
# setwd(DIRECTORY)
# FILE = "H3K18a_rep1_L-tron.t2.tiny.enrich.0.05.bed"

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No arguments provided. Please provide an argument.")
}

# Print the argument
cat("The provided argument is:", args[1], "\n")

FILE = args[1]
FILE_NAME = basename(FILE)
  
########################################################## a file in the format :
########################################################## CHR STARt END SCORE

df.file=read.table(FILE, header=F, sep="\t", stringsAsFactors=F)

colnames(df.file)[1] =  "chromosome"
colnames(df.file)[2] =  "start"
colnames(df.file)[3] =  "end"
colnames(df.file)[4] =  "score"

##########################################################
##########################################################

# We can read the peak files in one of these three ways : 
# peak = readPeakFile(FILE, as="GRanges")
# head(peak)
# We use peak2 in downstream analysis, although peak3 can be used as well.
peak2 = toGRanges(df.file)
# peak3 = makeGRangesFromDataFrame(df.file,keep.extra.columns=TRUE) 

# to make a Coverage Profile
# png(paste(FILE_NAME, "coverage.plot.png", sep="."), width = 2000, height = 1600, units = "px")
# covplot(peak2, weightCol="score")
# dev.off()

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

##########################################################
##########################################################

# Annotations based on Ensembl database
# These are almost identical with UCSC annotations.
# edb <- EnsDb.Hsapiens.v86
# seqlevelsStyle(edb) = "Ensembl"
# peakAnno_edb <- annotatePeak(peak2, tssRegion=c(-3000, 3000),
#                             TxDb=edb, annoDb="org.Hs.eg.db",   
#                             addFlankGeneInfo = TRUE,
   ## flankDistance = 1000000,
#  sameStrand = FALSE,
#  ignoreOverlap = FALSE,
#  ignoreUpstream = FALSE,
#  ignoreDownstream = FALSE,
#  overlap = "all",
#  verbose = TRUE)
# x = as.data.frame(peakAnno.edb@anno@elementMetadata)
# peakAnno_edb_GR <- as.GRanges(peakAnno_edb)
# peakAnno_edb_DF <- as.data.frame(peakAnno_edb)
# dim(df.file)
# length(peakAnno_edb_GR) 
# dim(peakAnno_edb_DF)[1]
# length(unique(peakAnno_edb_DF$SYMBOL))
# write.table(peakAnno_ucsc_DF, file=paste(FILE_NAME, "genome.ensembl.annotations.txt", sep="."),
#							    append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##########################################################
##########################################################

# Annotation Pie
plotAnnoPie(peakAnno_ucsc)
png(paste(FILE_NAME, "genome.ucsc.annotations.pie.png", sep="."), width = 400, height = 400, units = "px")
plotAnnoPie(peakAnno_ucsc)
dev.off()
# png(paste(FILE_NAME, "ensembl.annotations.pie.png", sep="."), width = 1000, height = 1000, units = "px")
# plotAnnoPie(peakAnno_edb)
# dev.off()
# png(paste(FILE_NAME, "ucsc.annotations.bar.png", sep="."), width = 600, height = 600, units = "px")
# plotAnnoBar(peakAnno_ucsc)
# dev.off()
# png(paste(FILE_NAME, "ensembl.annotations.bar.png", sep="."), width = 1000, height = 1000, units = "px")
# plotAnnoBar(peakAnno_edb)
# dev.off()
# png(paste(FILE_NAME, "ucsc.annotations.vennpie.png", sep="."), width = 600, height = 600, units = "px")
# vennpie(peakAnno_ucsc)
# dev.off()
# png(paste(FILE_NAME, "ensembl.annotations.vennpie.png", sep="."), width = 1000, height = 1000, units = "px")
# vennpie(peakAnno_edb)
# dev.off()

##########################################################
##########################################################

# Functional enrichment analysis
# Selecting the genes that we will work with
length(as.data.frame(peakAnno_ucsc)$geneId)
GENES = unique(as.data.frame(peakAnno_ucsc)$geneId)
length(GENES)

##########################################################
##########################################################

result = try({
  
egoo <- enrichGO(gene    = GENES,
                OrgDb     = org.Hs.eg.db,
                # keyType  = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1, readable = TRUE)
# png(paste(FILE_NAME, "genome.ucsc.enrich.GO.BP.png", sep="."), width = 800, height = 800, units = "px")
dotplot(egoo)								 
# dev.off()
ggsave(paste(FILE_NAME, "genome.ucsc.enrich.GO.BP.png", sep="."), 
			   width = 2800,
			   height = 2800,
			   units = "px",
			   dpi = 300 )
write.table(as.data.frame(egoo), file=paste(FILE_NAME, "genome.ucsc.enrich.GO.BP.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)	
}, silent = TRUE)

##########################################################
##########################################################	

result = try({
  
egoo2 <- enrichGO(gene    = GENES,
                OrgDb     = org.Hs.eg.db,
                # keyType  = 'ENSEMBL',
                ont           = "MF",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1, readable = TRUE)
#png(paste(FILE_NAME, "genome.ucsc.enrich.GO.MF.png", sep="."), width = 800, height = 800, units = "px")
dotplot(egoo2)								 
# dev.off()
ggsave(paste(FILE_NAME, "genome.ucsc.enrich.GO.MF.png", sep="."), 
			   width = 2800,
			   height = 2800,
			   units = "px",
			   dpi = 300 ) 
write.table(as.data.frame(egoo2), file=paste(FILE_NAME, "genome.ucsc.enrich.GO.MF.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}, silent = TRUE)	

##########################################################
##########################################################

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
		   
# png(paste(FILE_NAME, "genome.ucsc.kegg1.png", sep="."), width = 800, height = 800, units = "px") 
dotplot(mkk1)
# dev.off()
ggsave( paste(FILE_NAME, "genome.ucsc.kegg1.png", sep="."), 
			   width = 2800,
			   height = 2800,
			   units = "px",
			   dpi = 300 ) 
write.table(as.data.frame(mkk1), file=paste(FILE_NAME, "genome.ucsc.enrich.KEGG1.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
				   							    
}, silent = TRUE)	

##########################################################
##########################################################

result = try({
mkk2 <- gseKEGG(gene = GENES,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1)
# png(paste(FILE_NAME, "genome.ucsc.kegg2.png", sep="."), width = 800, height = 800, units = "px") 
dotplot(mkk2)
# dev.off()
ggsave( paste(FILE_NAME, "genome.ucsc.kegg2.png", sep="."), 
			   width = 2800,
			   height = 2800,
			   units = "px",
			   dpi = 300 ) 
write.table(as.data.frame(mkk2), file=paste(FILE_NAME, "genome.ucsc.gse.KEGG2.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
				   							    
}, silent = TRUE)	

##########################################################
##########################################################

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

##########################################################
##########################################################

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

##########################################################
##########################################################

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

##########################################################
##########################################################

# Over-representation analysis for the network of cancer genes (NCG)
result = try({
  
ncg <- enrichNCG(GENES) 
head(ncg)
if (dim(as.data.frame(ncg))[1] != 0)
{
																 
# png(paste(FILE_NAME, "genome.ucsc.enrich.cancer.genes.png", sep="."), width = 800, height = 800, units = "px")
dotplot(ncg)								 
# dev.off()
ggsave( paste(FILE_NAME, "genome.ucsc.enrich.cancer.genes.png", sep="."), 
			   width = 2800,
			   height = 2800,
			   units = "px",
			   dpi = 300 ) 
write.table(as.data.frame(ncg), file=paste(FILE_NAME, "genome.ucsc.enrich.NCG.txt", sep="."),
								  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}
}, silent = TRUE)

##########################################################
##########################################################

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

##########################################################
##########################################################

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

##########################################################
##########################################################

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

##########################################################
##########################################################
