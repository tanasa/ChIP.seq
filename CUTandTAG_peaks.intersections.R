library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(clusterProfiler)
library("ggplot2")
library("ChIPpeakAnno")
library(gprofiler2)
library(diffloop)
library(GenomicRanges)
library(Repitools)
library(DOSE)
library(grid)
library("EnsDb.Hsapiens.v86")
library("enrichplot")
library(msigdbr)
library("rtracklayer")
library(VennDiagram)
library(tidyverse)
library(Gviz)
library(trackViewer)

# in the script below, the code extends the peaks with +/- 1kb, 
# merges the overlaping peaks, computes the Venn Diagrams
# the locations of the peaks can be visualized with Gviz

###############################################################
###############################################################

h3k18.rep1 = "H3K18a_rep1_L-tron.t2.tiny.enrich.0.05.bed"
h3k18.rep1.name = basename(h3k18.rep1)
  
df.h3k18.rep1 = read.table(h3k18.rep1, header = F, sep="\t", stringsAsFactors = F)
colnames(df.h3k18.rep1)[1] =  "chr"
colnames(df.h3k18.rep1)[2] =  "start"
colnames(df.h3k18.rep1)[3] =  "end"
colnames(df.h3k18.rep1)[4] =  "score"

head(df.h3k18.rep1, 2)
dim(df.h3k18.rep1)

gr.h3k18.rep1 = import("H3K18a_rep1_L-tron.t2.tiny.enrich.0.05.bed", format="bed")
head(gr.h3k18.rep1, 2)

###############################################################
###############################################################

h3k18.rep2 = "H3K18a_rep2_L-tron.t2.tiny.enrich.0.05.bed"
h3k18.rep2.name = basename(h3k18.rep2)
  
df.h3k18.rep2 = read.table(h3k18.rep2, header = F, sep="\t", stringsAsFactors = F)
colnames(df.h3k18.rep2)[1] =  "chr"
colnames(df.h3k18.rep2)[2] =  "start"
colnames(df.h3k18.rep2)[3] =  "end"
colnames(df.h3k18.rep2)[4] =  "score"

head(df.h3k18.rep2, 2)
dim(df.h3k18.rep2)

gr.h3k18.rep2 = import("H3K18a_rep2_L-tron.t2.tiny.enrich.0.05.bed", format="bed")
head(gr.h3k18.rep2, 2)

###############################################################
###############################################################

h3k27.rep1 = "H3K27me3_rep1_L-tron.t2.tiny.enrich.0.05.bed"
h3k27.rep1.name = basename(h3k27.rep1)
  
df.h3k27.rep1 = read.table(h3k27.rep1, header = F, sep="\t", stringsAsFactors = F)
colnames(df.h3k27.rep1)[1] =  "chr"
colnames(df.h3k27.rep1)[2] =  "start"
colnames(df.h3k27.rep1)[3] =  "end"
colnames(df.h3k27.rep1)[4] =  "score"

head(df.h3k27.rep1, 2)
dim(df.h3k27.rep1)

gr.h3k27.rep1 = import("H3K27me3_rep1_L-tron.t2.tiny.enrich.0.05.bed", format="bed")
head(gr.h3k27.rep1, 2)

###############################################################
###############################################################

h3k27.rep2 = "H3K27me3_rep2_L-tron.t2.tiny.enrich.0.05.bed"
h3k27.rep2.name = basename(h3k27.rep2)
  
df.h3k27.rep2 = read.table(h3k27.rep2, header = F, sep="\t", stringsAsFactors = F)
colnames(df.h3k27.rep2)[1] =  "chr"
colnames(df.h3k27.rep2)[2] =  "start"
colnames(df.h3k27.rep2)[3] =  "end"
colnames(df.h3k27.rep2)[4] =  "score"

head(df.h3k27.rep2, 2)
dim(df.h3k27.rep2)

gr.h3k27.rep2 = import("H3K27me3_rep2_L-tron.t2.tiny.enrich.0.05.bed", format="bed")
head(gr.h3k27.rep2, 2)

###############################################################
###############################################################

h3k4.rep1 = "pos_rep1_control_H3K4me3_L-tron.t2.tiny.enrich.1.bed"
h3k27.rep1.name = basename(h3k4.rep1)
  
df.h3k4.rep1 = read.table(h3k4.rep1, header = F, sep="\t", stringsAsFactors = F)
colnames(df.h3k4.rep1)[1] =  "chr"
colnames(df.h3k4.rep1)[2] =  "start"
colnames(df.h3k4.rep1)[3] =  "end"
colnames(df.h3k4.rep1)[4] =  "score"

head(df.h3k4.rep1, 2)
dim(df.h3k4.rep1)

gr.h3k4.rep1 = import("pos_rep1_control_H3K4me3_L-tron.t2.tiny.enrich.1.bed", format="bed")
head(gr.h3k4.rep1, 2)

###############################################################
###############################################################

h3k4.rep2 = "pos_rep2_control_H3K4me3_L-tron.t2.tiny.enrich.1.bed"
h3k4.rep2.name = basename(h3k4.rep2)
  
df.h3k4.rep2 = read.table(h3k4.rep2, header = F, sep="\t", stringsAsFactors = F)
colnames(df.h3k4.rep2)[1] =  "chr"
colnames(df.h3k4.rep2)[2] =  "start"
colnames(df.h3k4.rep2)[3] =  "end"
colnames(df.h3k4.rep2)[4] =  "score"

head(df.h3k4.rep2, 2)
dim(df.h3k4.rep2)

gr.h3k4.rep2 = import("pos_rep2_control_H3K4me3_L-tron.t2.tiny.enrich.1.bed", format="bed")
head(gr.h3k4.rep2, 2)

###############################################################
###############################################################

length(gr.h3k18.rep1)
length(gr.h3k18.rep2)
length(gr.h3k27.rep1)
length(gr.h3k27.rep2)
length(gr.h3k4.rep1)
length(gr.h3k4.rep2)

###############################################################
###############################################################

left_extension <- 1000
right_extension <- 1000

extend_and_merge_peaks = function(gr_peaks, left_extension, right_extension) {
    
   extended_peaks <- GRanges(seqnames = seqnames(gr_peaks),
                             ranges = IRanges(start = start(gr_peaks) - left_extension,
                                               end = end(gr_peaks) + right_extension))
    
   merged_peaks <- GenomicRanges::reduce(extended_peaks)
}


gr.h3k18.rep1.1kem = extend_and_merge_peaks(gr.h3k18.rep1, left_extension, right_extension)
gr.h3k18.rep2.1kem = extend_and_merge_peaks(gr.h3k18.rep2, left_extension, right_extension)

gr.h3k27.rep1.1kem = extend_and_merge_peaks(gr.h3k27.rep1, left_extension, right_extension)
gr.h3k27.rep2.1kem = extend_and_merge_peaks(gr.h3k27.rep2, left_extension, right_extension)

gr.h3k4.rep1.1kem = extend_and_merge_peaks(gr.h3k4.rep1, left_extension, right_extension)
gr.h3k4.rep2.1kem = extend_and_merge_peaks(gr.h3k4.rep2, left_extension, right_extension)

# another way to extend the peaks

# gr.h3k18.rep1.1k = gr.h3k18.rep1 + 1000
# gr.h3k18.rep1.1k1m = GenomicRanges::reduce(gr.h3k18.rep1.1k)

# gr.h3k18.rep2.1k = gr.h3k18.rep2 + 1000
# gr.h3k18.rep2.1k1m = GenomicRanges::reduce(gr.h3k18.rep2.1k)

# gr.h3k27.rep1.1k = gr.h3k27.rep1 + 1000
# gr.h3k27.rep1.1k1m = GenomicRanges::reduce(gr.h3k27.rep1.1k)

# gr.h3k27.rep2.1k = gr.h3k27.rep2 + 1000
# gr.h3k27.rep2.1k1m = GenomicRanges::reduce(gr.h3k27.rep2.1k)

# gr.h3k4.rep1.1k = gr.h3k4.rep1 + 1000
# gr.h3k4.rep1.1k1m = GenomicRanges::reduce(gr.h3k4.rep1.1k)

# gr.h3k4.rep2.1k = gr.h3k4.rep2 + 1000
# gr.h3k4.rep2.1k1m = GenomicRanges::reduce(gr.h3k4.rep2.1k)

# the number of peaks after extending the peaks with 1kb on each side, and merging the peaks

length(gr.h3k18.rep1)
# length(gr.h3k18.rep1.1k1m)
length(gr.h3k18.rep1.1kem)

length(gr.h3k27.rep1)
# length(gr.h3k27.rep1.1k1m)
length(gr.h3k27.rep1.1kem)

length(gr.h3k4.rep1)
# length(gr.h3k4.rep1.1k1m)
length(gr.h3k4.rep1.1kem)

length(gr.h3k18.rep2)
# length(gr.h3k18.rep2.1k1m)
length(gr.h3k18.rep2.1kem)

length(gr.h3k27.rep2)
# length(gr.h3k27.rep2.1k1m)
length(gr.h3k27.rep2.1kem)

length(gr.h3k4.rep2)
# length(gr.h3k4.rep2.1k1m)
length(gr.h3k4.rep2.1kem)

###############################################################
###############################################################

# we would like to keep only the chromosomes 1 to 22, X, and Y

selected_chromosomes <- c(paste0(1:22), "X")
gr.h3k18.rep1.1kem.filter = keepSeqlevels(gr.h3k18.rep1.1kem, selected_chromosomes)
# gr.h3k18.rep1.1kem.filter

length(gr.h3k18.rep1.1kem)
length(gr.h3k18.rep1.1kem.filter)
gr.h3k18.rep1.1kem = gr.h3k18.rep1.1kem.filter

selected_chromosomes <- c(paste0(1:22), "X")
gr.h3k18.rep2.1kem.filter = keepSeqlevels(gr.h3k18.rep2.1kem, selected_chromosomes)
# gr.h3k18.rep2.1kem.filter

length(gr.h3k18.rep2.1kem)
length(gr.h3k18.rep2.1kem.filter)
gr.h3k18.rep2.1kem = gr.h3k18.rep2.1kem.filter

selected_chromosomes <- c(paste0(1:22), "X")
gr.h3k27.rep1.1kem.filter = keepSeqlevels(gr.h3k27.rep1.1kem, selected_chromosomes, pruning.mode="coarse")
# gr.h3k27.rep1.1kem.filter

length(gr.h3k27.rep1.1kem)
length(gr.h3k27.rep1.1kem.filter)
gr.h3k27.rep1.1kem = gr.h3k27.rep1.1kem.filter

selected_chromosomes <- c(paste0(1:22), "X")
gr.h3k27.rep2.1kem.filter = keepSeqlevels(gr.h3k27.rep2.1kem, selected_chromosomes, pruning.mode="coarse")
# gr.h3k27.rep2.1kem.filter

length(gr.h3k27.rep2.1kem)
length(gr.h3k27.rep2.1kem.filter)
gr.h3k27.rep2.1kem = gr.h3k27.rep2.1kem.filter

selected_chromosomes <- c(paste0(1:22), "X")
gr.h3k4.rep1.1kem.filter = keepSeqlevels(gr.h3k4.rep1.1kem, selected_chromosomes, pruning.mode="coarse")
# gr.h3k4.rep1.1kem.filter

length(gr.h3k4.rep1.1kem)
length(gr.h3k4.rep1.1kem.filter)
gr.h3k4.rep1.1kem = gr.h3k4.rep1.1kem.filter

selected_chromosomes <- c(paste0(1:22), "X", "Y")
gr.h3k4.rep2.1kem.filter = keepSeqlevels(gr.h3k4.rep2.1kem, selected_chromosomes, pruning.mode="coarse")
# gr.h3k4.rep2.1kem.filter

length(gr.h3k4.rep2.1kem)
length(gr.h3k4.rep2.1kem.filter)
gr.h3k4.rep2.1kem = gr.h3k4.rep2.1kem.filter

# the number of the peaks : H3K18ac
cat("the number of H3K18ac peaks in each replicate and in the overlap :")

length(gr.h3k18.rep1.1kem)
length(gr.h3k18.rep2.1kem)

gr.h3k18.overlap <- findOverlaps(gr.h3k18.rep1.1kem, gr.h3k18.rep2.1kem)
length(gr.h3k18.overlap)
 
# makeVennDiagram(gr.h3k18.overlap)

gr.h3k18.venn = draw.pairwise.venn(
  area1 = length(gr.h3k18.rep1.1kem),
  area2 = length(gr.h3k18.rep2.1kem),
  cross.area = length(gr.h3k18.overlap), # ov,
  category = c("rep1", "rep2"),
  fill = c("orange", "darkgreen"),
  cat.cex = 0.9,
  
  # Add these lines to show numbers:
  # 1. Include label texts:
  
  labels = c(paste0("rep1 (", length(gr.h3k18.rep1.1kem), ")"),
            paste0("rep2 (", length(gr.h3k18.rep2.1kem), ")"),
            paste0("Overlap (", ov, ")")),
  
  # 2. Specify text positions and justification:
  text.col = "black",
  text.font = 2,
  text.pos = c("top", "top", "center"),
  text.adjust = c(0.5, 0.5, 0.5)
)

grid.draw(gr.h3k18.venn)
grid.text("Venn Diagram of H3K18ac Peaks in Replicates 1 and 2", 
           x = 0.5, y = 0.9, gp = gpar(col = "black", font = 2))
# grid.newpage()

# the number of the peaks : H3K27me3
cat("the number of H3K27me3 peaks in each replicate and in the overlap :")

length(gr.h3k27.rep1.1kem)
length(gr.h3k27.rep2.1kem)

gr.h3k27.overlap <- findOverlaps(gr.h3k27.rep1.1kem, gr.h3k27.rep2.1kem)
length(gr.h3k27.overlap)

# makeVennDiagram(gr.h3k27.overlap)

gr.h3k27.venn = draw.pairwise.venn(
  area1 = length(gr.h3k27.rep1.1kem),
  area2 = length(gr.h3k27.rep2.1kem),
  cross.area = length(gr.h3k27.overlap), # ov,
  category = c("rep1", "rep2"),
  fill = c("orange", "darkgreen"),
  cat.cex = 0.9,
  
  # Add these lines to show numbers:
  # 1. Include label texts:
  
  labels = c(paste0("rep1 (", length(gr.h3k27.rep1.1kem), ")"),
            paste0("rep2 (", length(gr.h3k27.rep2.1kem), ")"),
            paste0("Overlap (", ov, ")")),
  
  # 2. Specify text positions and justification:
  text.col = "black",
  text.font = 2,
  text.pos = c("top", "top", "center"),
  text.adjust = c(0.5, 0.5, 0.5)
)

grid.draw(gr.h3k27.venn)
grid.text("Venn Diagram of H3K27me3 Peaks in Replicates 1 and 2", 
           x = 0.5, y = 0.9, gp = gpar(col = "black", font = 2))
# grid.newpage())


###############################################################
###############################################################

# the number of the peaks : H3K4me3
cat("the number of H3K4me3 peaks in each replicate and in the overlap :")

length(gr.h3k4.rep1.1kem)
length(gr.h3k4.rep2.1kem)

gr.h3k4.overlap <- findOverlaps(gr.h3k4.rep1.1kem, gr.h3k4.rep2.1kem)
length(gr.h3k4.overlap)

# makeVennDiagram(gr.h3k4.overlap)

gr.h3k4.venn = draw.pairwise.venn(
  area1 = length(gr.h3k4.rep1.1kem),
  area2 = length(gr.h3k4.rep2.1kem),
  cross.area = length(gr.h3k4.overlap), # ov,
  category = c("rep1", "rep2"),
  fill = c("orange", "darkgreen"),
  cat.cex = 0.9,
  
  # Add these lines to show numbers:
  # 1. Include label texts:
  
  labels = c(paste0("rep1 (", length(gr.h3k4.rep1.1kem), ")"),
            paste0("rep2 (", length(gr.h3k4.rep2.1kem), ")"),
            paste0("Overlap (", ov, ")")),
  
  # 2. Specify text positions and justification:
  text.col = "black",
  text.font = 2,
  text.pos = c("top", "top", "center"),
  text.adjust = c(0.5, 0.5, 0.5)
)

grid.draw(gr.h3k4.venn)
grid.text("Venn Diagram of H3K4me3 Peaks in Replicates 1 and 2", 
           x = 0.5, y = 0.9, gp = gpar(col = "black", font = 2))
# grid.newpage())

###############################################################
###############################################################

# Visualize the peaks with Gviz

GENOME <- "hg38"
txdb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
CHR ="chr1"

# chromosome(obj) <- 'chrId' to change the active chromosome 

# gr.h3k4.rep1.1kem.track 
# gr.h3k4.rep2.1kem.track 

# GENES <- trackViewer::geneModelFromTxdb("TxDb.Hsapiens.UCSC.hg38.knownGene", org.Hs.eg.db )

GENES = GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene, 
                 transcriptAnnotation = "symbol", 
                 geneSymbol = TRUE,
                 name = "Genes")



gtrack <- GenomeAxisTrack(genome=GENOME)
itrack <- IdeogramTrack(genome = GENOME)
strack = SequenceTrack(Hsapiens) # the sequence tracks

atrack_gr.h3k4.rep1.1kem <- AnnotationTrack(gr.h3k4.rep1.1kem, name = "h3k4.rep1", genome = GENOME)
atrack_gr.h3k4.rep2.1kem <- AnnotationTrack(gr.h3k4.rep2.1kem, name = "h3k4.rep1", genome = GENOME)
atrack_gr.h3k27.rep1.1kem <- AnnotationTrack(gr.h3k27.rep1.1kem, name = "h3k27.rep1", genome = GENOME)
atrack_gr.h3k27.rep2.1kem <- AnnotationTrack(gr.h3k27.rep2.1kem, name = "h3k27.rep1", genome = GENOME)
atrack_gr.h3k18.rep1.1kem <- AnnotationTrack(gr.h3k18.rep1.1kem, name = "h3k18.rep1", genome = GENOME)
atrack_gr.h3k18.rep2.1kem <- AnnotationTrack(gr.h3k18.rep2.1kem, name = "h3k18.rep1", genome = GENOME)

options(repr.plot.width = 16, repr.plot.height = 6) 
plotTracks(list(gtrack,
                itrack, 
                GENES,
                atrack_gr.h3k4.rep1.1kem,
                atrack_gr.h3k4.rep2.1kem,
                atrack_gr.h3k27.rep1.1kem,
                atrack_gr.h3k27.rep2.1kem,
                atrack_gr.h3k18.rep1.1kem,
                atrack_gr.h3k18.rep2.1kem ),
           add53 = TRUE, 
           add35 = TRUE,
           littleTicks = TRUE,
         # labelPos="below",
           chromosome = "chr1", from = 4000000, to = 5000000)



plotTracks(list(gtrack,
                itrack, 
                strack, 
                GENES,
                atrack_gr.h3k4.rep1.1kem,
                atrack_gr.h3k4.rep2.1kem,
                atrack_gr.h3k27.rep1.1kem,
                atrack_gr.h3k27.rep2.1kem,
                atrack_gr.h3k18.rep1.1kem,
                atrack_gr.h3k18.rep2.1kem ),
           add53 = TRUE, 
           add35 = TRUE,
           littleTicks = TRUE,
         # labelPos="below",
           chromosome = "chr1", from = 4000000, to = 5000000)

plotTracks(list(gtrack,
                itrack, 
                strack, 
                GENES,
                atrack_gr.h3k4.rep1.1kem,
                atrack_gr.h3k4.rep2.1kem,
                atrack_gr.h3k27.rep1.1kem,
                atrack_gr.h3k27.rep2.1kem,
                atrack_gr.h3k18.rep1.1kem,
                atrack_gr.h3k18.rep2.1kem ),
           add53 = TRUE, 
           add35 = TRUE,
           littleTicks = TRUE,
         # labelPos="below",
           chromosome = "chr1", from = 4999900, to = 5000000)


###############################################################
###############################################################
###############################################################
