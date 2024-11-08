## ----install_G4SNVHunter, eval = FALSE----------------------------------------
#  
#  # Currently can only be installed by Github
#  devtools::install_github("rongxinzh/G4SNVHunter")
#  

## ----load_G4SNVHunter, message = FALSE----------------------------------------

library(G4SNVHunter)


## ----install_required_pkg, eval = FALSE---------------------------------------
#  
#  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
#      BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
#  }
#  
#  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
#      BiocManager::install("GenomicRanges")
#  }
#  
#  if (!requireNamespace("DT", quietly = TRUE)) {
#    install.packages("DT")
#  }
#  
#  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
#      BiocManager::install("rtracklayer")
#  }
#  

## ----load_required_pkg, message = FALSE, warning = FALSE----------------------

library(BSgenome.Hsapiens.UCSC.hg19)

library(GenomicRanges)

library(DT)

library(rtracklayer)


## ----load_seq_from_df_using_loadSequence--------------------------------------

seq_df <- data.frame(chr = c("seq1", "seq2"),
                     sequence = c(paste0(rep("G", 100), collapse = ""), 
                                  paste0(rep("A", 100), collapse = "")))
seq <- loadSequence(genome_seq = seq_df)


## ----load_seq_from_fa_using_loadSequence--------------------------------------

# File path to the sequences in fasta format
fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
seq <- loadSequence(seq_path = fa_path)


## ----load_seq_from_txt_using_loadSequence-------------------------------------

# File path to the sequences in txt format
txt_path <- system.file("extdata", "seq.txt", package = "G4SNVHunter")
seq <- loadSequence(seq_path = txt_path)


## ----install_hg19_refseq, message = FALSE, warning = FALSE--------------------

# Load sequence for chromosome 21 (hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
chr21_seq <- DNAStringSet(hg19$chr21)
# Chromosome names are needed for analysis
names(chr21_seq) <- "chr21"


## ----load_SNV, message = FALSE------------------------------------------------

# Path to SNPs
snp_path <- system.file("extdata", "snp.txt", package = "G4SNVHunter")
# Load SNPs into memory
snp <- read.table(snp_path, sep = "\t", header = FALSE)
# Convert snp to GRanges
snp <- GRanges(seqnames = snp$V1,
               ranges = IRanges(start = snp$V2, width = 1),
               rsid = snp$V3,
               ref = snp$V4,
               alt = snp$V5)
print(snp)


## ----check_SNVs---------------------------------------------------------------

gr <- GRanges(
  seqnames = Rle("seq1"),
  ranges = IRanges(c(100, 200, 300), width=1),
    ref = c("A", "C", "G"),
    alt = c("T", "T", "A")
  )
# Check width ('w'), ref ('r'), and alt ('a')
# Returns TRUE if it passed the validation
checkSNV(gr, mode = "wra", ref_col = "ref", alt_col = "alt")


## ----predict_G4s--------------------------------------------------------------

# Sequence file in fasta file format 
fa_path <- system.file("extdata", "seq.fa", package = "G4SNVHunter")
# Load sequences
seq <- loadSequence(seq_path = fa_path)

# Predict G4s
G4_detected <- G4HunterDetect(seq)


## ----print_G4s, warning = FALSE, message = FALSE------------------------------

datatable(as.data.frame(G4_detected), options = list(scrollX = TRUE))


## ----G4_score_maxscore, echo=FALSE--------------------------------------------

knitr::include_graphics("images/G4SNVHunter_01.svg")


## ----predict_G4s_parameters---------------------------------------------------

# Predict G4 by customizing other parameters
G4_detected <- G4HunterDetect(seq, threshold = 1.5, window_size = 20)


## ----export_G4_as_csv, eval = FALSE-------------------------------------------
#  
#  # export as csv format
#  write.csv(as.data.frame(G4_detected), "/path/to/your/file.csv")
#  

## ----export_G4_as_other_formats, eval = FALSE---------------------------------
#  
#  # export as bed format
#  export(G4_detected, "/path/to/your/file.bed", format = "bed")
#  # export as bigwig format
#  export(G4_detected, "/path/to/your/file.bw", format = "bigWig")
#  

## ----plot_G4, message = FALSE, fig.width = 8.5, fig.height = 8----------------

plotG4Info(G4_detected)


## ----SNV_effect_interpretation, echo = FALSE----------------------------------

knitr::include_graphics("images/G4SNVHunter_02.svg")


## ----load_example_gr_data-----------------------------------------------------

data(snp_gr)

data(snv_gr)


## ----predict_chr21_G4s--------------------------------------------------------

# Predict the G4s in human chr 21 (hg19)
chr21_G4 <- G4HunterDetect(chr21_seq)


## ----SNV_effects_s_mode, message = FALSE--------------------------------------

snp_eff <- SNVImpactG4(chr21_G4, 
                       snp_gr, 
                       alt_col = "alt")


## ----SNV_effects_s_mode_header3-----------------------------------------------

datatable(as.data.frame(snp_eff[1:3]), options = list(scrollX = TRUE))


## ----SNV_effects_m_mode, message = FALSE--------------------------------------

# Column names of the Sample ID and SNV ID must be specified
snv_eff <- SNVImpactG4(chr21_G4, 
                       snv_gr, 
                       alt_col = "alt", 
                       mode = "m", 
                       sampleid_col = "sampleid", 
                       snvid_col = "snv_id")


## ----SNV_effects_m_mode_header3-----------------------------------------------

datatable(as.data.frame(snv_eff[1:3]), options = list(scrollX = TRUE))


## ----SNV_effects_m_mode_example-----------------------------------------------

datatable(as.data.frame(snv_eff[528]), options = list(scrollX = TRUE))


## ----filter_SNV_impact--------------------------------------------------------

filtered_snv_eff <- filterSNVImpact(snv_eff,
                                    mut_score_threshold = 1.2,
                                    score_diff_threshold = -0.35)


## ----view_SNV_significant_impact----------------------------------------------

datatable(as.data.frame(filtered_snv_eff), options = list(scrollX = TRUE))


## ----eff_plot, fig.height = 3, message = FALSE--------------------------------

plotSNVImpact(snv_eff)


## ----eff_plot2, fig.height = 3, message = FALSE-------------------------------

plotSNVImpact(filtered_snv_eff)


## ----seqlogo_plot, fig.height = 8, fig.width = 10, message = FALSE------------

top5_snv_eff <- filtered_snv_eff[order(filtered_snv_eff$score.diff)[1:5]]
plotImpactSeq(top5_snv_eff, ncol = 2)


## ----coherent_example, eval = FALSE-------------------------------------------
#  
#  # First step, predict G4s
#  chr21_G4 <- G4HunterDetect(chr21_seq)
#  
#  # Second step, evaluate the impact of SNVs on G4s in 's' mode
#  snv_eff <- SNVImpactG4(chr21_G4, snv_gr, alt_col = "alt")
#  
#  # Filter the results under default parameters
#  filtered_snv_eff <- filterSNVImpact(snv_eff, mut_score_threshold = 1.2)
#  
#  # export as csv format
#  write.csv(as.data.frame(filtered_snv_eff), "/path/to/your/file.csv")
#  

## ----session_info-------------------------------------------------------------

sessionInfo()


