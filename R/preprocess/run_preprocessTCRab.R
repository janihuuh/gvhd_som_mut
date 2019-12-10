
message("Preprocess TCRab...")
## Read in the sequencing results; the most important files are all_contig and clonotypes -files

t2013_contig     <- read.delim("data/TCRseq/2013/all_contig_annotations.csv", sep = ",")
t2017_contig     <- read.delim("data/TCRseq/2017/all_contig_annotations.csv", sep = ",")
fimm_contig      <- read.delim("data/TCRseq/fimm_reanalysis/all_contig_annotations.csv", sep = ",")

t2013_clonotype  <- read.delim("data/TCRseq/2013/clonotypes.csv", sep = ",")
t2017_clonotype  <- read.delim("data/TCRseq/2017/clonotypes.csv", sep = ",")



## Pre-process the data. The function produces two files, barcoded and clonotyped
t2013_contig <- t2013_contig %>% dplyr::mutate(cdr3_aa = cdr3)
t2017_contig <- t2017_contig %>% dplyr::mutate(cdr3_aa = cdr3)

preprocess_10X_TCR(contig_file = t2013_contig, clonotype_file = t2013_contig, prefix = "data/TCRseq/preprocessed/2013")
preprocess_10X_TCR(contig_file = t2017_contig, clonotype_file = t2017_contig, prefix = "data/TCRseq/preprocessed/2017")



## Read in the just processed files
t2013_clonotyped <- read.delim("data/TCRseq/preprocessed/2013_clonotyped.txt")
t2017_clonotyped <- read.delim("data/TCRseq/preprocessed/2017_clonotyped.txt")

t2013_barcode    <- read.delim("data/TCRseq/preprocessed/2013_barcoded.txt") %>% mutate(barcode_uniq = paste0("2013_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
t2017_barcode    <- read.delim("data/TCRseq/preprocessed/2017_barcoded.txt") %>% mutate(barcode_uniq = paste0("2017_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))



t2013_clonotyped <- read.delim("data/TCRseq/preprocessed/s1a1_clonotyped.txt")
t2017_clonotyped <- read.delim("data/TCRseq/preprocessed/s1a2_clonotyped.txt")

t2013_barcode    <- read.delim("data/TCRseq/preprocessed/s1a1_barcoded.txt") %>% mutate(barcode_uniq = paste0("2013_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))
t2017_barcode    <- read.delim("data/TCRseq/preprocessed/s1a2_barcoded.txt") %>% mutate(barcode_uniq = paste0("2017_", barcode)) %>% mutate(barcode_uniq = substr(barcode_uniq, 1, nchar(barcode_uniq) - 2))



## Make new clonotype id:s for each of the patient profiled
gvhd_tot_barcode <- rbind(t2013_barcode, t2017_barcode) %>% newClonotype_df() %>% mutate(new_clonotypes_id = paste0("gvhd_", new_clonotypes_id))

## Write down results
write.table(gvhd_tot_barcode, "data/TCRseq/preprocessed/fgvhd_tcrab.txt",  sep = "\t", row.names = F, quote = F)
gvhd_tot_barcode[!duplicated(gvhd_tot_barcode$new_clonotypes_id), ] %>% write.table("data/TCRseq/preprocessed/gvhd_uniq_tcrab.txt", sep = "\t", row.names = F, quote = F)
