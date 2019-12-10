

filterContigFiles <- function(contig_file){
  
  # @ param:
  # input: needs to outputted by read.delim df (and not e.g. fread)
  
  # Filters contig files in such way that removes:
  ## - sequences with low confidence
  ## - non-productive
  ## - no cdr3 gathered
  ## - also could consider to remove non-TCR (i.e. BCR) reads
  
  require(dplyr)
  
  contig_file %>%
    filter(high_confidence == "True") %>%
    filter(chain %in% c("TRA", "TRB", "TRD", "TRG", "multi")) %>%
    filter(productive != "False") %>%
    filter(cdr3 != "None")
  
  
}

getVDJforFile <- function(file_temp){
  
  # @ params: file_temp = contig file, which is hopefully filtered
  # @ output: file with consensus V, D and J annotations
  
  get_vdj <- function(file_temp, raw_clonotype_id, chain){
    
    file_temp       <- filter(file_temp, chain == chain_temp)
    clonotype_index <- which(as.character(file_temp$raw_clonotype_id) %in% raw_clonotype_id)
    
    conv_temp       <- file_temp[clonotype_index, ]
    
    ## Get consensus cdr3aa, v, d and j genes and their respective frequencies
    conv_temp <- conv_temp %>%
      mutate(cdr3_nt = as.factor(cdr3_nt), cdr3_aa = as.factor(cdr3_aa), raw_clonotype_id = as.factor(raw_clonotype_id)) %>%
      mutate(cdr3_nt = droplevels(cdr3_nt), cdr3_aa = droplevels(cdr3_aa), raw_clonotype_id = droplevels(raw_clonotype_id))
    
    cdr3_nt_all    <- names(sort(table(conv_temp$cdr3_nt), decreasing = T))
    cdr3_nt        <- cdr3_nt_all[1]
    cdr3_nt_amount <- as.numeric(sort(table(conv_temp$cdr3_nt), decreasing = T)[1])
    cdr3_nt_freq   <- cdr3_nt_amount / sum(table(conv_temp$cdr3_nt))
    
    cdr3_aa_all    <- names(sort(table(conv_temp$cdr3_aa), decreasing = T))
    cdr3_aa        <- cdr3_aa_all[1]
    cdr3_aa_amouaa <- as.numeric(sort(table(conv_temp$cdr3_aa), decreasing = T)[1])
    cdr3_aa_freq   <- cdr3_aa_amouaa / sum(table(conv_temp$cdr3_aa))
    
    v_all           <- names(sort(table(conv_temp$v_gene), decreasing = T))
    v               <- v_all[1]
    v_amount        <- as.numeric(sort(table(conv_temp$v), decreasing = T)[1])
    v_freq          <- v_amount / sum(table(conv_temp$v))
    
    d               <- names(sort(table(conv_temp$d_gene), decreasing = T)[1])
    d_amount        <- as.numeric(sort(table(conv_temp$d), decreasing = T)[1])
    d_freq          <- d_amount / sum(table(conv_temp$d))
    
    j_all           <- names(sort(table(conv_temp$j_gene), decreasing = T))
    j               <- j_all[1]
    j_amount        <- as.numeric(sort(table(conv_temp$j), decreasing = T)[1])
    j_freq          <- j_amount / sum(table(conv_temp$j))
    
    
    # Combine the results
    tot <- data.frame(raw_clonotype_id, cdr3_nt, cdr3_aa, chain, v, d, j, cdr3_nt_freq, cdr3_aa_freq, v_freq, d_freq, j_freq)
    return(tot)
    
  }
  
  ## Analyse the data
  file_temp$chain        <- as.factor(file_temp$chain)
  file_temp$raw_clonotype_id <- as.factor(file_temp$raw_clonotype_id)
  
  file_temp$chain        <- droplevels(file_temp$chain)
  file_temp$raw_clonotype_id <- droplevels(file_temp$raw_clonotype_id)
  
  tcr_list <- list()
  i <- 1
  
  for(chain_temp in levels(file_temp$chain)){
    
    message(paste("Analysing chain", chain_temp, "..."))
    
    ## Analyse only spesific chain at once
    file_temp2              <- filter(file_temp, chain == chain_temp)
    file_temp2$raw_clonotype_id <- droplevels(file_temp2$raw_clonotype_id)
    
    ## Go through every clonotype
    vdj <- lapply(levels(file_temp2$raw_clonotype_id), get_vdj, chain = chain_temp, file_temp = file_temp2)
    vdj <- do.call(rbind, vdj)
    tcr_list[[i]] <- vdj
    i <- i + 1
    
  }
  
  ## Combine TRA & TRB or TRD & TRG; add NAs if only TRA or TRB is found
  
  tra_ind   = NA
  trb_ind   = NA
  trg_ind   = NA
  trd_ind   = NA
  tot_tcrab = NA
  tot_tcrgd = NA
  
  for(i in 1:length(tcr_list)){
    
    if(tcr_list[[i]]$chain == "TRA"){tra_ind = i}
    if(tcr_list[[i]]$chain == "TRB"){trb_ind = i}
    if(tcr_list[[i]]$chain == "TRG"){trd_ind = i}
    if(tcr_list[[i]]$chain == "TRD"){trg_ind = i}
    
  }
  
  if(!is.na(tra_ind) & !is.na(trb_ind)){
    tot_tcrab <- merge(tcr_list[[tra_ind]], tcr_list[[trb_ind]], by = "raw_clonotype_id", all = T)
  }
  
  if(!is.na(trg_ind) & !is.na(trd_ind)){
    tot_tcrgd <- merge(tcr_list[[trg_ind]], tcr_list[[trd_ind]], by = "raw_clonotype_id", all = T)
  }
  
  tot_tcr <- rbind(tot_tcrab, tot_tcrgd)
  
  return(tot_tcr)
  
}



splitCDR3nt <- function(cdr3_nt){
  
  require(seqinr)
  
  tra_cdr3_nt = ""
  trb_cdr3_nt = ""
  
  # print(as.character(cdr3_nt))
  temp_nt <- s2c(as.character(cdr3_nt))
  
  if(sum(temp_nt == ";") > 0){
    
    tra_cdr3_nt <- c2s(temp_nt[6:which(temp_nt == ";") - 1])
    trb_cdr3_nt <- c2s(temp_nt[c(which(temp_nt == ";")+5):length(temp_nt)])
    
  }
  
  else{
    
    name = c2s(temp_nt[1:3])
    
    if(name == "TRA"){
      tra_cdr3_nt = c2s(temp_nt[6:length(temp_nt)])
      
    }
    
    if(name == "TRB"){
      trb_cdr3_nt = c2s(temp_nt[5:length(temp_nt)])
      
    }
    
  }
  
  total_cdr3 <- data.frame(cdr3_nt, tra_cdr3_nt, trb_cdr3_nt)
  return(total_cdr3)
  
}

splitCDR3aa <- function(cdr3_aa){
  
  
  # cdr3_aa = as.character(s1a1_tcr_clonotype$cdr3_aa[996])
  
  require(seqinr)
  tra_cdr3_aa = ""
  trb_cdr3_aa = ""
  
  # priaa(as.character(cdr3_aa))
  temp_aa <- s2c(as.character(cdr3_aa))
  
  if(sum(temp_aa == ";") == 1){
    
    tra_cdr3_aa <- c2s(temp_aa[6:which(temp_aa == ";") - 1])
    trb_cdr3_aa <- c2s(temp_aa[c(which(temp_aa == ";")+5):length(temp_aa)])
    
  }
  
  if(sum(temp_aa == ";") == 0){
    
    name = c2s(temp_aa[1:3])
    
    if(name == "TRA"){
      tra_cdr3_aa = c2s(temp_aa[6:length(temp_aa)])
    }
    
    if(name == "TRB"){
      trb_cdr3_aa = c2s(temp_aa[5:length(temp_aa)])
    }
    
  }
  
  total_cdr3 <- data.frame(cdr3_aa, tra_cdr3_aa, trb_cdr3_aa)
  return(total_cdr3)
  
}







## Add new clonotype_id info, because the clonotypes are sample spesific from 10X

generateNewClonotypes <- function(old_clonotypes_df){
  
  temp <- data.frame(table(old_clonotypes_df$cdr3_nt))
  temp <- temp[order(temp$Freq, decreasing = T), ]
  temp$clonotype_new <- paste0("clonotype", 1:nrow(temp))
  
  colnames(temp) <- c("cdr3_nt", "freq", "new_clonotypes_id")
  return(temp)
  
}

newClonotype_df <- function(old_clonotypes_df){
  
  # basically just merge but faster (and neater)
  
  # @ params: old clonotypes df
  # @ output: df with new clonotypes, merged by the same cdr3s_net
  
  
  old_clonotypes_df$new_clonotypes_id <- "NA"
  new_clonotypes_df <- generateNewClonotypes(old_clonotypes_df)
  
  for(cdr3 in levels(old_clonotypes_df$cdr3_nt)){
    
    old_clonotypes_df[old_clonotypes_df$cdr3_nt %in% cdr3, "new_clonotypes_id"] <- new_clonotypes_df[new_clonotypes_df$cdr3_nt %in% cdr3, "new_clonotypes_id"]
    
  }
  
  return(old_clonotypes_df)
  
}

newClonotype_trb_df <- function(old_clonotypes_df){
  
  # basically just merge but faster (and neater)
  
  # @ params: old clonotypes df
  # @ output: df with new clonotypes, merged by the same cdr3s_net
  
  
  
  old_clonotypes_df$new_clonotypes_id <- "NA"
  new_clonotypes_df <- generateNewClonotypes(old_clonotypes_df)
  
  for(cdr3 in levels(old_clonotypes_df$cdr3_trb)){
    
    old_clonotypes_df[old_clonotypes_df$cdr3_trb %in% cdr3, "new_clonotypes_id"] <- new_clonotypes_df[new_clonotypes_df$cdr3_trb %in% cdr3, "new_clonotypes_id"]
    
  }
  
  return(old_clonotypes_df)
  
}


## Add conservative cysteine into every cdr3 determined
add_C <- function(cdr3){
  
  cdr3 <- as.character(cdr3)
  if(nchar(cdr3) > 0){
    
    if(substr(cdr3, 1, 1) != "C"){
      cdr3 <- paste0("C", cdr3)
    }
  }
  
  return(cdr3)
  
}



preprocess_10X_TCR <- function(contig_file, clonotype_file, prefix){
  
  ###################
  # This is the workhorse function
  ###################
  
  
  # 1) Filter and merge the datas to gether
  contig_file                     <- filterContigFiles(contig_file)
  contig_file                     <- contig_file[,c("barcode", "raw_clonotype_id", "chain", "v_gene", "d_gene", "j_gene", "c_gene")]
  
  master_tcr_raw                  <- inner_join(clonotype_file, contig_file, by.x = "raw_clonotype_id", by.y = "raw_clonotype_id")
  master_tcr_raw                  <- master_tcr_raw[!duplicated(master_tcr_raw), ]
  
  
  # 2) Get the consensus TCRab and TCRgd files
  master_tcr_raw_detailed             <- getVDJforFile(file_temp = master_tcr_raw)
  colnames(master_tcr_raw_detailed)   <- gsub("\\.x", "_tra", colnames(master_tcr_raw_detailed))
  colnames(master_tcr_raw_detailed)   <- gsub("\\.y", "_trb", colnames(master_tcr_raw_detailed))
  master_tcr_raw_detailed$tcr_type    <- ifelse(master_tcr_raw_detailed$chain_tra %in% c("TRA", "TRB"), "TCRab", "TCRgd")
  
  # 3) Combine the data
  master_tcr_barcoded             <- inner_join(master_tcr_raw, master_tcr_raw_detailed, by.x = "raw_clonotype_id", by.y = "raw_clonotype_id")
  master_tcr_barcoded             <- dplyr::select(master_tcr_barcoded, -chain, -v_gene, -d_gene, -j_gene, -c_gene)
  
  # 4) Split the CDR3 information into pieces
  nt <- lapply(master_tcr_barcoded$cdr3_nt, FUN = splitCDR3nt); nt <- do.call(rbind, nt)
  aa <- lapply(master_tcr_barcoded$cdr3_aa, FUN = splitCDR3aa); aa <- do.call(rbind, aa)
  
  master_tcr_barcoded <- data.frame(nt, aa, dplyr::select(master_tcr_barcoded, -cdr3_nt_tra, -cdr3_nt_trb, -cdr3_aa_tra, -cdr3_aa_trb))
  
  # 5) Create file just by clonotypes
  master_tcr_clonotype            <- master_tcr_barcoded %>% dplyr::select(-raw_clonotype_id, -barcode)
  master_tcr_clonotype            <- master_tcr_clonotype[!duplicated(master_tcr_clonotype), ]
  master_tcr_clonotype$proportion <- master_tcr_clonotype$reads / sum(master_tcr_clonotype$reads)
  
  # 6) Remove duplicate entries by barcode
  master_tcr_barcoded <- master_tcr_barcoded[!duplicated(master_tcr_barcoded$barcode), ]
  
  ## Add conservative C if it is missing
  master_tcr_barcoded$tra_cdr3_aa <- lapply(master_tcr_barcoded$tra_cdr3_aa, add_C)
  master_tcr_barcoded$trb_cdr3_aa <- lapply(master_tcr_barcoded$trb_cdr3_aa, add_C)
  
  master_tcr_clonotype$tra_cdr3_aa <- lapply(master_tcr_clonotype$tra_cdr3_aa, add_C)
  master_tcr_clonotype$trb_cdr3_aa <- lapply(master_tcr_clonotype$trb_cdr3_aa, add_C)
  
  
  ### ==== Write down ====
  data.table::fwrite(master_tcr_barcoded,  paste0(prefix, "_barcoded.txt"),   sep = "\t", row.names = F)
  data.table::fwrite(master_tcr_clonotype, paste0(prefix, "_clonotyped.txt"), sep = "\t", row.names = F)
  
}

