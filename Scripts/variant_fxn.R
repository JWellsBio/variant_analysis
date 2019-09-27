## functions needed in variant processing
'%!in%' <- function(x,y)!('%in%'(x,y))
## function for trimming mutect calls to trusight genes ----

mutect_process <- function(mutect_calls) {
  #ann_idx <- mutect_short$ANN == '' # index of mutations without annotations
  #mutect_shorter <- mutect_short[!ann_idx, ] # remove any without annotation
  
  mutect_info <- mutect_calls$ANN # just the annotation info from snpEff
  info_split <- strsplit(mutect_info, '\\|') # split them by | and now a list
  
  mutant_allele <- rep(NA, length(info_split))
  for (i in 1:length(mutant_allele)) {
    mutant_allele[i] <- info_split[[i]][1]
  }
  mutect_calls$mutant_allele <- mutant_allele
  
  effect <- rep(NA, length(info_split))
  for (i in 1:length(effect)) {
    effect[i] <- info_split[[i]][2]
  }
  mutect_calls$effect <- effect
  
  impact <- rep(NA, length(info_split))
  for (i in 1:length(impact)) {
    impact[i] <- info_split[[i]][3]
  }
  mutect_calls$impact <- impact
  
  gene_name <- rep(NA, length(info_split))
  for (i in 1:length(gene_name)) {
    gene_name[i] <- info_split[[i]][4]
  }
  mutect_calls$gene_name <- gene_name
  
  gene_id <- rep(NA, length(info_split))
  for (i in 1:length(gene_id)) {
    gene_id[i] <- info_split[[i]][5]
  }
  mutect_calls$gene_id <- gene_id
  
  feature_type <- rep(NA, length(info_split))
  for (i in 1:length(feature_type)) {
    feature_type[i] <- info_split[[i]][6]
  }
  mutect_calls$feature_type <- feature_type
  
  feature_id <- rep(NA, length(info_split))
  for (i in 1:length(feature_id)) {
    feature_id[i] <- info_split[[i]][7]
  }
  mutect_calls$feature_id <- feature_id
  
  transcript_biotype <- rep(NA, length(info_split))
  for (i in 1:length(transcript_biotype)) {
    transcript_biotype[i] <- info_split[[i]][8]
  }
  mutect_calls$transcript_biotype <- transcript_biotype
  
  rank_total <- rep(NA, length(info_split))
  for (i in 1:length(rank_total)) {
    rank_total[i] <- info_split[[i]][9]
  }
  mutect_calls$rank_total <- rank_total
  
  hgvs_c <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_c)) {
    hgvs_c[i] <- info_split[[i]][10]
  }
  mutect_calls$hgvs_c <- hgvs_c
  
  hgvs_p <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_p)) {
    hgvs_p[i] <- info_split[[i]][11]
  }
  mutect_calls$hgvs_p <- hgvs_p
  
  cdna_pos <- rep(NA, length(info_split))
  for (i in 1:length(cdna_pos)) {
    cdna_pos[i] <- info_split[[i]][12]
  }
  mutect_calls$cdna_pos <- cdna_pos
  
  cds_pos_cds_len <- rep(NA, length(info_split))
  for (i in 1:length(cds_pos_cds_len)) {
    cds_pos_cds_len[i] <- info_split[[i]][13]
  }
  mutect_calls$cds_pos_cds_len <- cds_pos_cds_len
  
  protein_pos <- rep(NA, length(info_split))
  for (i in 1:length(protein_pos)) {
    protein_pos[i] <- info_split[[i]][14]
  }
  mutect_calls$protein_pos <- protein_pos
  
  dist_to_feat <- rep(NA, length(info_split))
  for (i in 1:length(dist_to_feat)) {
    dist_to_feat[i] <- info_split[[i]][15]
  }
  mutect_calls$dist_to_feat <- dist_to_feat
  
  #all calls with trusight genes and add location and remove ANN column
  mutect_calls$X.CHROM <- gsub('chr', '', mutect_calls$X.CHROM)
  mutect_calls$location <- paste0(mutect_calls$X.CHROM, ':', mutect_calls$POS)
  mutect_calls <- mutect_calls[, -which(names(mutect_calls) %in% c('ANN', 'FILTER', 'gene_id', 'feature_id', 'cdna_pos'))]
  mutect_calls <- mutect_calls[mutect_calls$effect != 'intron_variant', ]
  mutect_calls <- mutect_calls[mutect_calls$effect != 'intergenic_region', ]
  mutect_calls$AF <- as.numeric(mutect_calls$AF)
  mutect_calls <- mutect_calls[!is.na(mutect_calls$AF), ]
  mutect_calls <- mutect_calls[nchar(mutect_calls$REF) == 1, ]
  mutect_calls <- mutect_calls[nchar(mutect_calls$ALT) == 1, ]
  #mutect_trusight <- mutect_trusight[mutect_trusight$hgvs_p != '', ]
  #mutect_trusight <- mutect_trusight[mutect_trusight$effect != 'synonymous_variant', ]
  return(mutect_calls)
}