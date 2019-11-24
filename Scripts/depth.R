## check read depth per exon -----
average_depths <- function(depth_table) {
  colnames(depth_table) <- c('chr', 'start', 'stop', 'exon', 'depth', 'bases_at_depth', 'exon_size', 'percent_at_depth')
  depth_list <- split(depth_table, depth_table$exon)
  depth_list_depths <- rep(NA, length(depth_list))
  for (i in 1:length(depth_list_depths)) {
    exon_depths <- depth_list[[i]]
    exon_depths$totals <- exon_depths$depth*exon_depths$bases_at_depth
    total_depth <- sum(exon_depths$totals)
    depth_list_depths[i] <- total_depth/exon_depths$exon_size[1]
  }
  exon_names <- names(depth_list)
  names(depth_list_depths) <- exon_names
  return(depth_list_depths)
}


## patient 9 ----
pat_9_buffy_depth <- read.delim('G:/depth_files/pat_9_buffy_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_9_buffy_depths <- average_depths(pat_9_buffy_depth)

pat_9_ln_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_ln_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_9_ln_depths <- average_depths(pat_9_ln_depth)

pat_9_oment_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_oment_TST170_all_depth.txt', header = TRUE, stringsAsFactors = TRUE, sep = '\t')
pat_9_oment_depths <- average_depths(pat_9_oment_depth)

pat_9_ovary_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_ovary_TST170_all_depth.txt', header = TRUE, stringsAsFactors = TRUE, sep = '\t')
pat_9_ovary_depths <- average_depths(pat_9_ovary_depth)

pat_9_plasma_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_plasma_TST170_all_depth.txt', header = TRUE, stringsAsFactors = TRUE, sep = '\t')
pat_9_plasma_depths <- average_depths(pat_9_plasma_depth)

## look at specific genes not found in these samples ----
# TP53
pat_9_buffy_TP53 <- pat_9_buffy_depths[grep('TP53', names(pat_9_buffy_depths))]
pat_9_ln_TP53 <- pat_9_ln_depths[grep('TP53', names(pat_9_ln_depths))]
pat_9_oment_TP53 <- pat_9_oment_depths[grep('TP53', names(pat_9_oment_depths))]
pat_9_ovary_TP53 <- pat_9_ovary_depths[grep('TP53', names(pat_9_ovary_depths))]
pat_9_plasma_TP53 <- pat_9_plasma_depths[grep('TP53', names(pat_9_plasma_depths))]

pat_9_TP53 <- data.frame(pat_9_buffy_TP53, pat_9_ln_TP53, pat_9_oment_TP53, pat_9_ovary_TP53, pat_9_plasma_TP53)

# BRCA1
pat_9_buffy_BRCA1 <- pat_9_buffy_depths[grep('BRCA1', names(pat_9_buffy_depths))]
pat_9_ln_BRCA1 <- pat_9_ln_depths[grep('BRCA1', names(pat_9_ln_depths))]
pat_9_oment_BRCA1 <- pat_9_oment_depths[grep('BRCA1', names(pat_9_oment_depths))]
pat_9_ovary_BRCA1 <- pat_9_ovary_depths[grep('BRCA1', names(pat_9_ovary_depths))]
pat_9_plasma_BRCA1 <- pat_9_plasma_depths[grep('BRCA1', names(pat_9_plasma_depths))]

pat_9_BRCA1 <- data.frame(pat_9_buffy_BRCA1, pat_9_ln_BRCA1, pat_9_oment_BRCA1, pat_9_ovary_BRCA1, pat_9_plasma_BRCA1)
