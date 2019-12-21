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
mean(pat_9_buffy_depths, na.rm = TRUE) #2224.41

pat_9_ln_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_ln_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_9_ln_depths <- average_depths(pat_9_ln_depth)
mean(pat_9_ln_depths, na.rm = TRUE) #433.87328

pat_9_oment_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_oment_TST170_all_depth.txt', header = TRUE, stringsAsFactors = TRUE, sep = '\t')
pat_9_oment_depths <- average_depths(pat_9_oment_depth)
mean(pat_9_oment_depths, na.rm = TRUE) #735.9079

pat_9_ovary_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_ovary_TST170_all_depth.txt', header = TRUE, stringsAsFactors = TRUE, sep = '\t')
pat_9_ovary_depths <- average_depths(pat_9_ovary_depth)
mean(pat_9_ovary_depths, na.rm = TRUE) #825.1706

pat_9_plasma_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_9_plasma_TST170_all_depth.txt', header = TRUE, stringsAsFactors = TRUE, sep = '\t')
pat_9_plasma_depths <- average_depths(pat_9_plasma_depth)
mean(pat_9_plasma_depths, na.rm = TRUE) #1142.788

# look at specific genes not found in these samples ----
#TP53
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

## pat ema ----
pat_ema_heart_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_heart_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_heart_depths <- average_depths(pat_ema_heart_depth)
mean(pat_ema_heart_depths, na.rm = TRUE) #1875.578

pat_ema_l_kidney_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_l_kidney_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_l_kidney_depths <- average_depths(pat_ema_l_kidney_depth)
mean(pat_ema_l_kidney_depths, na.rm = TRUE) #1836.295

pat_ema_r_kidney_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_r_kidney_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_r_kidney_depths <- average_depths(pat_ema_r_kidney_depth)
mean(pat_ema_r_kidney_depths, na.rm = TRUE) #1737.077

pat_ema_liver_1_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_liver_1_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_liver_1_depths <- average_depths(pat_ema_liver_1_depth)
mean(pat_ema_liver_1_depths, na.rm = TRUE) #1351.407

pat_ema_liver_2_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_liver_2_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_liver_2_depths <- average_depths(pat_ema_liver_2_depth)
mean(pat_ema_liver_2_depths, na.rm = TRUE) #1557.616

pat_ema_oment_1_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_oment_1_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_oment_1_depths <- average_depths(pat_ema_oment_1_depth)
mean(pat_ema_oment_1_depths, na.rm = TRUE) #2220.201

pat_ema_oment_2_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_oment_2_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_oment_2_depths <- average_depths(pat_ema_oment_2_depth)
mean(pat_ema_oment_2_depths, na.rm = TRUE) #1858.725

pat_ema_plasma_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_ema_plasma_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_ema_plasma_depths <- average_depths(pat_ema_plasma_depth)
mean(pat_ema_plasma_depths, na.rm = TRUE) #908.9195

## pat 10 ----
pat_10_buffy_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_10_buffy_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_10_buffy_depths <- average_depths(pat_10_buffy_depth)
mean(pat_10_buffy_depths, na.rm = TRUE) #1506.219

pat_10_liver_1_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_10_liver_1_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_10_liver_1_depths <- average_depths(pat_10_liver_1_depth)
mean(pat_10_liver_1_depths, na.rm = TRUE) #2086.211

pat_10_liver_2a_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_10_liver_2a_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_10_liver_2a_depths <- average_depths(pat_10_liver_2a_depth)
mean(pat_10_liver_2a_depths, na.rm = TRUE) #2040.209

pat_10_liver_5_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_10_liver_5_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_10_liver_5_depths <- average_depths(pat_10_liver_5_depth)
mean(pat_10_liver_5_depths, na.rm = TRUE) #936.972

pat_10_plasma_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_10_plasma_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_10_plasma_depths <- average_depths(pat_10_plasma_depth)
mean(pat_10_plasma_depths, na.rm = TRUE) #1150.357

## pat 8 ----
pat_8_buffy_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_8_buffy_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_8_buffy_depths <- average_depths(pat_8_buffy_depth)
mean(pat_8_buffy_depths, na.rm = TRUE) #1228.899

pat_8_axillary_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_8_axillary_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_8_axillary_depths <- average_depths(pat_8_axillary_depth)
mean(pat_8_axillary_depths, na.rm = TRUE) #2410.584

pat_8_breast_1_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_8_breast_1_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_8_breast_1_depths <- average_depths(pat_8_breast_1_depth)
mean(pat_8_breast_1_depths, na.rm = TRUE) #739.5229

pat_8_breast_2_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_8_breast_2_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_8_breast_2_depths <- average_depths(pat_8_breast_2_depth)
mean(pat_8_breast_2_depths, na.rm = TRUE) #274.7677

pat_8_plasma_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_8_plasma_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_8_plasma_depths <- average_depths(pat_8_plasma_depth)
mean(pat_8_plasma_depths, na.rm = TRUE) #129.4204

## pat 2 ----
pat_2_liver_norm_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_2_liver_norm_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_2_liver_norm_depths <- average_depths(pat_2_liver_norm_depth)
mean(pat_2_liver_norm_depths, na.rm = TRUE) #763.6067

pat_2_liver_1_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_2_liver_1_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_2_liver_1_depths <- average_depths(pat_2_liver_1_depth)
mean(pat_2_liver_1_depths, na.rm = TRUE) #1747.352

pat_2_liver_2_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_2_liver_2_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_2_liver_2_depths <- average_depths(pat_2_liver_2_depth)
mean(pat_2_liver_2_depths, na.rm = TRUE) #1767.042

pat_2_breast_1_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_2_breast_1_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_2_breast_1_depths <- average_depths(pat_2_breast_1_depth)
mean(pat_2_breast_1_depths, na.rm = TRUE) #1627.685

pat_2_breast_2_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_2_breast_2_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_2_breast_2_depths <- average_depths(pat_2_breast_2_depth)
mean(pat_2_breast_2_depths, na.rm = TRUE) #1907.079

pat_2_plasma_depth <- read.delim('/Volumes/T5_2TB/depth_files/pat_2_plasma_TST170_all_depth.txt', header = FALSE, stringsAsFactors = TRUE, sep = '\t')
pat_2_plasma_depths <- average_depths(pat_2_plasma_depth)
mean(pat_2_plasma_depths, na.rm = TRUE) #263.342


plasma_means <- c(mean(pat_9_plasma_depths, na.rm = TRUE), mean(pat_ema_plasma_depths, na.rm = TRUE), 
                  mean(pat_10_plasma_depths, na.rm = TRUE), mean(pat_8_plasma_depths, na.rm = TRUE), 
                  mean(pat_2_plasma_depths, na.rm = TRUE))
mean(plasma_means) #718.9495

buffy_means <- c(mean(pat_9_buffy_depths, na.rm = TRUE), mean(pat_10_buffy_depths, na.rm = TRUE), 
                 mean(pat_8_buffy_depths, na.rm = TRUE), mean(pat_2_liver_norm_depths, na.rm = TRUE))
mean(buffy_means)

tumor_means <- c(mean(pat_9_ln_depths, na.rm = TRUE), mean(pat_9_oment_depths, na.rm = TRUE), 
                 mean(pat_9_ovary_depths, na.rm = TRUE), mean(pat_ema_heart_depths, na.rm = TRUE), 
                 mean(pat_ema_l_kidney_depths, na.rm = TRUE), mean(pat_ema_r_kidney_depths, na.rm = TRUE), 
                 mean(pat_ema_liver_1_depths, na.rm = TRUE), mean(pat_ema_liver_2_depths, na.rm = TRUE), 
                 mean(pat_ema_oment_1_depths, na.rm = TRUE), mean(pat_ema_oment_2_depths, na.rm = TRUE), 
                 mean(pat_10_liver_1_depths, na.rm = TRUE), mean(pat_10_liver_2a_depths, na.rm = TRUE), 
                 mean(pat_10_liver_5_depths, na.rm = TRUE), mean(pat_8_axillary_depths, na.rm = TRUE), 
                 mean(pat_8_breast_1_depths, na.rm = TRUE), mean(pat_8_breast_2_depths, na.rm = TRUE), 
                 mean(pat_2_liver_1_depths, na.rm = TRUE), mean(pat_2_liver_2_depths, na.rm = TRUE), 
                 mean(pat_2_breast_1_depths, na.rm = TRUE), mean(pat_2_breast_2_depths, na.rm = TRUE))
tumor_means[1] <- 433.9
mean(tumor_means)
summary(tumor_means)
