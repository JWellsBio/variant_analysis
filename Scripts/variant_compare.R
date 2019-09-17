# clean mutation analyses

library(readxl)
library(ggsci)
library(UpSetR)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(gplots)
library(plotrix)
library(stringr)

'%!in%' <- function(x,y)!('%in%'(x,y))
## function for trimming mutect calls to trusight genes ----

trusight_from_mutect <- function(mutect_short) {
  ann_idx <- mutect_short$ANN == '' # index of mutations without annotations
  mutect_shorter <- mutect_short[!ann_idx, ] # remove any without annotation
  
  mutect_info <- mutect_shorter$ANN # just the annotation info from snpEff
  info_split <- strsplit(mutect_info, '\\|') # split them by | and now a list
  
  mutant_allele <- rep(NA, length(info_split))
  for (i in 1:length(mutant_allele)) {
    mutant_allele[i] <- info_split[[i]][1]
  }
  mutect_shorter$mutant_allele <- mutant_allele
  
  effect <- rep(NA, length(info_split))
  for (i in 1:length(effect)) {
    effect[i] <- info_split[[i]][2]
  }
  mutect_shorter$effect <- effect
  
  impact <- rep(NA, length(info_split))
  for (i in 1:length(impact)) {
    impact[i] <- info_split[[i]][3]
  }
  mutect_shorter$impact <- impact
  
  gene_name <- rep(NA, length(info_split))
  for (i in 1:length(gene_name)) {
    gene_name[i] <- info_split[[i]][4]
  }
  mutect_shorter$gene_name <- gene_name
  
  gene_id <- rep(NA, length(info_split))
  for (i in 1:length(gene_id)) {
    gene_id[i] <- info_split[[i]][5]
  }
  mutect_shorter$gene_id <- gene_id
  
  feature_type <- rep(NA, length(info_split))
  for (i in 1:length(feature_type)) {
    feature_type[i] <- info_split[[i]][6]
  }
  mutect_shorter$feature_type <- feature_type
  
  feature_id <- rep(NA, length(info_split))
  for (i in 1:length(feature_id)) {
    feature_id[i] <- info_split[[i]][7]
  }
  mutect_shorter$feature_id <- feature_id
  
  transcript_biotype <- rep(NA, length(info_split))
  for (i in 1:length(transcript_biotype)) {
    transcript_biotype[i] <- info_split[[i]][8]
  }
  mutect_shorter$transcript_biotype <- transcript_biotype
  
  rank_total <- rep(NA, length(info_split))
  for (i in 1:length(rank_total)) {
    rank_total[i] <- info_split[[i]][9]
  }
  mutect_shorter$rank_total <- rank_total
  
  hgvs_c <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_c)) {
    hgvs_c[i] <- info_split[[i]][10]
  }
  mutect_shorter$hgvs_c <- hgvs_c
  
  hgvs_p <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_p)) {
    hgvs_p[i] <- info_split[[i]][11]
  }
  mutect_shorter$hgvs_p <- hgvs_p
  
  cdna_pos <- rep(NA, length(info_split))
  for (i in 1:length(cdna_pos)) {
    cdna_pos[i] <- info_split[[i]][12]
  }
  mutect_shorter$cdna_pos <- cdna_pos
  
  cds_pos_cds_len <- rep(NA, length(info_split))
  for (i in 1:length(cds_pos_cds_len)) {
    cds_pos_cds_len[i] <- info_split[[i]][13]
  }
  mutect_shorter$cds_pos_cds_len <- cds_pos_cds_len
  
  protein_pos <- rep(NA, length(info_split))
  for (i in 1:length(protein_pos)) {
    protein_pos[i] <- info_split[[i]][14]
  }
  mutect_shorter$protein_pos <- protein_pos
  
  dist_to_feat <- rep(NA, length(info_split))
  for (i in 1:length(dist_to_feat)) {
    dist_to_feat[i] <- info_split[[i]][15]
  }
  mutect_shorter$dist_to_feat <- dist_to_feat
  #how many genes in trusight panel?
  trusight <- intersect(mutect_shorter$gene_name, trusight_genes)
  
  #all calls with trusight genes and add location and remove ANN column
  mutect_trusight <- mutect_shorter[mutect_shorter$gene_name %in% trusight, ]
  mutect_trusight$CHROM <- gsub('chr', '', mutect_trusight$CHROM)
  mutect_trusight$location <- paste0(mutect_trusight$CHROM, ':', mutect_trusight$POS)
  mutect_trusight <- mutect_trusight[, -which(names(mutect_trusight) %in% c('ANN', 'gene_id', 'feature_id', 'cdna_pos'))]
  mutect_trusight <- mutect_trusight[mutect_trusight$effect != 'intron_variant', ]
  return(mutect_trusight)
}

## trusight 170 panel genes ----
trusight_genes <- c('AKT1', 'BRIP1', 'CREBBP', 'FANCI', 'FGFR2', 'JAK3', 'MSH3', 'PALB2', 'RAD51D', 'TSC1', 'AKT2', 'BTK', 'CSF1R', 'FANCL', 'FGFR3', 'KDR', 'MSH6', 'PDGFRA', 'RAD54L', 'TSC2',
                    'AKT3', 'CARD11', 'CTNNB1', 'FBXW7', 'FGFR4', 'KIT', 'MTOR', 'PDGFRB', 'RB1', 'VHL', 'ALK', 'CCND1', 'DDR2', 'FGF1', 'FLT1', 'KMT2A', 'MLL', 'MUTYH', 'PIK3CA', 'RET', 'XRCC2',
                    'APC', 'CCND2', 'DNMT3A', 'FGF2', 'FLT3', 'KRAS', 'MYC', 'PIK3CB', 'RICTOR', 'AR', 'CCNE1', 'EGFR', 'FGF3', 'FOXL2', 'MAP2K1', 'MYCL1', 'PIK3CD', 'ROS1',
                    'ARID1A', 'CD79A', 'EP300', 'FGF4', 'GEN1', 'MAP2K2', 'MYCN', 'PIK3CG', 'RPS6KB1', 'ATM', 'CD79B', 'ERBB2', 'FGF5', 'GNA11', 'MCL1', 'MYD88', 'PIK3R1', 'SLX4',
                    'ATR', 'CDH1', 'ERBB3', 'FGF6', 'GNAQ', 'MDM2', 'NBN', 'PMS2', 'SMAD4', 'BAP1', 'CDK12', 'ERBB4', 'FGF7', 'GNAS', 'MDM4', 'NF1', 'PPP2R2A', 'SMARCB1',
                    'BARD1', 'CDK4', 'ERCC1', 'FGF8', 'HNF1A', 'MET', 'NOTCH1', 'PTCH1', 'SMO', 'BCL2', 'CDK6', 'ERCC2', 'FGF9', 'HRAS', 'MLH1', 'NOTCH2', 'PTEN', 'SRC',
                    'BCL6', 'CDKN2A', 'ERG', 'FGF10', 'IDH1', 'MLLT3', 'NOTCH3', 'PTPN11', 'STK11', 'BRAF', 'CEBPA', 'ESR1', 'FGF14', 'IDH2', 'MPL', 'NPM1', 'RAD51', 'TERT',
                    'BRCA1', 'CHEK1', 'EZH2', 'FGF23', 'INPP4B', 'MRE11A', 'NRAS', 'RAD51B', 'TET2', 'BRCA2', 'CHEK2', 'FAM175A', 'FGFR1', 'JAK2', 'MSH2', 'NRG1', 'RAD51C', 'TP53')

## PATIENT 10 ----
## liver 1----
# import data
pat_10_liver_1 <- read.delim('../pat_10_liver_1_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_10_liver_1 <- trusight_from_mutect(pat_10_liver_1) #39

## liver 2a----
# import data
pat_10_liver_2a <- read.delim('../pat_10_liver_2a_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_10_liver_2a <- trusight_from_mutect(pat_10_liver_2a) #23

## liver 5----
# import data
pat_10_liver_5 <- read.delim('../pat_10_liver_5_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_10_liver_5 <- trusight_from_mutect(pat_10_liver_5) #39

## plasma----
# import data
pat_10_plasma <- read.delim('../pat_10_plasma_bwamem_bqsr_mutect_ann_hg19.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_10_plasma <- trusight_from_mutect(pat_10_plasma) #1378


## looking at how well plasma detects tumor mutations ----

length(intersect(pat_10_liver_1$location, pat_10_plasma$location)) #21/39
length(intersect(pat_10_liver_2a$location, pat_10_plasma$location)) #17/23
length(intersect(pat_10_liver_5$location, pat_10_plasma$location)) #21/39


#pool mutations from 3 mets
pat_10_met_pool <- unique(c(pat_10_liver_1$location, pat_10_liver_2a$location, pat_10_liver_5$location)) #60 unique mutations


pat_10_plasma_found <- pat_10_plasma[pat_10_plasma$location %in% pat_10_met_pool, ] #28 mutations
pat_10_plasma_found_vars <- (pat_10_plasma_found$location) # in 19 different genes

pat_10_plasma_not_found <- pat_10_plasma[pat_10_plasma$location %!in% pat_10_met_pool, ]
pat_10_plasma_not_found <- pat_10_plasma_not_found[pat_10_plasma_not_found$effect == 'missense_variant', ]
pat_10_plasma_not_found <- pat_10_plasma_not_found[pat_10_plasma_not_found$DP.1 > 953, ]
pat_10_plasma_not_found$hgvs_p[pat_10_plasma_not_found$gene_name == 'MTOR']

pat_10_liver_1_not_found <- pat_10_liver_1[pat_10_liver_1$location %!in% pat_10_plasma_found_vars, ]
pat_10_liver_2a_not_found <- pat_10_liver_2a[pat_10_liver_2a$location %!in% pat_10_plasma_found_vars, ]
pat_10_liver_5_not_found <- pat_10_liver_5[pat_10_liver_5$location %!in% pat_10_plasma_found_vars, ]
Reduce(intersect, list(pat_10_liver_1_not_found$lo, pat_10_liver_2a_not_found$location, pat_10_liver_5_not_found$location))



# subset to pooled plasma and take a look
pat_10_liver_1_pooled <- pat_10_liver_1[pat_10_liver_1$location %in% pat_10_plasma_found_vars, ]

pat_10_liver_2a_pooled <- pat_10_liver_2a[pat_10_liver_2a$location %in% pat_10_plasma_found_vars, ]

pat_10_liver_5_pooled <- pat_10_liver_5[pat_10_liver_5$location %in% pat_10_plasma_found_vars, ]

pat_10_plasma_pooled <- pat_10_plasma[pat_10_plasma$location %in% pat_10_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_10_liver_1_pooled <- pat_10_liver_1_pooled[, c('location', 'AF')]
colnames(pat_10_liver_1_pooled) <- c('location', 'AF_liver_1')
pat_10_liver_1_pooled$AF_liver_1[c(2, 6, 9, 10, 11, 12, 19)] <- c(0.015, 0.207, 0.238, 0.213, 0.319, 0.350, 0.250)
pat_10_liver_1_pooled$AF_liver_1 <- as.numeric(pat_10_liver_1_pooled$AF_liver_1)

pat_10_liver_2a_pooled <- pat_10_liver_2a_pooled[, c('location', 'AF')]
colnames(pat_10_liver_2a_pooled) <- c('location', 'AF_liver_2a')
pat_10_liver_2a_pooled$AF_liver_2a[3] <- 0.219
pat_10_liver_2a_pooled$AF_liver_2a[9] <- 0.352
pat_10_liver_2a_pooled$AF_liver_2a[10] <- 0.340
pat_10_liver_2a_pooled$AF_liver_2a <- as.numeric(pat_10_liver_2a_pooled$AF_liver_2a)

pat_10_liver_5_pooled <- pat_10_liver_5_pooled[, c('location', 'AF')]
colnames(pat_10_liver_5_pooled) <- c('location', 'AF_liver_5')
pat_10_liver_5_pooled$AF_liver_5[6] <- 0.267
pat_10_liver_5_pooled$AF_liver_5[9] <- 0.365
pat_10_liver_5_pooled$AF_liver_5[12] <- 0.231
pat_10_liver_5_pooled$AF_liver_5[13] <- 0.274
pat_10_liver_5_pooled$AF_liver_5 <- as.numeric(pat_10_liver_5_pooled$AF_liver_5)

pat_10_plasma_pooled <- pat_10_plasma_pooled[, c('location', 'AF')]
colnames(pat_10_plasma_pooled) <- c('location', 'AF_plasma')
pat_10_plasma_pooled$AF_plasma[c(3, 4, 9, 11, 12, 13, 15, 16, 17, 19, 20, 21, 25)] <- c(0.087, 0.041, 0.222, 0.049, 0.398, 0.214, 0.288, 0.216, 0.268, 0.299, 0.135, 0.039, 0.218) 
pat_10_plasma_pooled$AF_plasma <- as.numeric(pat_10_plasma_pooled$AF_plasma)

# put them together
pat_10_muts_pooled <- merge(pat_10_liver_1_pooled, pat_10_liver_2a_pooled, by = 'location', all = TRUE)
pat_10_muts_pooled <- merge(pat_10_muts_pooled, pat_10_liver_5_pooled, by = 'location', all = TRUE)
pat_10_muts_pooled <- merge(pat_10_muts_pooled, pat_10_plasma_pooled, by = 'location', all = TRUE)
pat_10_muts_pooled <- pat_10_muts_pooled[order(pat_10_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_10_muts_pooled) <- pat_10_muts_pooled$location
row_muts <- rownames(pat_10_muts_pooled)
pat_10_muts_pooled <- pat_10_muts_pooled[, -1]

pat_10_muts_pooled <- as.data.frame(pat_10_muts_pooled)
rownames(pat_10_muts_pooled) <- row_muts

#plot 

pat_10_3 <- c(pat_10_muts_pooled$AF_liver_1, pat_10_muts_pooled$AF_liver_2a, pat_10_muts_pooled$AF_liver_5)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_10_muts_pooled)[1] * dim(pat_10_muts_pooled)[2])
# plasma reds
cell_cols[85:112] <- color.scale(pat_10_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:84] <- color.scale(pat_10_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 28, byrow = FALSE)
pat_10_pooled_t <- data.frame(t(pat_10_muts_pooled))
pat_10_pooled_t <- pat_10_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]

# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_10_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_10_muts_pooled[, 4], na.rm = TRUE),max(pat_10_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(1,-0.9,11,-0.5,round(c(min(pat_10_muts_pooled[, 4], na.rm = TRUE), max(pat_10_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=2.2, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE),max(pat_10_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(12,-0.9,23,-0.5,round(c(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE), max(pat_10_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=13, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 11.5, cex = 1.1, font = 2)

# add NA legend
color.legend(27, -0.9, 28, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.8, at=24.9, cex = 1.1, font = 2)
legend(x=27.1,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_10_muts_pooled)
mut_col_labels[1:5] <- c("FGF9\n5\' UTR", 'EP300\nmissense', "MAP2K1\n5\' UTR", 'SLX4\nmissense', 'SLX4\nmissense')
mut_col_labels[6:10] <- c('ROS1\nsplice region', "FGF4\n5\' UTR", 'PIK3CB\nsplice region', "CCND1\n3\' UTR", 'AR\nsynonymous')
mut_col_labels[11:15] <- c('BRAF\nsplice region', "FLT1\n3\' UTR", 'ATM\nsplice region', 'ATM\nsplice region', 'EP300\nsplice region')
mut_col_labels[16:20] <- c('GNAQ\nsplice region', 'EP300\nsplice region', 'ATM\nsplice acceptor', "CCND1\n3\' UTR", 'BRCA2\nsplice region')
mut_col_labels[21:25] <- c("NRAS\n3\' UTR", 'ATR\nmissense', 'CCND1\nmissense', "CCND1\n3\' UTR", "BRAF\n3\' UTR")
mut_col_labels[26:28] <- c('PIK3CA\nsplice region', 'RB1\nsplice region', 'MSH3\nmissense')
# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_10_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver 1', 'Liver 2a', 'Liver 5')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_10_muts_pooled$AF_liver_1)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_10_muts_pooled$AF_liver_1))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_10_muts_pooled$AF_liver_2a)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_10_muts_pooled$AF_liver_2a))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_10_muts_pooled$AF_liver_5)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_10_muts_pooled$AF_liver_5))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))



## PATIENT 9 ----
## lymph node----
pat_9_ln <- read.delim('../pat_9_ln_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_9_ln <- trusight_from_mutect(pat_9_ln) #48

## oment----
# import data
pat_9_oment <- read.delim('../pat_9_oment_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_9_oment <- trusight_from_mutect(pat_9_oment) #47

## ovary----
# import data
pat_9_ovary <- read.delim('../pat_9_ovary_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_9_ovary <- trusight_from_mutect(pat_9_ovary) #34


## plasma----
# import data
pat_9_plasma <- read.delim('../pat_9_plasma_bwamem_bqsr_mutect_ann_hg19.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_9_plasma <- trusight_from_mutect(pat_9_plasma) #230


## looking at how well plasma detects tumor mutations ----

length(intersect(pat_9_ln$location, pat_9_plasma$location)) #8/48
length(intersect(pat_9_oment$location, pat_9_plasma$location)) #20/47
length(intersect(pat_9_ovary$location, pat_9_plasma$location)) #21/34


#pool mutations from 3 mets
pat_9_met_pool <- unique(c(pat_9_ln$location, pat_9_oment$location, pat_9_ovary$location)) #108 unique mutations


pat_9_plasma_found <- pat_9_plasma[pat_9_plasma$location %in% pat_9_met_pool, ] #30 mutations
pat_9_plasma_found_vars <- (pat_9_plasma_found$location) # in 14 different genes
pat_9_ln_not_found <- pat_9_ln[pat_9_ln$location %!in% pat_9_plasma_found_vars, ]
pat_9_oment_not_found <- pat_9_oment[pat_9_oment$location %!in% pat_9_plasma_found_vars, ]
pat_9_ovary_not_found <- pat_9_ovary[pat_9_ovary$location %!in% pat_9_plasma_found_vars, ]

Reduce(intersect, list(pat_9_ln_not_found$location, pat_9_oment$location, pat_9_ovary_not_found$location))

# subset to pooled plasma and take a look
pat_9_ln_pooled <- pat_9_ln[pat_9_ln$location %in% pat_9_plasma_found_vars, ]

pat_9_oment_pooled <- pat_9_oment[pat_9_oment$location %in% pat_9_plasma_found_vars, ]

pat_9_ovary_pooled <- pat_9_ovary[pat_9_ovary$location %in% pat_9_plasma_found_vars, ]

pat_9_plasma_pooled <- pat_9_plasma[pat_9_plasma$location %in% pat_9_plasma_found_vars, ]


# subset these to tumor AF for new figure
pat_9_ln_pooled <- pat_9_ln_pooled[, c('location', 'AF')]
colnames(pat_9_ln_pooled) <- c('location', 'AF_ln')
pat_9_ln_pooled$AF_ln <- as.numeric(pat_9_ln_pooled$AF_ln)

pat_9_oment_pooled <- pat_9_oment_pooled[, c('location', 'AF')]
colnames(pat_9_oment_pooled) <- c('location', 'AF_oment')
pat_9_oment_pooled$AF_oment[11] <- 0.343
pat_9_oment_pooled$AF_oment[13] <- 0.046
pat_9_oment_pooled$AF_oment[16] <- 0.203
pat_9_oment_pooled$AF_oment[19] <- 0.295
pat_9_oment_pooled$AF_oment <- as.numeric(pat_9_oment_pooled$AF_oment)

pat_9_ovary_pooled <- pat_9_ovary_pooled[, c('location', 'AF')]
colnames(pat_9_ovary_pooled) <- c('location', 'AF_ovary')
pat_9_ovary_pooled$AF_ovary[11] <- 0.330
pat_9_ovary_pooled$AF_ovary[21] <- 0.165
pat_9_ovary_pooled$AF_ovary[20] <- 0.209
pat_9_ovary_pooled$AF_ovary[16] <- 0.246
pat_9_ovary_pooled$AF_ovary[15] <- 0.465
pat_9_ovary_pooled$AF_ovary[14] <- 0.143
pat_9_ovary_pooled$AF_ovary <- as.numeric(pat_9_ovary_pooled$AF_ovary)

pat_9_plasma_pooled <- pat_9_plasma_pooled[, c('location', 'AF')]
colnames(pat_9_plasma_pooled) <- c('location', 'AF_plasma')
pat_9_plasma_pooled$AF_plasma[5] <- 0.020 
pat_9_plasma_pooled$AF_plasma[14] <- 0.333
pat_9_plasma_pooled$AF_plasma[17] <- 0.057
pat_9_plasma_pooled$AF_plasma[30] <- 0.177
pat_9_plasma_pooled$AF_plasma[28] <- 0.209
pat_9_plasma_pooled$AF_plasma[19:22] <- c(0.230, 0.332, 0.260, 0.281)
pat_9_plasma_pooled$AF_plasma <- as.numeric(pat_9_plasma_pooled$AF_plasma)

# put them together
pat_9_muts_pooled <- merge(pat_9_ln_pooled, pat_9_oment_pooled, by = 'location', all = TRUE)
pat_9_muts_pooled <- merge(pat_9_muts_pooled, pat_9_ovary_pooled, by = 'location', all = TRUE)
pat_9_muts_pooled <- merge(pat_9_muts_pooled, pat_9_plasma_pooled, by = 'location', all = TRUE)
pat_9_muts_pooled <- pat_9_muts_pooled[order(pat_9_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_9_muts_pooled) <- pat_9_muts_pooled$location
row_muts <- rownames(pat_9_muts_pooled)
pat_9_muts_pooled <- pat_9_muts_pooled[, -1]

pat_9_muts_pooled <- as.data.frame(pat_9_muts_pooled)
rownames(pat_9_muts_pooled) <- row_muts
#plot 

pat_9_3 <- c(pat_9_muts_pooled$AF_ln, pat_9_muts_pooled$AF_oment, pat_9_muts_pooled$AF_ovary)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_9_muts_pooled)[1] * dim(pat_9_muts_pooled)[2])
# plasma reds
cell_cols[91:120] <- color.scale(pat_9_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:90] <- color.scale(pat_9_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 30, byrow = FALSE)
pat_9_pooled_t <- data.frame(t(pat_9_muts_pooled))
pat_9_pooled_t <- pat_9_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]
# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_9_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_9_muts_pooled[, 4], na.rm = TRUE),max(pat_9_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(1,-0.9,11,-0.5,round(c(min(pat_9_muts_pooled[, 4], na.rm = TRUE), max(pat_9_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=2.2, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE),max(pat_9_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(12,-0.9,23,-0.5,round(c(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE), max(pat_9_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=13, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 11.5, cex = 1.1, font = 2)

# add NA legend
color.legend(29, -0.9, 30, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.8, at=26.9, cex = 1.1, font = 2)
legend(x=29.1,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_9_muts_pooled)
mut_col_labels[1:5] <- c('MTOR\nsynonymous', "CCND1\n3\' UTR", 'NOTCH3\nsplice region', 'EGFR\nsynonymous', 'NOTCH1\nsynonymous')
mut_col_labels[6:10] <- c("ALK\n3\' UTR", 'PIK3CB\nsplice region', 'ALK\nsynonymous', 'BRIP1\nmissense', 'EP300\nmissense')
mut_col_labels[11:15] <- c('BRIP1\nsynonymous', 'GNAQ\nsplice region', 'ATM\nsplice region', "MAP2K1\n5\' UTR", 'ATM\nsplice region')
mut_col_labels[16:20] <- c('ATM\nsplice region', "CCND1\n3\' UTR", 'EP300\nsplice region', 'EP300\nsplice region', 'NOTCH1\nmissense')
mut_col_labels[21:25] <- c("CCND1\n3\' UTR", 'BRCA2\nsplice region', 'EGFR\nmissense', 'EGFR\nmissense', 'EGFR\nmissense')
mut_col_labels[26:30] <- c('EGFR\nsynonymous', 'EGFR\nmissense', 'EGFR\nmissense', 'PIK3CA\nsplice region', 'EGFR\nmissense')
# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_9_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Lymph\nMet', 'Omental\nMet', 'Ovary\nMet')
axis(2, at = c(0.6, 1.6, 2.6, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_9_muts_pooled$AF_ln)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_9_muts_pooled$AF_ln))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_9_muts_pooled$AF_oment)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_9_muts_pooled$AF_oment))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_9_muts_pooled$AF_ovary)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_9_muts_pooled$AF_ovary))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))

## ADD FUNCTION HISTOGRAM

##upset
pat_9_lymph_met_muts <- pat_9_ln$location
pat_9_lymph_met_muts <- as.data.frame(pat_9_lymph_met_muts)
pat_9_lymph_met_muts$mutation <- rep(1, nrow(pat_9_lymph_met_muts))
colnames(pat_9_lymph_met_muts) <- c('mutation', 'Lymph_Met')

pat_9_oment_met_muts <- pat_9_oment$location
pat_9_oment_met_muts <- as.data.frame(pat_9_oment_met_muts)
pat_9_oment_met_muts$mutation <- rep(1, nrow(pat_9_oment_met_muts))
colnames(pat_9_oment_met_muts) <- c('mutation', 'Omental_Met')

pat_9_ovary_met_muts <- pat_9_ovary$location
pat_9_ovary_met_muts <- as.data.frame(pat_9_ovary_met_muts)
pat_9_ovary_met_muts$mutation <- rep(1, nrow(pat_9_ovary_met_muts))
colnames(pat_9_ovary_met_muts) <- c('mutation', 'Ovary_Met')

pat_9_plasma_muts <- pat_9_plasma$location
pat_9_plasma_muts <- as.data.frame(pat_9_plasma_muts)
pat_9_plasma_muts$mutation <- rep(1, nrow(pat_9_plasma_muts))
colnames(pat_9_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_9_lymph_met_muts, pat_9_oment_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_ovary_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

rownames(upset_df) <- upset_df$mutation
row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:4])
randomnumber <- round(runif(4, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 6, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(4)
orange_pal <- pal_material('orange', n = 8, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(4)
purple_pal <- pal_material('purple', n= 4, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(2)
bar_colors <- c(blue_pal[1:3], orange_pal[1:4], purple_pal[1], 'darkgreen')

upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray60', Lymph_Met = 'white', Omental_Met = 'white', 
                                                                 Ovary_Met = 'white')))), 
      intersections = list(list('Lymph_Met'), 
                           list('Ovary_Met'), 
                           list('Omental_Met'), 
                           list('Ovary_Met', 'Plasma'), 
                           list('Lymph_Met', 'Plasma'), 
                           list('Omental_Met', 'Plasma'), 
                           list('Omental_Met', 'Ovary_Met'), 
                           list('Omental_Met', 'Ovary_Met', 'Plasma'), 
                           list(colnames(upset_df))),
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray60', 'goldenrod4', 'aquamarine3', 'chocolate3'), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
      text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 2.0))

pat_9_oment_met_muts <- pat_9_oment$location
pat_9_oment_met_muts <- as.data.frame(pat_9_oment_met_muts)
pat_9_oment_met_muts$mutation <- rep(1, nrow(pat_9_oment_met_muts))
colnames(pat_9_oment_met_muts) <- c('mutation', 'Omental_Met')

pat_9_ovary_met_muts <- pat_9_ovary$location
pat_9_ovary_met_muts <- as.data.frame(pat_9_ovary_met_muts)
pat_9_ovary_met_muts$mutation <- rep(1, nrow(pat_9_ovary_met_muts))
colnames(pat_9_ovary_met_muts) <- c('mutation', 'Ovary_Met')

pat_9_lymph_met_muts <- pat_9_ln$location
pat_9_lymph_met_muts <- as.data.frame(pat_9_lymph_met_muts)
pat_9_lymph_met_muts$mutation <- rep(1, nrow(pat_9_lymph_met_muts))
colnames(pat_9_lymph_met_muts) <- c('mutation', 'Lymph_Met')

pat_9_plasma_muts <- pat_9_plasma$location
pat_9_plasma_muts <- as.data.frame(pat_9_plasma_muts)
pat_9_plasma_muts$mutation <- rep(1, nrow(pat_9_plasma_muts))
colnames(pat_9_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_9_oment_met_muts, pat_9_ovary_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_lymph_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
upset_df$mutation <- as.character(upset_df$mutation)

pat_9_ln_only <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_ln_only <- pat_9_ln_only$mutation
pat_9_ln_af <- pat_9_ln[pat_9_ln$location %in% pat_9_ln_only, ]
pat_9_ln_af <- pat_9_ln_af$AF
pat_9_ln_sd <- sd(as.numeric(pat_9_ln_af))
pat_9_ln_af <- median(as.numeric(pat_9_ln_af))

pat_9_om_only <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_om_only <- pat_9_om_only$mutation
pat_9_om_af <- pat_9_oment[pat_9_oment$location %in% pat_9_om_only, ]
pat_9_om_af <- pat_9_om_af$AF
pat_9_om_sd <- sd(as.numeric(pat_9_om_af))
pat_9_om_af <- median(as.numeric(pat_9_om_af))

pat_9_ov_only <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_ov_only <- pat_9_ov_only$mutation
pat_9_ov_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_only, ]
pat_9_ov_af <- pat_9_ov_af$AF
pat_9_ov_sd <- sd(as.numeric(pat_9_ov_af))
pat_9_ov_af <- median(as.numeric(pat_9_ov_af))


pat_9_pl_ov <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov <- pat_9_pl_ov$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ovary_met_af[2] <- 0.246
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_plasma_af[2] <- 0.281
pat_9_pl_ov_af <- c(pat_9_ovary_met_af, pat_9_plasma_af)
pat_9_pl_ov_sd <- sd(as.numeric(pat_9_pl_ov_af))
pat_9_pl_ov_af <- median(as.numeric(pat_9_pl_ov_af))

pat_9_pl_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_pl_ly <- pat_9_pl_ly$mutation
pat_9_lymph_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ly, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ly, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_pl_ly_af <- c(pat_9_lymph_met_af, pat_9_plasma_af)
pat_9_pl_ly_sd <- sd(as.numeric(pat_9_pl_ly_af))
pat_9_pl_ly_af <- median(as.numeric(pat_9_pl_ly_af))


pat_9_pl_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_pl_om <- pat_9_pl_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_oment_met_af[2] <- 0.046
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_plasma_af[2] <- 0.057
pat_9_pl_om_af <- c(pat_9_oment_met_af, pat_9_plasma_af)
pat_9_pl_om_sd <- sd(as.numeric(pat_9_pl_om_af))
pat_9_pl_om_af <- median(as.numeric(pat_9_pl_om_af))

pat_9_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_ov_om <- pat_9_ov_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ov_om_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af)
pat_9_ov_om_sd <- sd(as.numeric(pat_9_ov_om_af))
pat_9_ov_om_af <- median(as.numeric(pat_9_ov_om_af))
#FIX THE REST OF THESE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_9_pl_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov_om <- pat_9_pl_ov_om$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ovary_met_af[9] <- 0.330
pat_9_ovary_met_af[12] <- 0.143
pat_9_ovary_met_af[13] <- 0.465
pat_9_ovary_met_af[16] <- 0.209
pat_9_ovary_met_af[17] <- 0.165
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_plasma_af[c(1, 9, 11, 12, 13, 16, 17)] <- c(0.20, 0.333, 0.230, 0.332, 0.260, 0.209, 0.177)


pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_oment_met_af[c(9, 13, 16)] <- c(0.343, 0.203, 0.295)
pat_9_pl_ov_om_af <- c(pat_9_ovary_met_af, pat_9_plasma_af, pat_9_oment_met_af)
pat_9_pl_ov_om_sd <- sd(as.numeric(pat_9_pl_ov_om_af))
pat_9_pl_ov_om_af <- median(as.numeric(pat_9_pl_ov_om_af))

pat_9_all <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
pat_9_all <- pat_9_all$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_all, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_all, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_all, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_lymph_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_all, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$AF
pat_9_all_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af, pat_9_plasma_af, pat_9_lymph_met_af)
pat_9_all_sd <- sd(as.numeric(pat_9_all_af))
pat_9_all_af <- median(as.numeric(pat_9_all_af))


pat_9_afs <- c(pat_9_ln_af, pat_9_om_af, pat_9_ov_af, pat_9_pl_ly_af, pat_9_pl_ov_af, pat_9_pl_om_af, pat_9_ov_om_af, pat_9_pl_ov_om_af, pat_9_all_af)
pat_9_sds <- c(pat_9_ln_sd, pat_9_om_sd, pat_9_ov_sd, pat_9_pl_ly_sd, pat_9_pl_ov_sd, pat_9_pl_om_sd, pat_9_ov_om_sd, pat_9_pl_ov_om_sd, pat_9_all_sd)
pat_9_labels <- c('Lymph Met', 'Omental Met', 'Ovary Met', 'Lymph Met\n+ Plasma', 'Ovary Met\n+ Plasma', 'Omental Met\n+ Plasma', 'Omental Met\n+ Ovary Met', 
                  'Omental Met\n+ Ovary Met\n+ Plasma', 'All Samples')
pat_9_afs <- data.frame(pat_9_labels, pat_9_afs, pat_9_sds, bar_colors)
p<-ggplot(data=pat_9_afs, aes(x=pat_9_labels, y=pat_9_afs, fill = bar_colors)) +
  geom_bar(stat="identity", width = 0.5, color = 'black') + scale_x_discrete(limits=pat_9_afs$pat_9_labels) + ylab('Median Mutant Allele Frequency') +
  xlab('Intersections') + theme_bw() + 
  theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) + 
  scale_fill_manual("legend", values = bar_colors[c(1:3,8,4:6,7,9)]) +
  geom_errorbar(aes(ymin=pat_9_afs, ymax=pat_9_afs+pat_9_sds), width=.2,
                position=position_dodge(.9)) + theme(axis.text=element_text(size=14, face = 'bold'),
                                                     axis.title=element_text(size=16,face="bold"))

p



# patient 9
pat_9_lymph_met_stats <- pat_9_ln[, c('location', 'AF')]
pat_9_oment_met_stats <- pat_9_oment[, c('location', 'AF')]
pat_9_oment_met_stats$AF[c(21, 24, 32, 43)] <- c(0.343, 0.046, 0.203, 0.295)
pat_9_ovary_met_stats <- pat_9_ovary[, c('location', 'AF')]
pat_9_ovary_met_stats$AF[c(17, 22, 23, 24, 32, 34)] <- c(0.330, 0.143, 0.465, 0.287, 0.209, 0.165)
pat_9_all <- rbind(pat_9_lymph_met_stats, pat_9_oment_met_stats)
pat_9_all <- rbind(pat_9_all, pat_9_ovary_met_stats)
pat_9_all$AF <- as.numeric(pat_9_all$AF)
pat_9_all <- pat_9_all[order(pat_9_all$AF, decreasing = TRUE), ]
pat_9_all$color <- ifelse(pat_9_all$location %in% pat_9_plasma$location, 'red', 'black')
barplot(pat_9_all$AF, col = pat_9_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 9\n(3 tumors)', ylim = c(0,1.0))

## CHECKING HOW TO DRAW EVOL CLADE ----

#how many in all? this is the length of the main branch moving from normal
pat_9_all_common <- Reduce(intersect, list(pat_9_ln_no_introns$location, pat_9_oment_no_introns$location, pat_9_ovary_no_introns$location)) #1

#which ones?
pat_9_all <- pat_9_ln_no_introns[pat_9_ln_no_introns$location %in% pat_9_all_common, ]
#which in plasma?
length(intersect(pat_9_all_common, pat_9_plasma_no_introns$location)) #all 1

#exclude those, these can be the lengths of the next branches if split solo
pat_9_ln_not_all <- pat_9_ln_no_introns[pat_9_ln_no_introns$location %!in% pat_9_all_common, ] #47
pat_9_oment_not_all <- pat_9_oment_no_introns[pat_9_oment_no_introns$location %!in% pat_9_all_common, ] #46
pat_9_ovary_not_all <- pat_9_ovary_no_introns[pat_9_ovary_no_introns$location %!in% pat_9_all_common, ] #33

#gather all locations left after all common
pat_9_2nd_locations <- c(pat_9_ln_not_all$location, pat_9_oment_not_all$location, pat_9_ovary_not_all$location)
pat_9_2nd_locations <- table(pat_9_2nd_locations)
pat_9_2nd_locations <- pat_9_2nd_locations[order(pat_9_2nd_locations, decreasing = TRUE)]

#only keep those that hit 2 out of 3
pat_9_2nd_locations_2 <- pat_9_2nd_locations[pat_9_2nd_locations == 2]
pat_9_2nd_locations_2 <- names(pat_9_2nd_locations_2) #19

#subset to those to see how many are shared at least twice in each
pat_9_ln_2_of_3 <- pat_9_ln_not_all[pat_9_ln_not_all$location %in% pat_9_2nd_locations_2, ] #1/19
pat_9_oment_2_of_3 <- pat_9_oment_not_all[pat_9_oment_not_all$location %in% pat_9_2nd_locations_2, ] #19/19
pat_9_ovary_2_of_3 <- pat_9_ovary_not_all[pat_9_ovary_not_all$location %in% pat_9_2nd_locations_2, ] #18/19
#ln is clearlythe least like the others, branch it off here and continue with om and ov

#ln excluding commons, length of ln solo branch
pat_9_ln_only <- pat_9_ln_not_all[pat_9_ln_not_all$location %!in% pat_9_2nd_locations_2, ] #46
# how many of these in plasma?
length(intersect(pat_9_ln_only$location, pat_9_plasma_no_introns$location)) #7


#other 2 but not ln
pat_9_oment_not_ln <- pat_9_oment_2_of_3[pat_9_oment_2_of_3$location %!in% pat_9_ln_2_of_3$location, ] #18
pat_9_ovary_not_ln <- pat_9_ovary_2_of_3[pat_9_ovary_2_of_3$location %!in% pat_9_ln_2_of_3$location, ] #18

#in common? (should be length of the 2 just made)
length(intersect(pat_9_oment_not_ln$location, pat_9_ovary_not_ln$location)) #all 18
# how many of these in plasma?
length(Reduce(intersect, list(pat_9_oment_not_ln$location, pat_9_ovary_not_ln$location, pat_9_plasma_no_introns$location))) #17

#subset out ln 3
pat_9_oment_not_ln_long <- pat_9_oment_not_all[pat_9_oment_not_all$location %!in% pat_9_ln_not_all$location, ] #45
pat_9_ovary_not_ln_long <- pat_9_ovary_not_all[pat_9_ovary_not_all$location %!in% pat_9_ln_not_all$location, ] #33

#length of common branch
pat_9_ovary_oment <- intersect(pat_9_oment_not_ln_long$location, pat_9_ovary_not_ln_long$location) #18

#subset out those, lengths of solo branches
pat_9_oment_final_length <- pat_9_oment_not_ln_long[pat_9_oment_not_ln_long$location %!in% pat_9_ovary_oment, ] #27
pat_9_ovary_final_length <- pat_9_ovary_not_ln_long[pat_9_ovary_not_ln_long$location %!in% pat_9_ovary_oment, ] #15
# how many of these in plasma?
length(intersect(pat_9_oment_final_length$location, pat_9_plasma_no_introns$location)) #2
length(intersect(pat_9_ovary_final_length$location, pat_9_plasma_no_introns$location)) #3








## PATIENT 8 ----
## axillary----
# import data
pat_8_axillary <- read.delim('../pat_8_axillary_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_8_axillary <- trusight_from_mutect(pat_8_axillary) #29

## breast 1----
# import data
pat_8_breast_1 <- read.delim('../pat_8_breast_1_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_8_breast_1 <- trusight_from_mutect(pat_8_breast_1) #40

## breast 2----
# import data
pat_8_breast_2 <- read.delim('../pat_8_breast_2_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_8_breast_2 <- trusight_from_mutect(pat_8_breast_2) #1367

## plasma----
# import data
pat_8_plasma <- read.delim('../pat_8_plasma_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_8_plasma <- trusight_from_mutect(pat_8_plasma) #54

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_8_axillary$location, pat_8_plasma$location)) #6/29
length(intersect(pat_8_breast_1$location, pat_8_plasma$location)) #5/40
length(intersect(pat_8_breast_2$location, pat_8_plasma$location)) #18/1367



#pool mutations from 3 mets
pat_8_met_pool <- unique(c(pat_8_axillary$location, pat_8_breast_1$location, pat_8_breast_2$location)) #1402 unique mutations


pat_8_plasma_found <- pat_8_plasma[pat_8_plasma$location %in% pat_8_met_pool, ] #20 mutations
pat_8_plasma_found_vars <- (pat_8_plasma_found$location) # in 13 different genes

pat_8_axillary_not_found <- pat_8_axillary[pat_8_axillary$location %!in% pat_8_plasma_found_vars, ]
pat_8_breast_1_not_found <- pat_8_breast_1[pat_8_breast_1$location %!in% pat_8_plasma_found_vars, ]
pat_8_breast_2_not_found <- pat_8_breast_2[pat_8_breast_2$location %!in% pat_8_plasma_found_vars, ]
Reduce(intersect, list(pat_8_axillary_not_found$location, pat_8_breast_1_not_found$location, pat_8_breast_2_not_found$location))

# subset to pooled plasma and take a look
pat_8_axillary_pooled <- pat_8_axillary[pat_8_axillary$location %in% pat_8_plasma_found_vars, ]

pat_8_breast_1_pooled <- pat_8_breast_1[pat_8_breast_1$location %in% pat_8_plasma_found_vars, ]

pat_8_breast_2_pooled <- pat_8_breast_2[pat_8_breast_2$location %in% pat_8_plasma_found_vars, ]

pat_8_plasma_pooled <- pat_8_plasma[pat_8_plasma$location %in% pat_8_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_8_axillary_pooled <- pat_8_axillary_pooled[, c('location', 'AF')]
colnames(pat_8_axillary_pooled) <- c('location', 'AF_axillary')
pat_8_axillary_pooled$AF_axillary[c(2, 3)] <- c(0.292, 0.328)
pat_8_axillary_pooled$AF_axillary <- as.numeric(pat_8_axillary_pooled$AF_axillary)

pat_8_breast_1_pooled <- pat_8_breast_1_pooled[, c('location', 'AF')]
colnames(pat_8_breast_1_pooled) <- c('location', 'AF_breast_1')
pat_8_breast_1_pooled$AF_breast_1[3] <- 0.335
pat_8_breast_1_pooled$AF_breast_1 <- as.numeric(pat_8_breast_1_pooled$AF_breast_1)

pat_8_breast_2_pooled <- pat_8_breast_2_pooled[, c('location', 'AF')]
colnames(pat_8_breast_2_pooled) <- c('location', 'AF_breast_2')
pat_8_breast_2_pooled$AF_breast_2 <- as.numeric(pat_8_breast_2_pooled$AF_breast_2)

pat_8_plasma_pooled <- pat_8_plasma_pooled[, c('location', 'AF')]
colnames(pat_8_plasma_pooled) <- c('location', 'AF_plasma')
pat_8_plasma_pooled$AF_plasma[c(3, 6)] <- c(0.313, 0.328) 
pat_8_plasma_pooled$AF_plasma <- as.numeric(pat_8_plasma_pooled$AF_plasma)

# put them together
pat_8_muts_pooled <- merge(pat_8_axillary_pooled, pat_8_breast_1_pooled, by = 'location', all = TRUE)
pat_8_muts_pooled <- merge(pat_8_muts_pooled, pat_8_breast_2_pooled, by = 'location', all = TRUE)
pat_8_muts_pooled <- merge(pat_8_muts_pooled, pat_8_plasma_pooled, by = 'location', all = TRUE)
pat_8_muts_pooled <- pat_8_muts_pooled[order(pat_8_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_8_muts_pooled) <- pat_8_muts_pooled$location
row_muts <- rownames(pat_8_muts_pooled)
pat_8_muts_pooled <- pat_8_muts_pooled[, -1]

pat_8_muts_pooled <- as.data.frame(pat_8_muts_pooled)
rownames(pat_8_muts_pooled) <- row_muts

#plot 

pat_8_3 <- c(pat_8_muts_pooled$AF_axillary, pat_8_muts_pooled$AF_breast_1, pat_8_muts_pooled$AF_breast_2)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_8_muts_pooled)[1] * dim(pat_8_muts_pooled)[2])
# plasma reds
cell_cols[61:80] <- color.scale(pat_8_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:60] <- color.scale(pat_8_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 20, byrow = FALSE)
pat_8_pooled_t <- data.frame(t(pat_8_muts_pooled))
pat_8_pooled_t <- pat_8_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]

# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_8_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_8_muts_pooled[, 4], na.rm = TRUE),max(pat_8_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(1,-0.9,8,-0.5,round(c(min(pat_8_muts_pooled[, 4], na.rm = TRUE), max(pat_8_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=1.8, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_8_muts_pooled[, 1:3], na.rm = TRUE),max(pat_8_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(9,-0.9,16,-0.5,round(c(min(pat_8_muts_pooled[, 1:3], na.rm = TRUE), max(pat_8_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=9.7, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 8.5, cex = 1.1, font = 2)

# add NA legend
color.legend(19, -0.9, 20, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.8, at=17.7, cex = 1.1, font = 2)
legend(x=19.2,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_8_muts_pooled)
mut_col_labels[1:5] <- c('MTOR\nmissense', 'FGF6\nsynonymous', 'SLX4\nmissense', 'SLX4\nmissense', 'ALK\nmissense')
mut_col_labels[6:10] <- c('SLX4\nmissense', 'NOTCH3\nsplice region', 'MSH3\nmissense', 'MAP2K2\nsynonymous', 'SLX4\nsynonymous')
mut_col_labels[11:15] <- c('ATM\nsplice region', 'EP300\nsplice region', 'BRCA2\nsynonymous', 'ATM\nmissense', 'NOTCH3\nmissense')
mut_col_labels[16:20] <- c('BRIP1\nsynonymous', 'NOTCH3\nsynonymous', 'ATM\nsplice region', 'PIK3CA\nsplice region', 'ROS1\nmissense')

# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_8_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Axillary', 'Breast 1', 'Breast 2')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_8_muts_pooled$AF_axillary)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_8_muts_pooled$AF_axillary))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_8_muts_pooled$AF_breast_1)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_8_muts_pooled$AF_breast_1))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_8_muts_pooled$AF_breast_2)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_8_muts_pooled$AF_breast_2))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))


## PATIENT 2 ----
## breast 1----
# import data
pat_2_breast_1 <- read.delim('../pat_2_breast_1_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_2_breast_1 <- trusight_from_mutect(pat_2_breast_1) #21

## breast 2----
# import data
pat_2_breast_2 <- read.delim('../pat_2_breast_2_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_2_breast_2 <- trusight_from_mutect(pat_2_breast_2) #24

## liver 1----
# import data
pat_2_liver_1 <- read.delim('../pat_2_liver_1_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_2_liver_1 <- trusight_from_mutect(pat_2_liver_1) #28

## liver 2----
# import data
pat_2_liver_2 <- read.delim('../pat_2_liver_2_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_2_liver_2 <- trusight_from_mutect(pat_2_liver_2) #34

## plasma----
pat_2_plasma <- read.delim('../pat_2_plasma_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_2_plasma <- trusight_from_mutect(pat_2_plasma) #65

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_2_breast_1$location, pat_2_plasma$location)) #6/21
length(intersect(pat_2_breast_2$location, pat_2_plasma$location)) #4/24
length(intersect(pat_2_liver_1$location, pat_2_plasma$location)) #7/28
length(intersect(pat_2_liver_2$location, pat_2_plasma$location)) #7/34



#pool mutations from 3 mets
pat_2_met_pool <- unique(c(pat_2_liver_1$location, pat_2_liver_2$location, pat_2_breast_1$location, pat_2_breast_2$location)) #59 unique mutations


pat_2_plasma_found <- pat_2_plasma[pat_2_plasma$location %in% pat_2_met_pool, ] #8 mutations
pat_2_plasma_found_vars <- (pat_2_plasma_found$location) # in 6 different genes
pat_2_plasma_not_found <- pat_2_plasma[pat_2_plasma$location %!in% pat_2_met_pool, ]
pat_2_plasma_not_found <- pat_2_plasma_not_found[pat_2_plasma_not_found$effect != 'synonymous_variant', ]
pat_2_plasma_not_found <- pat_2_plasma_not_found[pat_2_plasma_not_found$effect == 'missense_variant', ]


# subset to pooled plasma and take a look
pat_2_liver_1_pooled <- pat_2_liver_1[pat_2_liver_1$location %in% pat_2_plasma_found_vars, ]

pat_2_liver_2_pooled <- pat_2_liver_2[pat_2_liver_2$location %in% pat_2_plasma_found_vars, ]

pat_2_breast_1_pooled <- pat_2_breast_1[pat_2_breast_1$location %in% pat_2_plasma_found_vars, ]

pat_2_breast_2_pooled <- pat_2_breast_2[pat_2_breast_2$location %in% pat_2_plasma_found_vars, ]

pat_2_plasma_pooled <- pat_2_plasma[pat_2_plasma$location %in% pat_2_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_2_liver_1_pooled <- pat_2_liver_1_pooled[, c('location', 'AF')]
colnames(pat_2_liver_1_pooled) <- c('location', 'AF_liver_1')
pat_2_liver_1_pooled$AF_liver_1[c(4, 5, 7)] <- c(0.237, 0.295, 0.091)
pat_2_liver_1_pooled$AF_liver_1 <- as.numeric(pat_2_liver_1_pooled$AF_liver_1)

pat_2_liver_2_pooled <- pat_2_liver_2_pooled[, c('location', 'AF')]
colnames(pat_2_liver_2_pooled) <- c('location', 'AF_liver_2')
pat_2_liver_2_pooled$AF_liver_2[c(4, 5)] <- c(0.337, 0.297)
pat_2_liver_2_pooled$AF_liver_2 <- as.numeric(pat_2_liver_2_pooled$AF_liver_2)

pat_2_breast_1_pooled <- pat_2_breast_1_pooled[, c('location', 'AF')]
colnames(pat_2_breast_1_pooled) <- c('location', 'AF_breast_1')
pat_2_breast_1_pooled$AF_breast_1[c(4,5,6)] <- c(0.251, 0.289, 0.080)
pat_2_breast_1_pooled$AF_breast_1 <- as.numeric(pat_2_breast_1_pooled$AF_breast_1)

pat_2_breast_2_pooled <- pat_2_breast_2_pooled[, c('location', 'AF')]
colnames(pat_2_breast_2_pooled) <- c('location', 'AF_breast_2')
pat_2_breast_2_pooled$AF_breast_2[4] <- 0.242
pat_2_breast_2_pooled$AF_breast_2 <- as.numeric(pat_2_breast_2_pooled$AF_breast_2)

pat_2_plasma_pooled <- pat_2_plasma_pooled[, c('location', 'AF')]
colnames(pat_2_plasma_pooled) <- c('location', 'AF_plasma')
pat_2_plasma_pooled$AF_plasma[c(8)] <- c(0.197) 
pat_2_plasma_pooled$AF_plasma <- as.numeric(pat_2_plasma_pooled$AF_plasma)

# put them together
pat_2_muts_pooled <- merge(pat_2_liver_1_pooled, pat_2_liver_2_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_breast_1_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_breast_2_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_plasma_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- pat_2_muts_pooled[order(pat_2_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_2_muts_pooled) <- pat_2_muts_pooled$location
row_muts <- rownames(pat_2_muts_pooled)
pat_2_muts_pooled <- pat_2_muts_pooled[, -1]

pat_2_muts_pooled <- as.data.frame(pat_2_muts_pooled)
rownames(pat_2_muts_pooled) <- row_muts

#plot 

pat_2_3 <- c(pat_2_muts_pooled$AF_liver_1, pat_2_muts_pooled$AF_liver_2, pat_2_muts_pooled$AF_breast_1, pat_2_muts_pooled$AF_breast_2)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_2_muts_pooled)[1] * dim(pat_2_muts_pooled)[2])
# plasma reds
cell_cols[33:40] <- color.scale(pat_2_muts_pooled[, 5], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:32] <- color.scale(pat_2_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 8, byrow = FALSE)
pat_2_pooled_t <- data.frame(t(pat_2_muts_pooled))
pat_2_pooled_t <- pat_2_pooled_t[c(5, 1:4), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(5, 1:4), ]

# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_2_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_2_muts_pooled[, 5], na.rm = TRUE),max(pat_2_muts_pooled[, 5], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(0.5,-0.9,2.5,-0.5,round(c(min(pat_2_muts_pooled[, 4], na.rm = TRUE), max(pat_2_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=1.7, at=0.8, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_2_muts_pooled[, 1:4], na.rm = TRUE),max(pat_2_muts_pooled[, 1:4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(3,-0.9,5,-0.5,round(c(min(pat_2_muts_pooled[, 1:3], na.rm = TRUE), max(pat_2_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=1.7, at=3.3, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 3.8, at=2.7, cex = 1.1, font = 2)

# add NA legend
color.legend(6.5, -0.9, 6.75, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.2, at=5.9, cex = 1.1, font = 2)
legend(x=6.525,y=-0.42,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_2_muts_pooled)
mut_col_labels[1:4] <- c("FGF23\n3\' UTR", 'ATM\nsplice region', "CCND1\n3\' UTR", "ALK\n3\' UTR")
mut_col_labels[5:8] <- c('ATM\nsplice region', 'SLX4\nmissense', 'ATM\nsplice region', 'EP300\nsplice region')


# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_2_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver 1', 'Liver 2', 'Breast 1', 'Breast 2')
axis(2, at = c(0.5, 1.5, 2.5, 3.5, 4.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_2_muts_pooled$AF_liver_1)) - 0.5, 
       y = rep(3.5, sum(is.na(pat_2_muts_pooled$AF_liver_1))), 
       pch = 16)

points(x = which(is.na(pat_2_muts_pooled$AF_liver_2)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_2_muts_pooled$AF_liver_2))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_2_muts_pooled$AF_breast_1)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_2_muts_pooled$AF_breast_1))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_2_muts_pooled$AF_breast_2)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_2_muts_pooled$AF_breast_2))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))



## PATIENT EMA ----
## heart----
# import data
pat_ema_heart <- read.delim('../pat_ema_heart_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_heart <- trusight_from_mutect(pat_ema_heart) #13

## oment 1----
# import data
pat_ema_oment_1 <- read.delim('../pat_ema_oment_1_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_oment_1 <- trusight_from_mutect(pat_ema_oment_1) #19

## oment 2----
# import data
pat_ema_oment_2 <- read.delim('../pat_ema_oment_2_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_oment_2 <- trusight_from_mutect(pat_ema_oment_2) #26

## liver 1----
# import data
pat_ema_liver_1 <- read.delim('../pat_ema_liver_1_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_liver_1 <- trusight_from_mutect(pat_ema_liver_1) #5

## liver 2----
# import data
pat_ema_liver_2 <- read.delim('../pat_ema_liver_2_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_liver_2 <- trusight_from_mutect(pat_ema_liver_2) #12

## left kidney----
# import data
pat_ema_l_kidney <- read.delim('../pat_ema_l_kidney_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_l_kidney <- trusight_from_mutect(pat_ema_l_kidney) #21

## right kidney----
# import data
pat_ema_r_kidney <- read.delim('../pat_ema_r_kidney_bwamem_bqsr_mutect_filtered_liftover_hg19_ann.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_r_kidney <- trusight_from_mutect(pat_ema_r_kidney) #10

## plasma----
# import data
pat_ema_plasma <- read.delim('../pat_ema_plasma_bwamem_bqsr_mutect_ann_hg19.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pat_ema_plasma <- trusight_from_mutect(pat_ema_plasma) #428

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_ema_heart$location, pat_ema_plasma$location)) #12/13
length(intersect(pat_ema_l_kidney$location, pat_ema_plasma$location)) #12/21
length(intersect(pat_ema_liver_1$location, pat_ema_plasma$location)) #3/5
length(intersect(pat_ema_liver_2$location, pat_ema_plasma$location)) #9/12
length(intersect(pat_ema_oment_1$location, pat_ema_plasma$location)) #10/19
length(intersect(pat_ema_oment_2$location, pat_ema_plasma$location)) #18/26
length(intersect(pat_ema_r_kidney$location, pat_ema_plasma$location)) #8/10

#pool mutations from all 7 mets
pat_ema_met_pool <- unique(c(pat_ema_heart$location, pat_ema_oment_1$location, pat_ema_oment_2$location, pat_ema_liver_1$location, pat_ema_liver_2$location, 
                             pat_ema_l_kidney$location, pat_ema_r_kidney$location)) #54 unique mutations

pat_ema_plasma_found <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_met_pool, ] #26 mutations
pat_ema_plasma_found_vars <- (pat_ema_plasma_found$location) # in 17 different genes


# subset to pooled plasma and take a look
pat_ema_heart_pooled <- pat_ema_heart[pat_ema_heart$location %in% pat_ema_plasma_found_vars, ]

pat_ema_l_kidney_pooled <- pat_ema_l_kidney[pat_ema_l_kidney$location %in% pat_ema_plasma_found_vars, ]

pat_ema_r_kidney_pooled <- pat_ema_r_kidney[pat_ema_r_kidney$location %in% pat_ema_plasma_found_vars, ]

pat_ema_liver_1_pooled <- pat_ema_liver_1[pat_ema_liver_1$location %in% pat_ema_plasma_found_vars, ]

pat_ema_liver_2_pooled <- pat_ema_liver_2[pat_ema_liver_2$location %in% pat_ema_plasma_found_vars, ]

pat_ema_oment_1_pooled <- pat_ema_oment_1[pat_ema_oment_1$location %in% pat_ema_plasma_found_vars, ]

pat_ema_oment_2_pooled <- pat_ema_oment_2[pat_ema_oment_2$location %in% pat_ema_plasma_found_vars, ]

pat_ema_plasma_pooled <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_plasma_found_vars, ]


# subset these to tumor AF for new figure
pat_ema_heart_pooled <- pat_ema_heart_pooled[, c('location', 'AF')]
colnames(pat_ema_heart_pooled) <- c('location', 'AF_heart')
pat_ema_heart_pooled$AF_heart[c(4, 7, 8, 9, 12)] <- c(0.279, 0.302, 0.241, 0.256, 0.164)
pat_ema_heart_pooled$AF_heart <- as.numeric(pat_ema_heart_pooled$AF_heart)

pat_ema_l_kidney_pooled <- pat_ema_l_kidney_pooled[, c('location', 'AF')]
colnames(pat_ema_l_kidney_pooled) <- c('location', 'AF_l_kidney')
pat_ema_l_kidney_pooled$AF_l_kidney[c(4, 12)] <- c(0.513, 0.063)
pat_ema_l_kidney_pooled$AF_l_kidney <- as.numeric(pat_ema_l_kidney_pooled$AF_l_kidney)

pat_ema_r_kidney_pooled <- pat_ema_r_kidney_pooled[, c('location', 'AF')]
colnames(pat_ema_r_kidney_pooled) <- c('location', 'AF_r_kidney')
pat_ema_r_kidney_pooled$AF_r_kidney <- as.numeric(pat_ema_r_kidney_pooled$AF_r_kidney)

pat_ema_liver_1_pooled <- pat_ema_liver_1_pooled[, c('location', 'AF')]
colnames(pat_ema_liver_1_pooled) <- c('location', 'AF_liver_1')
pat_ema_liver_1_pooled$AF_liver_1 <- as.numeric(pat_ema_liver_1_pooled$AF_liver_1)

pat_ema_liver_2_pooled <- pat_ema_liver_2_pooled[, c('location', 'AF')]
colnames(pat_ema_liver_2_pooled) <- c('location', 'AF_liver_2')
pat_ema_liver_2_pooled$AF_liver_2[c(6, 9)] <- c(0.332, 0.062)
pat_ema_liver_2_pooled$AF_liver_2 <- as.numeric(pat_ema_liver_2_pooled$AF_liver_2)

pat_ema_oment_1_pooled <- pat_ema_oment_1_pooled[, c('location', 'AF')]
colnames(pat_ema_oment_1_pooled) <- c('location', 'AF_oment_1')
pat_ema_oment_1_pooled$AF_oment_1[c(6, 7, 10)] <- c(0.283, 0.301, 0.088)
pat_ema_oment_1_pooled$AF_oment_1 <- as.numeric(pat_ema_oment_1_pooled$AF_oment_1)

pat_ema_oment_2_pooled <- pat_ema_oment_2_pooled[, c('location', 'AF')]
colnames(pat_ema_oment_2_pooled) <- c('location', 'AF_oment_2')
pat_ema_oment_2_pooled$AF_oment_2[c(3, 4, 6, 8, 10, 11, 16, 17)] <- c(0.189, 0.507, 0.029, 0.181, 0.345, 0.260, 0.147, 0.110)
pat_ema_oment_2_pooled$AF_oment_2 <- as.numeric(pat_ema_oment_2_pooled$AF_oment_2)

pat_ema_plasma_pooled <- pat_ema_plasma_pooled[, c('location', 'AF')]
colnames(pat_ema_plasma_pooled) <- c('location', 'AF_plasma')
pat_ema_plasma_pooled$AF_plasma[c(3, 4, 7, 11, 12, 13, 14, 15, 16, 17, 18, 24, 25)] <- c(0.049, 0.026, 0.424, 0.101, 0.664, 0.252, 0.320, 0.236, 0.258, 0.124, 0.094, 0.197, 0.197)
pat_ema_plasma_pooled$AF_plasma <- as.numeric(pat_ema_plasma_pooled$AF_plasma)


# put them together
pat_ema_muts_pooled <- merge(pat_ema_heart_pooled, pat_ema_l_kidney_pooled, by = 'location', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_r_kidney_pooled, by = 'location', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_liver_1_pooled, by = 'location', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_liver_2_pooled, by = 'location', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_oment_1_pooled, by = 'location', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_oment_2_pooled, by = 'location', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_plasma_pooled, by = 'location', all = TRUE)

pat_ema_muts_pooled <- pat_ema_muts_pooled[order(pat_ema_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_ema_muts_pooled) <- pat_ema_muts_pooled$location
row_muts <- rownames(pat_ema_muts_pooled)
pat_ema_muts_pooled <- pat_ema_muts_pooled[, -1]

pat_ema_muts_pooled <- as.data.frame(pat_ema_muts_pooled)
rownames(pat_ema_muts_pooled) <- row_muts
#plot 

pat_ema_7 <- c(pat_ema_muts_pooled$AF_heart, pat_ema_muts_pooled$AF_l_kidney, pat_ema_muts_pooled$AF_r_kidney, 
               pat_ema_muts_pooled$AF_liver_1, pat_ema_muts_pooled$AF_liver_2, pat_ema_muts_pooled$AF_oment_1, 
               pat_ema_muts_pooled$AF_oment_2)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_ema_muts_pooled)[1] * dim(pat_ema_muts_pooled)[2])
# plasma reds
cell_cols[183:208] <- color.scale(pat_ema_muts_pooled[, 8], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:182] <- color.scale(pat_ema_7, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 26, byrow = FALSE)
pat_ema_pooled_t <- data.frame(t(pat_ema_muts_pooled))
pat_ema_pooled_t <- pat_ema_pooled_t[c(8, 1:7), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(8, 1:7), ]
# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_ema_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_ema_muts_pooled[, 8], na.rm = TRUE),max(pat_ema_muts_pooled[, 8], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(1,-1.7,10,-0.9,round(c(min(pat_ema_muts_pooled[, 8], na.rm = TRUE), max(pat_ema_muts_pooled[, 8], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.2, at=2.0, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_ema_muts_pooled[, 1:7], na.rm = TRUE),max(pat_ema_muts_pooled[, 1:7], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(11,-1.7,20,-0.9,round(c(min(pat_ema_muts_pooled[, 1:7], na.rm = TRUE), max(pat_ema_muts_pooled[, 1:7], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.2, at=12, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 11, cex = 1.1, font = 2)

# add NA legend
color.legend(25, -1.7, 26, -0.9, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.6, at=23.2, cex = 1.1, font = 2)
legend(x=25.2,y=-0.83,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_ema_muts_pooled)
mut_col_labels[1:6] <- c('ROS\nsplice region', "GNA11\n3\' UTR", "CCND1\n3' UTR", "EP300\n3\' UTR", "ALK\n3\' UTR", 'GNAQ\nsplice region')
mut_col_labels[7:12] <- c('ATM\nsplice region', 'BRAF\nsplice region', "MAP2K1\n5\' UTR", 'ATM\nsplice region', "CCND1\n3\' UTR", 'ATM\nsplice region')
mut_col_labels[13:18] <- c('EP300\nsplice region', 'EP300\nsplice region', "BCL2\n3\' UTR", 'NOTCH3\nmissense', "FLT1\n3\' UTR", "CCND1\n3\' UTR")
mut_col_labels[19:24] <- c('BRCA2\nsplice region', 'PIK3CA\nsplice region', 'NOTCH1\nmissense', 'RB1\nsplice region', 'CCND1\nmissense', 'CCND1\nsplice region')
mut_col_labels[25:26] <- c('MSH3\nmissense', 'ALK\nmissense')
# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_ema_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Heart', 'L Kidney', 'R Kidney', 'Liver 1', 'Liver 2', 'Omental 1', 'Omental 2')
axis(2, at = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# heart points
points(x = which(is.na(pat_ema_muts_pooled$AF_heart)) - 0.5, 
       y = rep(6.5, sum(is.na(pat_ema_muts_pooled$AF_heart))), 
       pch = 16)
# l kidney points
points(x = which(is.na(pat_ema_muts_pooled$AF_l_kidney)) - 0.5, 
       y = rep(5.5, sum(is.na(pat_ema_muts_pooled$AF_l_kidney))), 
       pch = 16)
# r_kidney points
points(x = which(is.na(pat_ema_muts_pooled$AF_r_kidney)) - 0.5, 
       y = rep(4.5, sum(is.na(pat_ema_muts_pooled$AF_r_kidney))), 
       pch = 16)
# liver_1 points
points(x = which(is.na(pat_ema_muts_pooled$AF_liver_1)) - 0.5, 
       y = rep(3.5, sum(is.na(pat_ema_muts_pooled$AF_liver_1))), 
       pch = 16)
# liver_2 points
points(x = which(is.na(pat_ema_muts_pooled$AF_liver_2)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_ema_muts_pooled$AF_liver_2))), 
       pch = 16)
# oment_1 points
points(x = which(is.na(pat_ema_muts_pooled$AF_oment_1)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_ema_muts_pooled$AF_oment_1))), 
       pch = 16)
# oment_2 points
points(x = which(is.na(pat_ema_muts_pooled$AF_oment_2)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_ema_muts_pooled$AF_oment_2))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))

## ADD FUNCTION HISTOGRAM

##upset
pat_ema_heart_met_muts <- pat_ema_heart$location
pat_ema_heart_met_muts <- as.data.frame(pat_ema_heart_met_muts)
pat_ema_heart_met_muts$mutation <- rep(1, nrow(pat_ema_heart_met_muts))
colnames(pat_ema_heart_met_muts) <- c('mutation', 'Heart')

pat_ema_l_kidney_met_muts <- pat_ema_l_kidney$location
pat_ema_l_kidney_met_muts <- as.data.frame(pat_ema_l_kidney_met_muts)
pat_ema_l_kidney_met_muts$mutation <- rep(1, nrow(pat_ema_l_kidney_met_muts))
colnames(pat_ema_l_kidney_met_muts) <- c('mutation', 'L_Kidney')

pat_ema_r_kidney_met_muts <- pat_ema_r_kidney$location
pat_ema_r_kidney_met_muts <- as.data.frame(pat_ema_r_kidney_met_muts)
pat_ema_r_kidney_met_muts$mutation <- rep(1, nrow(pat_ema_r_kidney_met_muts))
colnames(pat_ema_r_kidney_met_muts) <- c('mutation', 'R_Kidney')

pat_ema_liver_1_met_muts <- pat_ema_liver_1$location
pat_ema_liver_1_met_muts <- as.data.frame(pat_ema_liver_1_met_muts)
pat_ema_liver_1_met_muts$mutation <- rep(1, nrow(pat_ema_liver_1_met_muts))
colnames(pat_ema_liver_1_met_muts) <- c('mutation', 'Liver_1')

pat_ema_liver_2_met_muts <- pat_ema_liver_2$location
pat_ema_liver_2_met_muts <- as.data.frame(pat_ema_liver_2_met_muts)
pat_ema_liver_2_met_muts$mutation <- rep(1, nrow(pat_ema_liver_2_met_muts))
colnames(pat_ema_liver_2_met_muts) <- c('mutation', 'Liver_2')

pat_ema_oment_1_met_muts <- pat_ema_oment_1$location
pat_ema_oment_1_met_muts <- as.data.frame(pat_ema_oment_1_met_muts)
pat_ema_oment_1_met_muts$mutation <- rep(1, nrow(pat_ema_oment_1_met_muts))
colnames(pat_ema_oment_1_met_muts) <- c('mutation', 'Oment_1')

pat_ema_oment_2_met_muts <- pat_ema_oment_2$location
pat_ema_oment_2_met_muts <- as.data.frame(pat_ema_oment_2_met_muts)
pat_ema_oment_2_met_muts$mutation <- rep(1, nrow(pat_ema_oment_2_met_muts))
colnames(pat_ema_oment_2_met_muts) <- c('mutation', 'Oment_2')

pat_ema_plasma_muts <- pat_ema_plasma$location
pat_ema_plasma_muts <- as.data.frame(pat_ema_plasma_muts)
pat_ema_plasma_muts$mutation <- rep(1, nrow(pat_ema_plasma_muts))
colnames(pat_ema_plasma_muts) <- c('mutation', 'Plasma')


#merge them all together
upset_df <- merge(pat_ema_heart_met_muts, pat_ema_l_kidney_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_r_kidney_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_liver_1_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_liver_2_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_oment_1_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_oment_2_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

rownames(upset_df) <- upset_df$mutation
row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:8])
randomnumber <- round(runif(8, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 6, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(4)
orange_pal <- pal_material('orange', n = 8, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(4)
purple_pal <- pal_material('purple', n= 4, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(2)
bar_colors <- c(orange_pal[1:4], purple_pal[1], 'darkgreen')


upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray50', Heart = 'white', L_Kidney = 'white', 
                                                                 Liver_1 = 'white', Liver_2 = 'white', 
                                                                 Oment_1 = 'white', Oment_2 = 'white', 
                                                                 R_Kidney = 'white')))), 
      intersections = list(list('Heart'), 
                           list('L_Kidney'), 
                           list('Liver_1'), 
                           list('Liver_2'), 
                           list('Oment_1'), 
                           list('Oment_2'), 
                           list('R_Kidney'), 
                           list('Oment_2', 'Plasma'), 
                           list('Liver_1', 'R_Kidney'), 
                           list('Liver_2', 'Oment_1'), 
                           list('Liver_2', 'Oment_2'), 
                           list('L_Kidney', 'Oment_2'), 
                           list('Liver_1', 'Plasma'), 
                           list('R_Kidney', 'Plasma'), 
                           list('Liver_2', 'Plasma'), 
                           list('Heart', 'Plasma'), 
                           list('L_Kidney', 'Plasma'), 
                           list('Heart', 'Oment_2', 'Plasma'), 
                           list('Heart', 'L_Kidney', 'Oment_2'), 
                           list('R_Kidney', 'Oment_1', 'Plasma'), 
                           list('Oment_1', 'Oment_2', 'Plasma'), 
                           list('R_Kidney', 'Oment_1', 'L_Kidney', 'Plasma'), 
                           list('Heart', 'L_Kidney', 'Oment_2', 'Plasma'), 
                           list('R_Kidney', 'Liver_2', 'Heart', 'L_Kidney', 'Plasma'), 
                           list('Liver_2', 'Heart', 'Oment_1', 'Oment_2', 'Plasma'), 
                           list('R_Kidney', 'Liver_2', 'L_Kidney', 'Oment_2', 'Plasma'), 
                           list('Liver_2', 'Heart', 'L_Kidney', 'Oment_2', 'Plasma'), 
                           list('Liver_2', 'Oment_1', 'L_Kidney', 'Oment_2', 'Plasma'), 
                           list('Heart', 'Oment_1', 'L_Kidney', 'Oment_2', 'Plasma'), 
                           list('Liver_1', 'Heart', 'Oment_1', 'L_Kidney', 'Oment_2', 'Plasma'),
                           list('R_Kidney', 'Liver_2', 'Heart', 'Oment_1', 'L_Kidney', 'Oment_2', 'Plasma'), 
                           list('Liver_1', 'R_Kidney', 'Liver_2', 'Heart', 'Oment_1', 'L_Kidney', 'Oment_2', 'Plasma')),
      nsets = 8, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', sets.bar.color = c('gray60', rep('firebrick3', 2), 'chocolate3', rep('lawngreen', 2), 'chocolate3', 'red'), 
      matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, main.bar.color = c(rep('black', 7), 'darkblue', rep('dodgerblue2', 9), 'darkred', rep('red', 3), rep('darkgreen', 2), rep('orange', 6), 'green', 'purple', 'black'), mainbar.y.label = 'Number of Mutations\nin Common', 
      set_size.show = TRUE, set_size.numbers_size = 5, text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 1.3))




pat_9_oment_met_muts <- pat_9_oment$location
pat_9_oment_met_muts <- as.data.frame(pat_9_oment_met_muts)
pat_9_oment_met_muts$mutation <- rep(1, nrow(pat_9_oment_met_muts))
colnames(pat_9_oment_met_muts) <- c('mutation', 'Omental_Met')

pat_9_ovary_met_muts <- pat_9_ovary$location
pat_9_ovary_met_muts <- as.data.frame(pat_9_ovary_met_muts)
pat_9_ovary_met_muts$mutation <- rep(1, nrow(pat_9_ovary_met_muts))
colnames(pat_9_ovary_met_muts) <- c('mutation', 'Ovary_Met')

pat_9_lymph_met_muts <- pat_9_ln$location
pat_9_lymph_met_muts <- as.data.frame(pat_9_lymph_met_muts)
pat_9_lymph_met_muts$mutation <- rep(1, nrow(pat_9_lymph_met_muts))
colnames(pat_9_lymph_met_muts) <- c('mutation', 'Lymph_Met')

pat_9_plasma_muts <- pat_9_plasma$location
pat_9_plasma_muts <- as.data.frame(pat_9_plasma_muts)
pat_9_plasma_muts$mutation <- rep(1, nrow(pat_9_plasma_muts))
colnames(pat_9_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_9_oment_met_muts, pat_9_ovary_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_lymph_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
upset_df$mutation <- as.character(upset_df$mutation)

pat_9_pl_ov <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov <- pat_9_pl_ov$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ovary_met_af[2] <- 0.246
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_plasma_af[2] <- 0.281
pat_9_pl_ov_af <- c(pat_9_ovary_met_af, pat_9_plasma_af)
pat_9_pl_ov_sd <- sd(as.numeric(pat_9_pl_ov_af))
pat_9_pl_ov_af <- median(as.numeric(pat_9_pl_ov_af))

pat_9_pl_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_pl_ly <- pat_9_pl_ly$mutation
pat_9_lymph_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ly, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ly, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_pl_ly_af <- c(pat_9_lymph_met_af, pat_9_plasma_af)
pat_9_pl_ly_sd <- sd(as.numeric(pat_9_pl_ly_af))
pat_9_pl_ly_af <- median(as.numeric(pat_9_pl_ly_af))


pat_9_pl_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_pl_om <- pat_9_pl_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_oment_met_af[2] <- 0.046
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_plasma_af[2] <- 0.057
pat_9_pl_om_af <- c(pat_9_oment_met_af, pat_9_plasma_af)
pat_9_pl_om_sd <- sd(as.numeric(pat_9_pl_om_af))
pat_9_pl_om_af <- median(as.numeric(pat_9_pl_om_af))

pat_9_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_ov_om <- pat_9_ov_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ov_om_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af)
pat_9_ov_om_sd <- sd(as.numeric(pat_9_ov_om_af))
pat_9_ov_om_af <- median(as.numeric(pat_9_ov_om_af))
#FIX THE REST OF THESE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_9_pl_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov_om <- pat_9_pl_ov_om$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ovary_met_af[9] <- 0.330
pat_9_ovary_met_af[12] <- 0.143
pat_9_ovary_met_af[13] <- 0.465
pat_9_ovary_met_af[16] <- 0.209
pat_9_ovary_met_af[17] <- 0.165
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_plasma_af[c(1, 9, 11, 12, 13, 16, 17)] <- c(0.20, 0.333, 0.230, 0.332, 0.260, 0.209, 0.177)


pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_oment_met_af[c(9, 13, 16)] <- c(0.343, 0.203, 0.295)
pat_9_pl_ov_om_af <- c(pat_9_ovary_met_af, pat_9_plasma_af, pat_9_oment_met_af)
pat_9_pl_ov_om_sd <- sd(as.numeric(pat_9_pl_ov_om_af))
pat_9_pl_ov_om_af <- median(as.numeric(pat_9_pl_ov_om_af))

pat_9_all <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
pat_9_all <- pat_9_all$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_all, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_all, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_all, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_lymph_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_all, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$AF
pat_9_all_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af, pat_9_plasma_af, pat_9_lymph_met_af)
pat_9_all_sd <- sd(as.numeric(pat_9_all_af))
pat_9_all_af <- median(as.numeric(pat_9_all_af))


pat_9_afs <- c(pat_9_pl_ly_af, pat_9_pl_ov_af, pat_9_pl_om_af, pat_9_ov_om_af, pat_9_pl_ov_om_af, pat_9_all_af)
pat_9_sds <- c(pat_9_pl_ly_sd, pat_9_pl_ov_sd, pat_9_pl_om_sd, pat_9_ov_om_sd, pat_9_pl_ov_om_sd, pat_9_all_sd)
pat_9_labels <- c('Lymph Met\n+ Plasma', 'Ovary Met\n+ Plasma', 'Omental Met\n+ Plasma', 'Omental Met\n+ Ovary Met', 
                  'Omental Met\n+ Ovary Met\n+ Plasma', 'All Samples')
pat_9_afs <- data.frame(pat_9_labels, pat_9_afs, pat_9_sds, bar_colors)
p<-ggplot(data=pat_9_afs, aes(x=pat_9_labels, y=pat_9_afs, fill = bar_colors)) +
  geom_bar(stat="identity", width = 0.5, color = 'black') + scale_x_discrete(limits=pat_9_afs$pat_9_labels) + ylab('Median Mutant Allele Frequency') +
  xlab('Intersections') + theme_bw() + 
  theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) + 
  scale_fill_manual("legend", values = bar_colors[c(5,1:3,4,6)]) +
  geom_errorbar(aes(ymin=pat_9_afs, ymax=pat_9_afs+pat_9_sds), width=.2,
                position=position_dodge(.9)) + theme(axis.text=element_text(size=14, face = 'bold'),
                                                     axis.title=element_text(size=16,face="bold"))

p



# patient 9
pat_ema_heart_met_stats <- pat_ema_heart[, c('location', 'AF')]
pat_ema_heart_met_stats$AF[c(4, 7, 8, 9, 12)] <- c(0.279, 0.302, 0.241, 0.256, 0.164)
pat_ema_l_kidney_met_stats <- pat_ema_l_kidney[, c('location', 'AF')]
pat_ema_l_kidney_met_stats$AF[c(7, 21)] <- c(0.513, 0.063)
pat_ema_r_kidney_met_stats <- pat_ema_r_kidney[, c('location', 'AF')]
pat_ema_liver_1_met_stats <- pat_ema_liver_1[, c('location', 'AF')]
pat_ema_liver_2_met_stats <- pat_ema_liver_2[, c('location', 'AF')]
pat_ema_liver_2_met_stats$AF[c(8, 12)] <- c(0.332, 0.062)
pat_ema_oment_1_met_stats <- pat_ema_oment_1[, c('location', 'AF')]
pat_ema_oment_1_met_stats$AF[c(10, 11, 19)] <- c(0.283, 0.301, 0.088)
pat_ema_oment_2_met_stats <- pat_ema_oment_2[, c('location', 'AF')]
pat_ema_oment_2_met_stats$AF[c(5, 7, 9, 11, 13, 14, 21, 24)] <- c(0.189, 0.507, 0.029, 0.181, 0.345, 0.260, 0.147, 0.110)

pat_ema_all <- rbind(pat_ema_heart_met_stats, pat_ema_l_kidney_met_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_r_kidney_met_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_liver_1_met_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_liver_2_met_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_oment_1_met_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_oment_2_met_stats)

pat_ema_all$AF <- as.numeric(pat_ema_all$AF)
pat_ema_all <- pat_ema_all[order(pat_ema_all$AF, decreasing = TRUE), ]
pat_ema_all$color <- ifelse(pat_ema_all$location %in% pat_ema_plasma$location, 'red', 'black')
barplot(pat_ema_all$AF, col = pat_ema_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient EMA\n(7 tumors)', ylim = c(0,1.0))






## CHECKING HOW TO DRAW EVOL CLADE IN EMA ----

#how many in all? this is the length of the main branch moving from normal
pat_ema_all_common <- Reduce(intersect, list(pat_ema_heart$location, pat_ema_l_kidney$location, pat_ema_r_kidney$location, pat_ema_liver_1$location, 
                                             pat_ema_liver_2$location, pat_ema_oment_1$location, pat_ema_oment_2$location)) #1
#which ones?
pat_ema_all <- pat_ema_heart[pat_ema_heart$location %in% pat_ema_all_common, ]
#which in plasma?
length(intersect(pat_ema_all_common, pat_ema_plasma$location)) #all 1

#exclude those, these can be the lengths of the next branches if split solo
pat_ema_heart_not_all <- pat_ema_heart[pat_ema_heart$location %!in% pat_ema_all_common, ] #12
pat_ema_l_kidney_not_all <- pat_ema_l_kidney[pat_ema_l_kidney$location %!in% pat_ema_all_common, ] #20
pat_ema_r_kidney_not_all <- pat_ema_r_kidney[pat_ema_r_kidney$location %!in% pat_ema_all_common, ] #9
pat_ema_liver_1_not_all <- pat_ema_liver_1[pat_ema_liver_1$location %!in% pat_ema_all_common, ] #4
pat_ema_liver_2_not_all <- pat_ema_liver_2[pat_ema_liver_2$location %!in% pat_ema_all_common, ] #11
pat_ema_oment_1_not_all <- pat_ema_oment_1[pat_ema_oment_1$location %!in% pat_ema_all_common, ] #18
pat_ema_oment_2_not_all <- pat_ema_oment_2[pat_ema_oment_2$location %!in% pat_ema_all_common, ] #25

#gather all locations left after all common
pat_ema_2nd_locations <- c(pat_ema_heart_not_all$location, pat_ema_l_kidney_not_all$location, pat_ema_r_kidney_not_all$location, 
                           pat_ema_liver_1_not_all$location, pat_ema_liver_2_not_all$location, pat_ema_oment_1_not_all$location, pat_ema_oment_2_not_all$location)
pat_ema_2nd_locations <- table(pat_ema_2nd_locations)
pat_ema_2nd_locations <- pat_ema_2nd_locations[order(pat_ema_2nd_locations, decreasing = TRUE)]

#only keep those that hit 6 out of 7
pat_ema_2nd_locations_6 <- pat_ema_2nd_locations[pat_ema_2nd_locations == 6]
pat_ema_2nd_locations_6 <- names(pat_ema_2nd_locations_6) #2

#subset to those to see how many are in each
pat_ema_heart_6_of_7 <- pat_ema_heart_not_all[pat_ema_heart_not_all$location %in% pat_ema_2nd_locations_6, ] #2/2
pat_ema_l_kidney_6_of_7 <- pat_ema_l_kidney_not_all[pat_ema_l_kidney_not_all$location %in% pat_ema_2nd_locations_6, ] #2/2
pat_ema_r_kidney_6_of_7 <- pat_ema_r_kidney_not_all[pat_ema_r_kidney_not_all$location %in% pat_ema_2nd_locations_6, ] #2/2
pat_ema_liver_1_6_of_7 <- pat_ema_liver_1_not_all[pat_ema_liver_1_not_all$location %in% pat_ema_2nd_locations_6, ] #0/2
pat_ema_liver_2_6_of_7 <- pat_ema_liver_2_not_all[pat_ema_liver_2_not_all$location %in% pat_ema_2nd_locations_6, ] #2/2
pat_ema_oment_1_6_of_7 <- pat_ema_oment_1_not_all[pat_ema_oment_1_not_all$location %in% pat_ema_2nd_locations_6, ] #2/2
pat_ema_oment_2_6_of_7 <- pat_ema_oment_2_not_all[pat_ema_oment_2_not_all$location %in% pat_ema_2nd_locations_6, ] #2/2


#liver 1 is clearly the least like the others, branch it off here and continue with others

#ln excluding commons, length of ln solo branch
pat_ema_liver_1_only <- pat_ema_liver_1_not_all[pat_ema_liver_1_not_all$location %!in% pat_ema_2nd_locations_6, ] #4
# how many of these in plasma?
length(intersect(pat_ema_liver_1_only$location, pat_ema_plasma$location)) #2


#other 6 but not liv 1
pat_ema_heart_not_liv_1 <- pat_ema_heart_6_of_7[pat_ema_heart_6_of_7$location %!in% pat_ema_liver_1_6_of_7$location, ] #2
pat_ema_l_kidney_not_liv_1 <- pat_ema_l_kidney_6_of_7[pat_ema_l_kidney_6_of_7$location %!in% pat_ema_liver_1_6_of_7$location, ] #2
pat_ema_r_kidney_not_liv_1 <- pat_ema_r_kidney_6_of_7[pat_ema_r_kidney_6_of_7$location %!in% pat_ema_liver_1_6_of_7$location, ] #2
pat_ema_liver_2_not_liv_1 <- pat_ema_liver_2_6_of_7[pat_ema_liver_2_6_of_7$location %!in% pat_ema_liver_1_6_of_7$location, ] #2
pat_ema_oment_1_not_liv_1 <- pat_ema_oment_1_6_of_7[pat_ema_oment_1_6_of_7$location %!in% pat_ema_liver_1_6_of_7$location, ] #2
pat_ema_oment_2_not_liv_1 <- pat_ema_oment_2_6_of_7[pat_ema_oment_2_6_of_7$location %!in% pat_ema_liver_1_6_of_7$location, ] #2

#in common? (should be length of the 2 just made)
length(Reduce(intersect, list(pat_ema_heart_not_liv_1$location, pat_ema_l_kidney_not_liv_1$location, pat_ema_r_kidney_not_liv_1$location, 
                              pat_ema_liver_2_not_liv_1$location, pat_ema_oment_1_not_liv_1$location, pat_ema_oment_2_not_liv_1$location))) #all 2
# how many of these in plasma?
length(Reduce(intersect, list(pat_ema_heart_not_liv_1$location, pat_ema_l_kidney_not_liv_1$location, pat_ema_r_kidney_not_liv_1$location, 
                              pat_ema_liver_2_not_liv_1$location, pat_ema_oment_1_not_liv_1$location, pat_ema_oment_2_not_liv_1$location, 
                              pat_ema_plasma$location))) #all 2

#subset out liver 1 singles
pat_ema_heart_not_liv_1_long <- pat_ema_heart_not_all[pat_ema_heart_not_all$location %!in% pat_ema_liver_1_not_all$location, ] #11
pat_ema_l_kidney_not_liv_1_long <- pat_ema_l_kidney_not_all[pat_ema_l_kidney_not_all$location %!in% pat_ema_liver_1_not_all$location, ] #19
pat_ema_r_kidney_not_liv_1_long <- pat_ema_r_kidney_not_all[pat_ema_r_kidney_not_all$location %!in% pat_ema_liver_1_not_all$location, ] #8
pat_ema_liver_2_not_liv_1_long <- pat_ema_liver_2_not_all[pat_ema_liver_2_not_all$location %!in% pat_ema_liver_1_not_all$location, ] #11
pat_ema_oment_1_not_liv_1_long <- pat_ema_oment_1_not_all[pat_ema_oment_1_not_all$location %!in% pat_ema_liver_1_not_all$location, ] #17
pat_ema_oment_2_not_liv_1_long <- pat_ema_oment_2_not_all[pat_ema_oment_2_not_all$location %!in% pat_ema_liver_1_not_all$location, ] #24

#length of common branch
pat_ema_not_liv_1 <- Reduce(intersect, list(pat_ema_heart_not_liv_1_long$location, pat_ema_l_kidney_not_liv_1_long$location, 
                                            pat_ema_r_kidney_not_liv_1_long$location, pat_ema_liver_2_not_liv_1_long$location, 
                                            pat_ema_oment_1_not_liv_1_long$location, pat_ema_oment_2_not_liv_1_long$location)) #2

#subset out those for next branch
pat_ema_heart_next_length <- pat_ema_heart_not_liv_1_long[pat_ema_heart_not_liv_1_long$location %!in% pat_ema_not_liv_1, ] #9
pat_ema_l_kidney_next_length <- pat_ema_l_kidney_not_liv_1_long[pat_ema_l_kidney_not_liv_1_long$location %!in% pat_ema_not_liv_1, ] #17
pat_ema_r_kidney_next_length <- pat_ema_r_kidney_not_liv_1_long[pat_ema_r_kidney_not_liv_1_long$location %!in% pat_ema_not_liv_1, ] #6
pat_ema_liver_2_next_length <- pat_ema_liver_2_not_liv_1_long[pat_ema_liver_2_not_liv_1_long$location %!in% pat_ema_not_liv_1, ] #9
pat_ema_oment_1_next_length <- pat_ema_oment_1_not_liv_1_long[pat_ema_oment_1_not_liv_1_long$location %!in% pat_ema_not_liv_1, ] #15
pat_ema_oment_2_next_length <- pat_ema_oment_2_not_liv_1_long[pat_ema_oment_2_not_liv_1_long$location %!in% pat_ema_not_liv_1, ] #22

#gather all locations left after all common
pat_ema_3rd_locations <- c(pat_ema_heart_next_length$location, pat_ema_l_kidney_next_length$location, pat_ema_r_kidney_next_length$location, 
                           pat_ema_liver_2_next_length$location, pat_ema_oment_1_next_length$location, pat_ema_oment_2_next_length$location)
pat_ema_3rd_locations <- table(pat_ema_3rd_locations)
pat_ema_3rd_locations <- pat_ema_3rd_locations[order(pat_ema_3rd_locations, decreasing = TRUE)]

#only keep those that hit 5 out of 6
pat_ema_3rd_locations_5 <- pat_ema_3rd_locations[pat_ema_3rd_locations == 4] ##SPLITTING OFF 2!!!!
pat_ema_3rd_locations_5 <- names(pat_ema_3rd_locations_5) #6

#subset to those to see how many are in each #FIX THIS NOTATION LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_ema_heart_5_of_6 <- pat_ema_heart_next_length[pat_ema_heart_next_length$location %in% pat_ema_3rd_locations_5, ] #4/6
pat_ema_l_kidney_5_of_6 <- pat_ema_l_kidney_next_length[pat_ema_l_kidney_next_length$location %in% pat_ema_3rd_locations_5, ] #5/6
pat_ema_r_kidney_5_of_6 <- pat_ema_r_kidney_next_length[pat_ema_r_kidney_next_length$location %in% pat_ema_3rd_locations_5, ] #2/6
pat_ema_liver_2_5_of_6 <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %in% pat_ema_3rd_locations_5, ] #5/6
pat_ema_oment_1_5_of_6 <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %in% pat_ema_3rd_locations_5, ] #3/6
pat_ema_oment_2_5_of_6 <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %in% pat_ema_3rd_locations_5, ] #5/6


#r kidney is clearlythe least like the others, branch it off here and continue with others

#r kidney excluding commons, length of r kidney solo branch
pat_ema_r_kidney_only <- pat_ema_r_kidney_next_length[pat_ema_r_kidney_next_length$location %!in% pat_ema_3rd_locations_5, ] #4
# how many of these in plasma?
length(intersect(pat_ema_r_kidney_only$location, pat_ema_plasma$location)) #14


#other 5 but not r kid FIX THIS NOTATION LATER!!!!!!!!!!!!!!!!!!!!!!!!!
pat_ema_heart_not_r_kid <- pat_ema_heart_5_of_6[pat_ema_heart_5_of_6$location %!in% pat_ema_r_kidney_5_of_6$location, ] #3
pat_ema_l_kidney_not_r_kid <- pat_ema_l_kidney_5_of_6[pat_ema_l_kidney_5_of_6$location %!in% pat_ema_r_kidney_5_of_6$location, ] #3
pat_ema_liver_2_not_r_kid <- pat_ema_liver_2_5_of_6[pat_ema_liver_2_5_of_6$location %!in% pat_ema_r_kidney_5_of_6$location, ] #3
pat_ema_oment_1_not_r_kid <- pat_ema_oment_1_5_of_6[pat_ema_oment_1_5_of_6$location %!in% pat_ema_r_kidney_5_of_6$location, ] #3
pat_ema_oment_2_not_r_kid <- pat_ema_oment_2_5_of_6[pat_ema_oment_2_5_of_6$location %!in% pat_ema_r_kidney_5_of_6$location, ] #3

#in common? (should be length of the 2 just made)
length(Reduce(intersect, list(pat_ema_heart_not_r_kid$location, pat_ema_l_kidney_not_r_kid$location, pat_ema_liver_2_not_r_kid$location, 
                              pat_ema_oment_1_not_r_kid$location, pat_ema_oment_2_not_r_kid$location))) #NONE -- DO PAIRWISE

length(intersect(pat_ema_heart_not_r_kid$location, pat_ema_l_kidney_not_r_kid$location)) #2
length(intersect(pat_ema_heart_not_r_kid$location, pat_ema_liver_2_not_r_kid$location)) #2
length(intersect(pat_ema_heart_not_r_kid$location, pat_ema_oment_1_not_r_kid$location)) #2
length(intersect(pat_ema_heart_not_r_kid$location, pat_ema_oment_2_not_r_kid$location)) #3
length(intersect(pat_ema_l_kidney_not_r_kid$location, pat_ema_liver_2_not_r_kid$location)) #2
length(intersect(pat_ema_l_kidney_not_r_kid$location, pat_ema_oment_1_not_r_kid$location)) #2
length(intersect(pat_ema_l_kidney_not_r_kid$location, pat_ema_oment_2_not_r_kid$location)) #3
length(intersect(pat_ema_liver_2_not_r_kid$location, pat_ema_oment_1_not_r_kid$location)) #2
length(intersect(pat_ema_liver_2_not_r_kid$location, pat_ema_oment_2_not_r_kid$location)) #3
length(intersect(pat_ema_oment_1_not_r_kid$location, pat_ema_oment_2_not_r_kid$location)) #3

# how many of these in plasma?
length(Reduce(intersect, list(pat_ema_heart_not_r_kid$location, pat_ema_l_kidney_not_r_kid$location, pat_ema_liver_2_not_r_kid$location, 
                              pat_ema_oment_1_not_r_kid$location, pat_ema_oment_2_not_r_kid$location, pat_ema_plasma$location))) #all 4

#subset out r kidney singles
pat_ema_heart_not_r_kid_long <- pat_ema_heart_next_length[pat_ema_heart_next_length$location %!in% pat_ema_r_kidney_next_length$location, ] #183
pat_ema_l_kidney_not_r_kid_long <- pat_ema_l_kidney_next_length[pat_ema_l_kidney_next_length$location %!in% pat_ema_r_kidney_next_length$location, ] #193
pat_ema_liver_2_not_r_kid_long <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %!in% pat_ema_r_kidney_next_length$location, ] #200
pat_ema_oment_1_not_r_kid_long <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %!in% pat_ema_r_kidney_next_length$location, ] #283
pat_ema_oment_2_not_r_kid_long <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %!in% pat_ema_r_kidney_next_length$location, ] #324

#length of common branch
pat_ema_not_r_kid <- Reduce(intersect, list(pat_ema_heart_not_r_kid_long$location, pat_ema_l_kidney_not_r_kid_long$location, 
                                            pat_ema_liver_2_not_r_kid_long$location, 
                                            pat_ema_oment_1_not_r_kid_long$location, pat_ema_oment_2_not_r_kid_long$location)) #4

#subset out those for next branch
pat_ema_heart_next_length <- pat_ema_heart_not_r_kid_long[pat_ema_heart_not_r_kid_long$location %!in% pat_ema_not_r_kid, ] #179
pat_ema_l_kidney_next_length <- pat_ema_l_kidney_not_r_kid_long[pat_ema_l_kidney_not_r_kid_long$location %!in% pat_ema_not_r_kid, ] #189
pat_ema_liver_2_next_length <- pat_ema_liver_2_not_r_kid_long[pat_ema_liver_2_not_r_kid_long$location %!in% pat_ema_not_r_kid, ] #196
pat_ema_oment_1_next_length <- pat_ema_oment_1_not_r_kid_long[pat_ema_oment_1_not_r_kid_long$location %!in% pat_ema_not_r_kid, ] #279
pat_ema_oment_2_next_length <- pat_ema_oment_2_not_r_kid_long[pat_ema_oment_2_not_r_kid_long$location %!in% pat_ema_not_r_kid, ] #320

#gather all locations left after all common
pat_ema_4th_locations <- c(pat_ema_heart_next_length$location, pat_ema_l_kidney_next_length$location, 
                           pat_ema_liver_2_next_length$location, pat_ema_oment_1_next_length$location, pat_ema_oment_2_next_length$location)
pat_ema_4th_locations <- table(pat_ema_4th_locations)
pat_ema_4th_locations <- pat_ema_4th_locations[order(pat_ema_4th_locations, decreasing = TRUE)]

#only keep those that hit 4 out of 5
pat_ema_4th_locations_4 <- pat_ema_4th_locations[pat_ema_4th_locations == 4]
pat_ema_4th_locations_4 <- names(pat_ema_4th_locations_4) #21

#subset to those to see how many are in each
pat_ema_heart_4_of_5 <- pat_ema_heart_next_length[pat_ema_heart_next_length$location %in% pat_ema_4th_locations_4, ] #14/21
pat_ema_l_kidney_4_of_5 <- pat_ema_l_kidney_next_length[pat_ema_l_kidney_next_length$location %in% pat_ema_4th_locations_4, ] #16/21
pat_ema_liver_2_4_of_5 <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %in% pat_ema_4th_locations_4, ] #15/21
pat_ema_oment_1_4_of_5 <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %in% pat_ema_4th_locations_4, ] #19/21
pat_ema_oment_2_4_of_5 <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %in% pat_ema_4th_locations_4, ] #20/21

#heart is clearlythe least like the others, branch it off here and continue with others

#heart excluding commons, length of heart solo branch
pat_ema_heart_only <- pat_ema_heart_next_length[pat_ema_heart_next_length$location %!in% pat_ema_4th_locations_4, ] #165
# how many of these in plasma?
length(intersect(pat_ema_heart_only$location, pat_ema_plasma$location)) #52

#other 4 but not heart
pat_ema_l_kidney_not_heart <- pat_ema_l_kidney_5_of_6[pat_ema_l_kidney_5_of_6$location %!in% pat_ema_heart_5_of_6$location, ] #2
pat_ema_liver_2_not_heart <- pat_ema_liver_2_5_of_6[pat_ema_liver_2_5_of_6$location %!in% pat_ema_heart_5_of_6$location, ] #2
pat_ema_oment_1_not_heart <- pat_ema_oment_1_5_of_6[pat_ema_oment_1_5_of_6$location %!in% pat_ema_heart_5_of_6$location, ] #2
pat_ema_oment_2_not_heart <- pat_ema_oment_2_5_of_6[pat_ema_oment_2_5_of_6$location %!in% pat_ema_heart_5_of_6$location, ] #2

#in common? (should be length of the 2 just made)
length(Reduce(intersect, list(pat_ema_l_kidney_not_heart$location, pat_ema_liver_2_not_heart$location, 
                              pat_ema_oment_1_not_heart$location, pat_ema_oment_2_not_heart$location))) #all 2
# how many of these in plasma?
length(Reduce(intersect, list(pat_ema_l_kidney_not_heart$location, pat_ema_liver_2_not_heart$location, 
                              pat_ema_oment_1_not_heart$location, pat_ema_oment_2_not_heart$location, pat_ema_plasma$location))) #all 2

#subset out heart singles
pat_ema_l_kidney_not_heart_long <- pat_ema_l_kidney_next_length[pat_ema_l_kidney_next_length$location %!in% pat_ema_heart_next_length$location, ] #189
pat_ema_liver_2_not_heart_long <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %!in% pat_ema_heart_next_length$location, ] #196
pat_ema_oment_1_not_heart_long <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %!in% pat_ema_heart_next_length$location, ] #279
pat_ema_oment_2_not_heart_long <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %!in% pat_ema_heart_next_length$location, ] #320

#length of common branch
pat_ema_not_heart <- Reduce(intersect, list(pat_ema_l_kidney_not_heart_long$location, 
                                            pat_ema_liver_2_not_heart_long$location, 
                                            pat_ema_oment_1_not_heart_long$location, pat_ema_oment_2_not_heart_long$location)) #7

#subset out those for next branch
pat_ema_l_kidney_next_length <- pat_ema_l_kidney_not_heart_long[pat_ema_l_kidney_not_heart_long$location %!in% pat_ema_not_heart, ] #182
pat_ema_liver_2_next_length <- pat_ema_liver_2_not_heart_long[pat_ema_liver_2_not_heart_long$location %!in% pat_ema_not_heart, ] #189
pat_ema_oment_1_next_length <- pat_ema_oment_1_not_heart_long[pat_ema_oment_1_not_heart_long$location %!in% pat_ema_not_heart, ] #272
pat_ema_oment_2_next_length <- pat_ema_oment_2_not_heart_long[pat_ema_oment_2_not_heart_long$location %!in% pat_ema_not_heart, ] #313

#gather all locations left after all common
pat_ema_5th_locations <- c(pat_ema_l_kidney_next_length$location, 
                           pat_ema_liver_2_next_length$location, pat_ema_oment_1_next_length$location, pat_ema_oment_2_next_length$location)
pat_ema_5th_locations <- table(pat_ema_5th_locations)
pat_ema_5th_locations <- pat_ema_5th_locations[order(pat_ema_5th_locations, decreasing = TRUE)]

#only keep those that hit 3 out of 4
pat_ema_5th_locations_3 <- pat_ema_5th_locations[pat_ema_5th_locations == 3]
pat_ema_5th_locations_3 <- names(pat_ema_5th_locations_3) #17

#subset to those to see how many are in each
pat_ema_l_kidney_3_of_4 <- pat_ema_l_kidney_next_length[pat_ema_l_kidney_next_length$location %in% pat_ema_5th_locations_3, ] #10/17
pat_ema_liver_2_3_of_4 <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %in% pat_ema_5th_locations_3, ] #10/17
pat_ema_oment_1_3_of_4 <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %in% pat_ema_5th_locations_3, ] #15/17
pat_ema_oment_2_3_of_4 <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %in% pat_ema_5th_locations_3, ] #16/17

#l kidney is clearlythe least like the others, branch it off here and continue with others

#l kidney excluding commons, length of l kidney solo branch
pat_ema_l_kidney_only <- pat_ema_l_kidney_next_length[pat_ema_l_kidney_next_length$location %!in% pat_ema_5th_locations_3, ] #172
# how many of these in plasma?
length(intersect(pat_ema_l_kidney_only$location, pat_ema_plasma$location)) #20

#other 3 but not l kid
pat_ema_liver_2_not_l_kid <- pat_ema_liver_2_3_of_4[pat_ema_liver_2_3_of_4$location %!in% pat_ema_l_kidney_3_of_4$location, ] #7
pat_ema_oment_1_not_l_kid <- pat_ema_oment_1_3_of_4[pat_ema_oment_1_3_of_4$location %!in% pat_ema_l_kidney_3_of_4$location, ] #7
pat_ema_oment_2_not_l_kid <- pat_ema_oment_2_3_of_4[pat_ema_oment_2_3_of_4$location %!in% pat_ema_l_kidney_3_of_4$location, ] #7

#in common? (should be length of the 2 just made)
length(Reduce(intersect, list(pat_ema_liver_2_not_l_kid$location, 
                              pat_ema_oment_1_not_l_kid$location, pat_ema_oment_2_not_l_kid$location))) #all 7
# how many of these in plasma?
length(Reduce(intersect, list(pat_ema_liver_2_not_l_kid$location, 
                              pat_ema_oment_1_not_l_kid$location, pat_ema_oment_2_not_l_kid$location, pat_ema_plasma$location))) #6

#subset out l kidney singles
pat_ema_liver_2_not_l_kid_long <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %!in% pat_ema_l_kidney_next_length$location, ] #182
pat_ema_oment_1_not_l_kid_long <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %!in% pat_ema_l_kidney_next_length$location, ] #259
pat_ema_oment_2_not_l_kid_long <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %!in% pat_ema_l_kidney_next_length$location, ] #296

#length of common branch
pat_ema_not_l_kid <- Reduce(intersect, list(pat_ema_liver_2_not_l_kid_long$location, 
                                            pat_ema_oment_1_not_l_kid_long$location, pat_ema_oment_2_not_l_kid_long$location)) #7

#subset out those for next branch
pat_ema_liver_2_next_length <- pat_ema_liver_2_not_l_kid_long[pat_ema_liver_2_not_l_kid_long$location %!in% pat_ema_not_l_kid, ] #175
pat_ema_oment_1_next_length <- pat_ema_oment_1_not_l_kid_long[pat_ema_oment_1_not_l_kid_long$location %!in% pat_ema_not_l_kid, ] #252
pat_ema_oment_2_next_length <- pat_ema_oment_2_not_l_kid_long[pat_ema_oment_2_not_l_kid_long$location %!in% pat_ema_not_l_kid, ] #289

#gather all locations left after all common
pat_ema_6th_locations <- c(pat_ema_liver_2_next_length$location, pat_ema_oment_1_next_length$location, pat_ema_oment_2_next_length$location)
pat_ema_6th_locations <- table(pat_ema_6th_locations)
pat_ema_6th_locations <- pat_ema_6th_locations[order(pat_ema_6th_locations, decreasing = TRUE)]

#only keep those that hit 2 out of 3
pat_ema_6th_locations_2 <- pat_ema_6th_locations[pat_ema_6th_locations == 2]
pat_ema_6th_locations_2 <- names(pat_ema_6th_locations_2) #34

#subset to those to see how many are in each
pat_ema_liver_2_2_of_3 <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %in% pat_ema_6th_locations_2, ] #16/34
pat_ema_oment_1_2_of_3 <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %in% pat_ema_6th_locations_2, ] #28/34
pat_ema_oment_2_2_of_3 <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %in% pat_ema_6th_locations_2, ] #24/34

#l kidney is clearlythe least like the others, branch it off here and continue with others

#liver_2 excluding commons, length of liver 2 solo branch
pat_ema_liver_2_only <- pat_ema_liver_2_next_length[pat_ema_liver_2_next_length$location %!in% pat_ema_6th_locations_2, ] #159
# how many of these in plasma?
length(intersect(pat_ema_liver_2_only$location, pat_ema_plasma$location)) #14

#other 2 but not liver 2
pat_ema_oment_1_not_liv_2 <- pat_ema_oment_1_2_of_3[pat_ema_oment_1_2_of_3$location %!in% pat_ema_liver_2_2_of_3$location, ] #18
pat_ema_oment_2_not_liv_2 <- pat_ema_oment_2_2_of_3[pat_ema_oment_2_2_of_3$location %!in% pat_ema_liver_2_2_of_3$location, ] #18

#in common? (should be length of the 2 just made)
length(Reduce(intersect, list(pat_ema_oment_1_not_liv_2$location, pat_ema_oment_2_not_liv_2$location))) #all 18
# how many of these in plasma?
length(Reduce(intersect, list(pat_ema_oment_1_not_liv_2$location, pat_ema_oment_2_not_liv_2$location, pat_ema_plasma$location))) #11

#subset out liver 2 singles
pat_ema_oment_1_not_liv_2_long <- pat_ema_oment_1_next_length[pat_ema_oment_1_next_length$location %!in% pat_ema_liver_2_next_length$location, ] #242
pat_ema_oment_2_not_liv_2_long <- pat_ema_oment_2_next_length[pat_ema_oment_2_next_length$location %!in% pat_ema_liver_2_next_length$location, ] #283

#length of common branch
pat_ema_not_liv_2 <- Reduce(intersect, list(pat_ema_oment_1_not_liv_2_long$location, pat_ema_oment_2_not_liv_2_long$location)) #18


#subset out those for next branch
pat_ema_oment_1_next_length <- pat_ema_oment_1_not_liv_2_long[pat_ema_oment_1_not_liv_2_long$location %!in% pat_ema_not_liv_2, ] #224
pat_ema_oment_2_next_length <- pat_ema_oment_2_not_liv_2_long[pat_ema_oment_2_not_liv_2_long$location %!in% pat_ema_not_liv_2, ] #265

#in plasma?
length(intersect(pat_ema_oment_1_next_length$location, pat_ema_plasma$location)) #15
length(intersect(pat_ema_oment_2_next_length$location, pat_ema_plasma$location)) #57



## regression figures -----
pat_ema_pl_h <- upset_df[upset_df$Heart == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_h <- rownames(pat_ema_pl_h)

pat_ema_heart_new <- pat_ema_heart[pat_ema_heart$location %in% pat_ema_pl_h, ]
pat_ema_heart_new <- pat_ema_heart_new[, c('location', 'AF')]
pat_ema_heart_new$AF[c(4, 7, 8, 9, 12)] <- c(0.279, 0.302, 0.241, 0.256, 0.164)

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_h, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[c(2, 3, 4, 5, 6, 7, 8, 9, 12)] <- c(0.049, 0.026, 0.424, 0.101, 0.664, 0.320, 0.236, 0.258, 0.197)

pat_ema_pl_h_bind <- merge(pat_ema_plasma_new, pat_ema_heart_new, by = 'location', all = TRUE)
pat_ema_pl_h_bind <- pat_ema_pl_h_bind[, -1]
pat_ema_pl_h_bind[, 1] <- as.numeric(pat_ema_pl_h_bind[, 1])
pat_ema_pl_h_bind[, 2] <- as.numeric(pat_ema_pl_h_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_h_bind)

pat_ema_pl_lk <- upset_df[upset_df$L_Kidney == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_lk <- rownames(pat_ema_pl_lk)

pat_ema_l_kidney_new <- pat_ema_l_kidney[pat_ema_l_kidney$location %in% pat_ema_pl_lk, ]
pat_ema_l_kidney_new <- pat_ema_l_kidney_new[, c('location', 'AF')]
pat_ema_l_kidney_new$AF[c(4, 12)] <- c(0.513, 0.063)

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_lk, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[c(4, 5, 6, 7, 11, 12)] <- c(0.424, 0.320, 0.236, 0.258, 0.197, 0.197)

pat_ema_pl_lk_bind <- merge(pat_ema_plasma_new, pat_ema_l_kidney_new, by = 'location', all = TRUE)
pat_ema_pl_lk_bind <- pat_ema_pl_lk_bind[, -1]
pat_ema_pl_lk_bind[, 1] <- as.numeric(pat_ema_pl_lk_bind[, 1])
pat_ema_pl_lk_bind[, 2] <- as.numeric(pat_ema_pl_lk_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_lk_bind)

pat_ema_pl_rk <- upset_df[upset_df$R_Kidney == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_rk <- rownames(pat_ema_pl_rk)

pat_ema_r_kidney_new <- pat_ema_r_kidney[pat_ema_r_kidney$location %in% pat_ema_pl_rk, ]
pat_ema_r_kidney_new <- pat_ema_r_kidney_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_rk, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[c(4, 6)] <- c(0.424, 0.236)

pat_ema_pl_rk_bind <- merge(pat_ema_plasma_new, pat_ema_r_kidney_new, by = 'location', all = TRUE)
pat_ema_pl_rk_bind <- pat_ema_pl_rk_bind[, -1]
pat_ema_pl_rk_bind[, 1] <- as.numeric(pat_ema_pl_rk_bind[, 1])
pat_ema_pl_rk_bind[, 2] <- as.numeric(pat_ema_pl_rk_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_rk_bind)

pat_ema_pl_l1 <- upset_df[upset_df$Liver_1 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_l1 <- rownames(pat_ema_pl_l1)

pat_ema_liver_1_new <- pat_ema_liver_1[pat_ema_liver_1$location %in% pat_ema_pl_l1, ]
pat_ema_liver_1_new <- pat_ema_liver_1_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_l1, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[2] <- 0.258

pat_ema_pl_l1_bind <- merge(pat_ema_plasma_new, pat_ema_liver_1_new, by = 'location', all = TRUE)
pat_ema_pl_l1_bind <- pat_ema_pl_l1_bind[, -1]
pat_ema_pl_l1_bind[, 1] <- as.numeric(pat_ema_pl_l1_bind[, 1])
pat_ema_pl_l1_bind[, 2] <- as.numeric(pat_ema_pl_l1_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_l1_bind)

pat_ema_l2 <- upset_df[upset_df$Liver_2 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_l2 <- rownames(pat_ema_l2)

pat_ema_liver_2_new <- pat_ema_liver_2[pat_ema_liver_2$location %in% pat_ema_l2, ]
pat_ema_liver_2_new <- pat_ema_liver_2_new[, c('location', 'AF')]
pat_ema_liver_2_new$AF[c(6, 9)] <- c(0.332, 0.062)

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_l2, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[c(3, 5, 6, 8, 9)] <- c(0.424, 0.664, 0.236, 0.197, 0.197)

pat_ema_l2_bind <- merge(pat_ema_plasma_new, pat_ema_liver_2_new, by = 'location', all = TRUE)
pat_ema_l2_bind <- pat_ema_l2_bind[, -1]
pat_ema_l2_bind[, 1] <- as.numeric(pat_ema_l2_bind[, 1])
pat_ema_l2_bind[, 2] <- as.numeric(pat_ema_l2_bind[, 2])
lm(AF.x ~ ., data = pat_ema_l2_bind)

pat_ema_pl_om1 <- upset_df[upset_df$Oment_1 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_om1 <- rownames(pat_ema_pl_om1)

pat_ema_oment_1_new <- pat_ema_oment_1[pat_ema_oment_1$location %in% pat_ema_pl_om1, ]
pat_ema_oment_1_new <- pat_ema_oment_1_new[, c('location', 'AF')]
pat_ema_oment_1_new$AF[c(6, 7, 10)] <- c(0.283, 0.301, 0.088)

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_om1, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[c(3, 4, 5, 6, 7, 8, 10)] <- c(0.424, 0.664, 0.320, 0.236, 0.258, 0.094, 0.197)

pat_ema_pl_om1_bind <- merge(pat_ema_plasma_new, pat_ema_oment_1_new, by = 'location', all = TRUE)
pat_ema_pl_om1_bind <- pat_ema_pl_om1_bind[, -1]
pat_ema_pl_om1_bind[, 1] <- as.numeric(pat_ema_pl_om1_bind[, 1])
pat_ema_pl_om1_bind[, 2] <- as.numeric(pat_ema_pl_om1_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_om1_bind)

pat_ema_pl_om2 <- upset_df[upset_df$Oment_2 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_om2 <- rownames(pat_ema_pl_om2)

pat_ema_oment_2_new <- pat_ema_oment_2[pat_ema_oment_2$location %in% pat_ema_pl_om2, ]
pat_ema_oment_2_new <- pat_ema_oment_2_new[, c('location', 'AF')]
pat_ema_oment_2_new$AF[c(3, 4, 6, 8, 10, 11, 16, 17)] <- c(0.189, 0.507, 0.029, 0.181, 0.345, 0.260, 0.147, 0.110)

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_om2, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]
pat_ema_plasma_new$AF[c(2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 16, 17)] <- c(0.049, 0.424, 0.101, 0.664, 0.252, 0.320, 0.236, 0.258, 0.364, 0.094, 0.197, 0.197)

pat_ema_pl_om2_bind <- merge(pat_ema_plasma_new, pat_ema_oment_2_new, by = 'location', all = TRUE)
pat_ema_pl_om2_bind <- pat_ema_pl_om2_bind[, -1]
pat_ema_pl_om2_bind[, 1] <- as.numeric(pat_ema_pl_om2_bind[, 1])
pat_ema_pl_om2_bind[, 2] <- as.numeric(pat_ema_pl_om2_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_om2_bind)
# all kidney
pat_ema_all_kidney <- rbind(pat_ema_pl_lk_bind, pat_ema_pl_rk_bind)
lm(pat_ema_all_kidney$AF.x ~ pat_ema_all_kidney$AF.y)

pat_ema_all_liver <- rbind(pat_ema_pl_l1_bind, pat_ema_l2_bind)
lm(pat_ema_all_liver$AF.x ~ pat_ema_all_liver$AF.y)

#need to add pat 9 to oment plus ovary and lymph THIS CAN ONLY BE DONE WITH PAT 9 UPSET DF IN ENVIRON!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_9_pl_om <- upset_df[upset_df$Omental_Met == 1 & upset_df$Plasma == 1, ] #12
pat_9_pl_om <- rownames(pat_9_pl_om)

pat_9_oment_new <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_om, ]
pat_9_oment_new <- pat_9_oment_new[, c('location', 'AF')]
pat_9_oment_new$AF[c(11, 13, 16, 19)] <- c(0.343, 0.046, 0.203, 0.295)

pat_9_plasma_new <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_om, ]
pat_9_plasma_new <- pat_9_plasma_new[, c('location', 'AF')]
pat_9_plasma_new$AF[c(3, 11, 13, 14, 15, 16, 19, 20)] <- c(0.020, 0.333, 0.057, 0.230, 0.332, 0.260, 0.209, 0.177)

pat_9_pl_om_bind <- merge(pat_9_plasma_new, pat_9_oment_new, by = 'location', all = TRUE)
pat_9_pl_om_bind <- pat_9_pl_om_bind[, -1]
pat_9_pl_om_bind[, 1] <- as.numeric(pat_9_pl_om_bind[, 1])
pat_9_pl_om_bind[, 2] <- as.numeric(pat_9_pl_om_bind[, 2])
lm(AF.x ~ ., data = pat_9_pl_om_bind)

pat_ema_all_oment <- rbind(pat_ema_pl_om1_bind, pat_ema_pl_om2_bind)
pat_ema_all_oment <- rbind(pat_ema_all_oment, pat_9_pl_om_bind)
lm(pat_ema_all_oment$AF.x ~ pat_ema_all_oment$AF.y)

pat_9_pl_ov <- upset_df[upset_df$Ovary_Met == 1 & upset_df$Plasma == 1, ] #12
pat_9_pl_ov <- rownames(pat_9_pl_ov)

pat_9_ovary_new <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov, ]
pat_9_ovary_new <- pat_9_ovary_new[, c('location', 'AF')]
pat_9_ovary_new$AF[c(11, 14, 15, 16, 20, 21)] <- c(0.330, 0.143, 0.465, 0.246, 0.209, 0.165)

pat_9_plasma_new <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov, ]
pat_9_plasma_new <- pat_9_plasma_new[, c('location', 'AF')]
pat_9_plasma_new$AF[c(3, 11, 13, 14, 15, 16, 20, 21)] <- c(0.020, 0.333, 0.230, 0.332, 0.260, 0.281, 0.209, 0.177)

pat_9_pl_ov_bind <- merge(pat_9_plasma_new, pat_9_ovary_new, by = 'location', all = TRUE)
pat_9_pl_ov_bind <- pat_9_pl_ov_bind[, -1]
pat_9_pl_ov_bind[, 1] <- as.numeric(pat_9_pl_ov_bind[, 1])
pat_9_pl_ov_bind[, 2] <- as.numeric(pat_9_pl_ov_bind[, 2])
lm(AF.x ~ ., data = pat_9_pl_ov_bind)

pat_9_pl_ln <- upset_df[upset_df$Lymph_Met == 1 & upset_df$Plasma == 1, ] #12
pat_9_pl_ln <- rownames(pat_9_pl_ln)

pat_9_ln_new <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ln, ]
pat_9_ln_new <- pat_9_ln_new[, c('location', 'AF')]

pat_9_plasma_new <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ln, ]
pat_9_plasma_new <- pat_9_plasma_new[, c('location', 'AF')]

pat_9_pl_ln_bind <- merge(pat_9_plasma_new, pat_9_ln_new, by = 'location', all = TRUE)
pat_9_pl_ln_bind <- pat_9_pl_ln_bind[, -1]
pat_9_pl_ln_bind[, 1] <- as.numeric(pat_9_pl_ln_bind[, 1])
pat_9_pl_ln_bind[, 2] <- as.numeric(pat_9_pl_ln_bind[, 2])
lm(AF.x ~ ., data = pat_9_pl_ln_bind)


plot(1, type="n", xlab="Mutant Allele Frequency in Tumor", ylab="Mutant Allele Frequency in Plasma", xlim=c(0, 1.0), ylim=c(0, 1.0))
abline(a = 0.06224, b = 1.01361, col = 'red') #heart
abline(a = 0.1177, b = 0.4953, col = 'blue') #all kidney
abline(a = 0.1152, b = 0.8010, col = 'green') #lymph
abline(a = 0.05954, b = 0.66369, col = 'orange') # all liver
#abline(a = 0.07176, b = 0.61002, col = 'purple') # liver 2
abline(a = 0.04582, b = 0.88119, col = 'purple') # all oment
abline(a = 0.1029, b = 0.5760, col = 'dodgerblue') #ovary
legend(x = 0.0, y = 0.9, legend = c('Heart', 'Kidney', 'Lymph', 'Liver', 'Omental', 'Ovary'), col = c('red', 'blue', 'green', 'orange', 'purple', 'dodgerblue'), lty = 1, bty = 'n')

## correlation plots ---
tumor_size <- c(1.8, 2.0, 8.0, 11.0, 10.0)
mean_plasma <- c(0.3102, 0.2987, 0.2013, 0.2285, 0.7834)
tumor_graph_df <- data.frame(tumor_size, mean_plasma)
plot(tumor_graph_df$tumor_size, tumor_graph_df$mean_plasma, ylim = c(0, max(tumor_graph_df$mean_plasma)), xlim = c(0, max(tumor_graph_df$tumor_size)), 
     col = c('blue', 'orange', 'dodgerblue', 'purple', 'green'), pch = 16, cex = 1.2, xlab = expression(paste('Tumor Size (cm' ^ 2, ')')), ylab = 'Mean MAF in Plasma')
fit <- lm(tumor_graph_df$mean_plasma ~ tumor_graph_df$tumor_size)
abline(a = 0.25429, b = 0.01514)
summary(fit)
legend(x = 0.0, y = 0.75, legend = c('Kidney', 'Lymph', 'Liver', 'Omental', 'Ovary'), pch = 16, cex = 1.2, col = c('blue', 'green', 'orange', 'purple', 'dodgerblue'), bty = 'n')
legend(x = 6.5, y = 0.6, legend = expression(paste('R'^2, '= 0.08')), bty = 'n')

perfusion <- c(0.25, 0.2, 0.0, 0.03, 0.0)
perf_graph_df <- data.frame(mean_plasma, perfusion)
plot(perf_graph_df$perfusion, perf_graph_df$mean_plasma)
fit2 <- lm(perf_graph_df$mean_plasma ~ perf_graph_df$perfusion)
fit2
abline(a = 0.4037, b = -0.5221)
summary(fit2)
