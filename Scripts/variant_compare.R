# clean mutation analyses

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
  
  #all calls with trusight genes
  mutect_trusight <- mutect_shorter[mutect_shorter$gene_name %in% trusight, ] 
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
pat_10_liver_1_ug <- read.delim('Data/pat_10_liver_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #12953
pat_10_liver_1_varscan <- read.delim('Data/pat_10_liver_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #6056
pat_10_liver_1_freebayes <- read.delim('Data/pat_10_liver_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #112330

pat_10_liver_1_mutect <- read.delim('Data/pat_10_liver_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #188803

# dummy variable for matching
pat_10_liver_1_ug$location <- paste0(pat_10_liver_1_ug$CHROM, ':', pat_10_liver_1_ug$POS)
pat_10_liver_1_varscan$location <- paste0(pat_10_liver_1_varscan$CHROM, ':', pat_10_liver_1_varscan$POS)
pat_10_liver_1_freebayes$location <- paste0(pat_10_liver_1_freebayes$CHROM, ':', pat_10_liver_1_freebayes$POS)

pat_10_liver_1_mutect$location <- paste0(pat_10_liver_1_mutect$CHROM, ':', pat_10_liver_1_mutect$POS)

#pairwise intersections for the 3
pat_10_liver_1_ug_fb <- intersect(pat_10_liver_1_ug$location, pat_10_liver_1_freebayes$location) #11581
pat_10_liver_1_ug_vs <- intersect(pat_10_liver_1_ug$location, pat_10_liver_1_varscan$location) #2733
pat_10_liver_1_vs_fb <- intersect(pat_10_liver_1_varscan$location, pat_10_liver_1_freebayes$location) #2673

#intersect b/w mutect calls and two out of the 3 others
pat_10_liver_1_mutect_ug_fb <- intersect(pat_10_liver_1_mutect$location, pat_10_liver_1_ug_fb) #6526
pat_10_liver_1_mutect_ug_vs <- intersect(pat_10_liver_1_mutect$location, pat_10_liver_1_ug_vs) #658
pat_10_liver_1_mutect_vs_fb <- intersect(pat_10_liver_1_mutect$location, pat_10_liver_1_vs_fb) #602

#combine them
pat_10_liver_1_all <- c(pat_10_liver_1_mutect_ug_fb, pat_10_liver_1_mutect_ug_vs, pat_10_liver_1_mutect_vs_fb)

#subset mutect calls
pat_10_liver_1_mutect_short <- pat_10_liver_1_mutect[pat_10_liver_1_mutect$location %in% pat_10_liver_1_all, ] #6680

#trim to trusight gene calls
pat_10_liver_1_mutect_trusight <- trusight_from_mutect(pat_10_liver_1_mutect_short) #17

## liver 2a----
# import data
pat_10_liver_2a_ug <- read.delim('Data/pat_10_liver_2a_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #9193
pat_10_liver_2a_varscan <- read.delim('Data/pat_10_liver_2a_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3947
pat_10_liver_2a_freebayes <- read.delim('Data/pat_10_liver_2a_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #85821

pat_10_liver_2a_mutect <- read.delim('Data/pat_10_liver_2a_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #124767

# dummy variable for matching
pat_10_liver_2a_ug$location <- paste0(pat_10_liver_2a_ug$CHROM, ':', pat_10_liver_2a_ug$POS)
pat_10_liver_2a_varscan$location <- paste0(pat_10_liver_2a_varscan$CHROM, ':', pat_10_liver_2a_varscan$POS)
pat_10_liver_2a_freebayes$location <- paste0(pat_10_liver_2a_freebayes$CHROM, ':', pat_10_liver_2a_freebayes$POS)

pat_10_liver_2a_mutect$location <- paste0(pat_10_liver_2a_mutect$CHROM, ':', pat_10_liver_2a_mutect$POS)

#pairwise intersections for the 3
pat_10_liver_2a_ug_fb <- intersect(pat_10_liver_2a_ug$location, pat_10_liver_2a_freebayes$location) #8146
pat_10_liver_2a_ug_vs <- intersect(pat_10_liver_2a_ug$location, pat_10_liver_2a_varscan$location) #2182
pat_10_liver_2a_vs_fb <- intersect(pat_10_liver_2a_varscan$location, pat_10_liver_2a_freebayes$location) #2073

#intersect b/w mutect calls and two out of the 3 others
pat_10_liver_2a_mutect_ug_fb <- intersect(pat_10_liver_2a_mutect$location, pat_10_liver_2a_ug_fb) #3639
pat_10_liver_2a_mutect_ug_vs <- intersect(pat_10_liver_2a_mutect$location, pat_10_liver_2a_ug_vs) #393
pat_10_liver_2a_mutect_vs_fb <- intersect(pat_10_liver_2a_mutect$location, pat_10_liver_2a_vs_fb) #354

#combine them
pat_10_liver_2a_all <- c(pat_10_liver_2a_mutect_ug_fb, pat_10_liver_2a_mutect_ug_vs, pat_10_liver_2a_mutect_vs_fb)

#subset mutect calls
pat_10_liver_2a_mutect_short <- pat_10_liver_2a_mutect[pat_10_liver_2a_mutect$location %in% pat_10_liver_2a_all, ] #3728

#trim to trusight gene calls
pat_10_liver_2a_mutect_trusight <- trusight_from_mutect(pat_10_liver_2a_mutect_short) #10

## liver 5----
# import data
pat_10_liver_5_ug <- read.delim('Data/pat_10_liver_5_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #12980
pat_10_liver_5_varscan <- read.delim('Data/pat_10_liver_5_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4821
pat_10_liver_5_freebayes <- read.delim('Data/pat_10_liver_5_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #118457

pat_10_liver_5_mutect <- read.delim('Data/pat_10_liver_5_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #193162

# dummy variable for matching
pat_10_liver_5_ug$location <- paste0(pat_10_liver_5_ug$CHROM, ':', pat_10_liver_5_ug$POS)
pat_10_liver_5_varscan$location <- paste0(pat_10_liver_5_varscan$CHROM, ':', pat_10_liver_5_varscan$POS)
pat_10_liver_5_freebayes$location <- paste0(pat_10_liver_5_freebayes$CHROM, ':', pat_10_liver_5_freebayes$POS)

pat_10_liver_5_mutect$location <- paste0(pat_10_liver_5_mutect$CHROM, ':', pat_10_liver_5_mutect$POS)

#pairwise intersections for the 3
pat_10_liver_5_ug_fb <- intersect(pat_10_liver_5_ug$location, pat_10_liver_5_freebayes$location) #11857
pat_10_liver_5_ug_vs <- intersect(pat_10_liver_5_ug$location, pat_10_liver_5_varscan$location) #2516
pat_10_liver_5_vs_fb <- intersect(pat_10_liver_5_varscan$location, pat_10_liver_5_freebayes$location) #2436

#intersect b/w mutect calls and two out of the 3 others
pat_10_liver_5_mutect_ug_fb <- intersect(pat_10_liver_5_mutect$location, pat_10_liver_5_ug_fb) #6711
pat_10_liver_5_mutect_ug_vs <- intersect(pat_10_liver_5_mutect$location, pat_10_liver_5_ug_vs) #528
pat_10_liver_5_mutect_vs_fb <- intersect(pat_10_liver_5_mutect$location, pat_10_liver_5_vs_fb) #495

#combine them
pat_10_liver_5_all <- c(pat_10_liver_5_mutect_ug_fb, pat_10_liver_5_mutect_ug_vs, pat_10_liver_5_mutect_vs_fb)

#subset mutect calls
pat_10_liver_5_mutect_short <- pat_10_liver_5_mutect[pat_10_liver_5_mutect$location %in% pat_10_liver_5_all, ] #6832

#trim to trusight gene calls
pat_10_liver_5_mutect_trusight <- trusight_from_mutect(pat_10_liver_5_mutect_short) #17

## plasma----
# import data
pat_10_plasma_ug <- read.delim('Data/pat_10_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #45442
pat_10_plasma_varscan <- read.delim('Data/pat_10_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #47419
pat_10_plasma_freebayes <- read.delim('Data/pat_10_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #356698

pat_10_plasma_mutect <- read.delim('Data/pat_10_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #575196

# dummy variable for matching
pat_10_plasma_ug$location <- paste0(pat_10_plasma_ug$CHROM, ':', pat_10_plasma_ug$POS)
pat_10_plasma_varscan$location <- paste0(pat_10_plasma_varscan$CHROM, ':', pat_10_plasma_varscan$POS)
pat_10_plasma_freebayes$location <- paste0(pat_10_plasma_freebayes$CHROM, ':', pat_10_plasma_freebayes$POS)

pat_10_plasma_mutect$location <- paste0(pat_10_plasma_mutect$CHROM, ':', pat_10_plasma_mutect$POS)

#pairwise intersections for the 3
pat_10_plasma_ug_fb <- intersect(pat_10_plasma_ug$location, pat_10_plasma_freebayes$location) #34782
pat_10_plasma_ug_vs <- intersect(pat_10_plasma_ug$location, pat_10_plasma_varscan$location) #5624
pat_10_plasma_vs_fb <- intersect(pat_10_plasma_varscan$location, pat_10_plasma_freebayes$location) #5250

#intersect b/w mutect calls and two out of the 3 others
pat_10_plasma_mutect_ug_fb <- intersect(pat_10_plasma_mutect$location, pat_10_plasma_ug_fb) #25324
pat_10_plasma_mutect_ug_vs <- intersect(pat_10_plasma_mutect$location, pat_10_plasma_ug_vs) #3249
pat_10_plasma_mutect_vs_fb <- intersect(pat_10_plasma_mutect$location, pat_10_plasma_vs_fb) #2915

#combine them
pat_10_plasma_all <- c(pat_10_plasma_mutect_ug_fb, pat_10_plasma_mutect_ug_vs, pat_10_plasma_mutect_vs_fb)

#subset mutect calls
pat_10_plasma_mutect_short <- pat_10_plasma_mutect[pat_10_plasma_mutect$location %in% pat_10_plasma_all, ] #25932

#trim to trusight gene calls
pat_10_plasma_mutect_trusight <- trusight_from_mutect(pat_10_plasma_mutect_short) #85

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 liver mets
pat_10_met_pool <- c(pat_10_liver_1_mutect_trusight$location, pat_10_liver_2a_mutect_trusight$location, pat_10_liver_5_mutect_trusight$location)
pat_10_met_pool <- unique(pat_10_met_pool) #23

pat_10_plasma_link <- pat_10_plasma_mutect_trusight[pat_10_plasma_mutect_trusight$location %in% pat_10_met_pool, ] #15 mutations in 9 genes









## PATIENT 9 ----
## lymph node----
# import data
pat_9_ln_ug <- read.delim('Data/pat_9_ln_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4463
pat_9_ln_varscan <- read.delim('Data/pat_9_ln_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #121
pat_9_ln_freebayes <- read.delim('Data/pat_9_ln_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3975

pat_9_ln_mutect <- read.delim('Data/pat_9_ln_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #2132

# dummy variable for matching
pat_9_ln_ug$location <- paste0(pat_9_ln_ug$CHROM, ':', pat_9_ln_ug$POS)
pat_9_ln_varscan$location <- paste0(pat_9_ln_varscan$CHROM, ':', pat_9_ln_varscan$POS)
pat_9_ln_freebayes$location <- paste0(pat_9_ln_freebayes$CHROM, ':', pat_9_ln_freebayes$POS)

pat_9_ln_mutect$location <- paste0(pat_9_ln_mutect$CHROM, ':', pat_9_ln_mutect$POS)

#pairwise intersections for the 3
pat_9_ln_ug_fb <- intersect(pat_9_ln_ug$location, pat_9_ln_freebayes$location) #1939
pat_9_ln_ug_vs <- intersect(pat_9_ln_ug$location, pat_9_ln_varscan$location) #66
pat_9_ln_vs_fb <- intersect(pat_9_ln_varscan$location, pat_9_ln_freebayes$location) #64

#intersect b/w mutect calls and two out of the 3 others
pat_9_ln_mutect_ug_fb <- intersect(pat_9_ln_mutect$location, pat_9_ln_ug_fb) #473
pat_9_ln_mutect_ug_vs <- intersect(pat_9_ln_mutect$location, pat_9_ln_ug_vs) #24
pat_9_ln_mutect_vs_fb <- intersect(pat_9_ln_mutect$location, pat_9_ln_vs_fb) #22

#combine them
pat_9_ln_all <- c(pat_9_ln_mutect_ug_fb, pat_9_ln_mutect_ug_vs, pat_9_ln_mutect_vs_fb)

#subset mutect calls
pat_9_ln_mutect_short <- pat_9_ln_mutect[pat_9_ln_mutect$location %in% pat_9_ln_all, ] #475

#trim to trusight gene calls
pat_9_ln_mutect_trusight <- trusight_from_mutect(pat_9_ln_mutect_short) #10

## oment----
# import data
pat_9_oment_ug <- read.delim('Data/pat_9_oment_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11359
pat_9_oment_varscan <- read.delim('Data/pat_9_oment_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4839
pat_9_oment_freebayes <- read.delim('Data/pat_9_oment_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #88156

pat_9_oment_mutect <- read.delim('Data/pat_9_oment_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #114997

# dummy variable for matching
pat_9_oment_ug$location <- paste0(pat_9_oment_ug$CHROM, ':', pat_9_oment_ug$POS)
pat_9_oment_varscan$location <- paste0(pat_9_oment_varscan$CHROM, ':', pat_9_oment_varscan$POS)
pat_9_oment_freebayes$location <- paste0(pat_9_oment_freebayes$CHROM, ':', pat_9_oment_freebayes$POS)

pat_9_oment_mutect$location <- paste0(pat_9_oment_mutect$CHROM, ':', pat_9_oment_mutect$POS)

#pairwise intersections for the 3
pat_9_oment_ug_fb <- intersect(pat_9_oment_ug$location, pat_9_oment_freebayes$location) #8155
pat_9_oment_ug_vs <- intersect(pat_9_oment_ug$location, pat_9_oment_varscan$location) #1761
pat_9_oment_vs_fb <- intersect(pat_9_oment_varscan$location, pat_9_oment_freebayes$location) #1643

#intersect b/w mutect calls and two out of the 3 others
pat_9_oment_mutect_ug_fb <- intersect(pat_9_oment_mutect$location, pat_9_oment_ug_fb) #3432
pat_9_oment_mutect_ug_vs <- intersect(pat_9_oment_mutect$location, pat_9_oment_ug_vs) #575
pat_9_oment_mutect_vs_fb <- intersect(pat_9_oment_mutect$location, pat_9_oment_vs_fb) #507

#combine them
pat_9_oment_all <- c(pat_9_oment_mutect_ug_fb, pat_9_oment_mutect_ug_vs, pat_9_oment_mutect_vs_fb)

#subset mutect calls
pat_9_oment_mutect_short <- pat_9_oment_mutect[pat_9_oment_mutect$location %in% pat_9_oment_all, ] #3536

#trim to trusight gene calls
pat_9_oment_mutect_trusight <- trusight_from_mutect(pat_9_oment_mutect_short) #40

## ovary----
# import data
pat_9_ovary_ug <- read.delim('Data/pat_9_ovary_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #8037
pat_9_ovary_varscan <- read.delim('Data/pat_9_ovary_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #7513
pat_9_ovary_freebayes <- read.delim('Data/pat_9_ovary_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #96219

pat_9_ovary_mutect <- read.delim('Data/pat_9_ovary_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #92251

# dummy variable for matching
pat_9_ovary_ug$location <- paste0(pat_9_ovary_ug$CHROM, ':', pat_9_ovary_ug$POS)
pat_9_ovary_varscan$location <- paste0(pat_9_ovary_varscan$CHROM, ':', pat_9_ovary_varscan$POS)
pat_9_ovary_freebayes$location <- paste0(pat_9_ovary_freebayes$CHROM, ':', pat_9_ovary_freebayes$POS)

pat_9_ovary_mutect$location <- paste0(pat_9_ovary_mutect$CHROM, ':', pat_9_ovary_mutect$POS)

#pairwise intersections for the 3
pat_9_ovary_ug_fb <- intersect(pat_9_ovary_ug$location, pat_9_ovary_freebayes$location) #6580
pat_9_ovary_ug_vs <- intersect(pat_9_ovary_ug$location, pat_9_ovary_varscan$location) #1391
pat_9_ovary_vs_fb <- intersect(pat_9_ovary_varscan$location, pat_9_ovary_freebayes$location) #1300

#intersect b/w mutect calls and two out of the 3 others
pat_9_ovary_mutect_ug_fb <- intersect(pat_9_ovary_mutect$location, pat_9_ovary_ug_fb) #2981
pat_9_ovary_mutect_ug_vs <- intersect(pat_9_ovary_mutect$location, pat_9_ovary_ug_vs) #519
pat_9_ovary_mutect_vs_fb <- intersect(pat_9_ovary_mutect$location, pat_9_ovary_vs_fb) #465

#combine them
pat_9_ovary_all <- c(pat_9_ovary_mutect_ug_fb, pat_9_ovary_mutect_ug_vs, pat_9_ovary_mutect_vs_fb)

#subset mutect calls
pat_9_ovary_mutect_short <- pat_9_ovary_mutect[pat_9_ovary_mutect$location %in% pat_9_ovary_all, ] #3055

#trim to trusight gene calls
pat_9_ovary_mutect_trusight <- trusight_from_mutect(pat_9_ovary_mutect_short) #20

## plasma----
# import data
pat_9_plasma_ug <- read.delim('Data/pat_9_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #59407
pat_9_plasma_varscan <- read.delim('Data/pat_9_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #20319
pat_9_plasma_freebayes <- read.delim('Data/pat_9_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #427734

pat_9_plasma_mutect <- read.delim('Data/pat_9_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #782680

# dummy variable for matching
pat_9_plasma_ug$location <- paste0(pat_9_plasma_ug$CHROM, ':', pat_9_plasma_ug$POS)
pat_9_plasma_varscan$location <- paste0(pat_9_plasma_varscan$CHROM, ':', pat_9_plasma_varscan$POS)
pat_9_plasma_freebayes$location <- paste0(pat_9_plasma_freebayes$CHROM, ':', pat_9_plasma_freebayes$POS)

pat_9_plasma_mutect$location <- paste0(pat_9_plasma_mutect$CHROM, ':', pat_9_plasma_mutect$POS)

#pairwise intersections for the 3
pat_9_plasma_ug_fb <- intersect(pat_9_plasma_ug$location, pat_9_plasma_freebayes$location) #52862
pat_9_plasma_ug_vs <- intersect(pat_9_plasma_ug$location, pat_9_plasma_varscan$location) #6420
pat_9_plasma_vs_fb <- intersect(pat_9_plasma_varscan$location, pat_9_plasma_freebayes$location) #6347

#intersect b/w mutect calls and two out of the 3 others
pat_9_plasma_mutect_ug_fb <- intersect(pat_9_plasma_mutect$location, pat_9_plasma_ug_fb) #40548
pat_9_plasma_mutect_ug_vs <- intersect(pat_9_plasma_mutect$location, pat_9_plasma_ug_vs) #2990
pat_9_plasma_mutect_vs_fb <- intersect(pat_9_plasma_mutect$location, pat_9_plasma_vs_fb) #2770

#combine them
pat_9_plasma_all <- c(pat_9_plasma_mutect_ug_fb, pat_9_plasma_mutect_ug_vs, pat_9_plasma_mutect_vs_fb)

#subset mutect calls
pat_9_plasma_mutect_short <- pat_9_plasma_mutect[pat_9_plasma_mutect$location %in% pat_9_plasma_all, ] #40970

#trim to trusight gene calls
pat_9_plasma_mutect_trusight <- trusight_from_mutect(pat_9_plasma_mutect_short) #59

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 mets
pat_9_met_pool <- c(pat_9_ln_mutect_trusight$location, pat_9_oment_mutect_trusight$location, pat_9_ovary_mutect_trusight$location)
pat_9_met_pool <- unique(pat_9_met_pool) #58

pat_9_plasma_link <- pat_9_plasma_mutect_trusight[pat_9_plasma_mutect_trusight$location %in% pat_9_met_pool, ] #32 mutations in 11 genes






## PATIENT 8 ----
## axillary----
# import data
pat_8_axillary_ug <- read.delim('Data/pat_8_axillary_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #7582
pat_8_axillary_varscan <- read.delim('Data/pat_8_axillary_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3213
pat_8_axillary_freebayes <- read.delim('Data/pat_8_axillary_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #54045

pat_8_axillary_mutect <- read.delim('Data/pat_8_axillary_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #68044

# dummy variable for matching
pat_8_axillary_ug$location <- paste0(pat_8_axillary_ug$CHROM, ':', pat_8_axillary_ug$POS)
pat_8_axillary_varscan$location <- paste0(pat_8_axillary_varscan$CHROM, ':', pat_8_axillary_varscan$POS)
pat_8_axillary_freebayes$location <- paste0(pat_8_axillary_freebayes$CHROM, ':', pat_8_axillary_freebayes$POS)

pat_8_axillary_mutect$location <- paste0(pat_8_axillary_mutect$CHROM, ':', pat_8_axillary_mutect$POS)

#pairwise intersections for the 3
pat_8_axillary_ug_fb <- intersect(pat_8_axillary_ug$location, pat_8_axillary_freebayes$location) #6239
pat_8_axillary_ug_vs <- intersect(pat_8_axillary_ug$location, pat_8_axillary_varscan$location) #1773
pat_8_axillary_vs_fb <- intersect(pat_8_axillary_varscan$location, pat_8_axillary_freebayes$location) #22

#intersect b/w mutect calls and two out of the 3 others
pat_8_axillary_mutect_ug_fb <- intersect(pat_8_axillary_mutect$location, pat_8_axillary_ug_fb) #1937
pat_8_axillary_mutect_ug_vs <- intersect(pat_8_axillary_mutect$location, pat_8_axillary_ug_vs) #271
pat_8_axillary_mutect_vs_fb <- intersect(pat_8_axillary_mutect$location, pat_8_axillary_vs_fb) #235

#combine them
pat_8_axillary_all <- c(pat_8_axillary_mutect_ug_fb, pat_8_axillary_mutect_ug_vs, pat_8_axillary_mutect_vs_fb)

#subset mutect calls
pat_8_axillary_mutect_short <- pat_8_axillary_mutect[pat_8_axillary_mutect$location %in% pat_8_axillary_all, ] #1987

#trim to trusight gene calls
pat_8_axillary_mutect_trusight <- trusight_from_mutect(pat_8_axillary_mutect_short) #12

## breast 1----
# import data
pat_8_breast_1_ug <- read.delim('Data/pat_8_breast_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11948
pat_8_breast_1_varscan <- read.delim('Data/pat_8_breast_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3669
pat_8_breast_1_freebayes <- read.delim('Data/pat_8_breast_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #74674

pat_8_breast_1_mutect <- read.delim('Data/pat_8_breast_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #76214

# dummy variable for matching
pat_8_breast_1_ug$location <- paste0(pat_8_breast_1_ug$CHROM, ':', pat_8_breast_1_ug$POS)
pat_8_breast_1_varscan$location <- paste0(pat_8_breast_1_varscan$CHROM, ':', pat_8_breast_1_varscan$POS)
pat_8_breast_1_freebayes$location <- paste0(pat_8_breast_1_freebayes$CHROM, ':', pat_8_breast_1_freebayes$POS)

pat_8_breast_1_mutect$location <- paste0(pat_8_breast_1_mutect$CHROM, ':', pat_8_breast_1_mutect$POS)

#pairwise intersections for the 3
pat_8_breast_1_ug_fb <- intersect(pat_8_breast_1_ug$location, pat_8_breast_1_freebayes$location) #8351
pat_8_breast_1_ug_vs <- intersect(pat_8_breast_1_ug$location, pat_8_breast_1_varscan$location) #1759
pat_8_breast_1_vs_fb <- intersect(pat_8_breast_1_varscan$location, pat_8_breast_1_freebayes$location) #1774

#intersect b/w mutect calls and two out of the 3 others
pat_8_breast_1_mutect_ug_fb <- intersect(pat_8_breast_1_mutect$location, pat_8_breast_1_ug_fb) #3080
pat_8_breast_1_mutect_ug_vs <- intersect(pat_8_breast_1_mutect$location, pat_8_breast_1_ug_vs) #502
pat_8_breast_1_mutect_vs_fb <- intersect(pat_8_breast_1_mutect$location, pat_8_breast_1_vs_fb) #470

#combine them
pat_8_breast_1_all <- c(pat_8_breast_1_mutect_ug_fb, pat_8_breast_1_mutect_ug_vs, pat_8_breast_1_mutect_vs_fb)

#subset mutect calls
pat_8_breast_1_mutect_short <- pat_8_breast_1_mutect[pat_8_breast_1_mutect$location %in% pat_8_breast_1_all, ] #3136

#trim to trusight gene calls
pat_8_breast_1_mutect_trusight <- trusight_from_mutect(pat_8_breast_1_mutect_short) #19

## breast 2----
# import data
pat_8_breast_2_ug <- read.delim('Data/pat_8_breast_2_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #62413
pat_8_breast_2_varscan <- read.delim('Data/pat_8_breast_2_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3669
pat_8_breast_2_freebayes <- read.delim('Data/pat_8_breast_2_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #70717

pat_8_breast_2_mutect <- read.delim('Data/pat_8_breast_2_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #61361

# dummy variable for matching
pat_8_breast_2_ug$location <- paste0(pat_8_breast_2_ug$CHROM, ':', pat_8_breast_2_ug$POS)
pat_8_breast_2_varscan$location <- paste0(pat_8_breast_2_varscan$CHROM, ':', pat_8_breast_2_varscan$POS)
pat_8_breast_2_freebayes$location <- paste0(pat_8_breast_2_freebayes$CHROM, ':', pat_8_breast_2_freebayes$POS)

pat_8_breast_2_mutect$location <- paste0(pat_8_breast_2_mutect$CHROM, ':', pat_8_breast_2_mutect$POS)

#pairwise intersections for the 3
pat_8_breast_2_ug_fb <- intersect(pat_8_breast_2_ug$location, pat_8_breast_2_freebayes$location) #9882
pat_8_breast_2_ug_vs <- intersect(pat_8_breast_2_ug$location, pat_8_breast_2_varscan$location) #1126
pat_8_breast_2_vs_fb <- intersect(pat_8_breast_2_varscan$location, pat_8_breast_2_freebayes$location) #1216

#intersect b/w mutect calls and two out of the 3 others
pat_8_breast_2_mutect_ug_fb <- intersect(pat_8_breast_2_mutect$location, pat_8_breast_2_ug_fb) #2995
pat_8_breast_2_mutect_ug_vs <- intersect(pat_8_breast_2_mutect$location, pat_8_breast_2_ug_vs) #850
pat_8_breast_2_mutect_vs_fb <- intersect(pat_8_breast_2_mutect$location, pat_8_breast_2_vs_fb) #795

#combine them
pat_8_breast_2_all <- c(pat_8_breast_2_mutect_ug_fb, pat_8_breast_2_mutect_ug_vs, pat_8_breast_2_mutect_vs_fb)

#subset mutect calls
pat_8_breast_2_mutect_short <- pat_8_breast_2_mutect[pat_8_breast_2_mutect$location %in% pat_8_breast_2_all, ] #3064

#trim to trusight gene calls
pat_8_breast_2_mutect_trusight <- trusight_from_mutect(pat_8_breast_2_mutect_short) #54

## plasma----
# import data
pat_8_plasma_ug <- read.delim('Data/pat_8_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #5277
pat_8_plasma_varscan <- read.delim('Data/pat_8_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #1837
pat_8_plasma_freebayes <- read.delim('Data/pat_8_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3975

pat_8_plasma_mutect <- read.delim('Data/pat_8_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #15390

# dummy variable for matching
pat_8_plasma_ug$location <- paste0(pat_8_plasma_ug$CHROM, ':', pat_8_plasma_ug$POS)
pat_8_plasma_varscan$location <- paste0(pat_8_plasma_varscan$CHROM, ':', pat_8_plasma_varscan$POS)
pat_8_plasma_freebayes$location <- paste0(pat_8_plasma_freebayes$CHROM, ':', pat_8_plasma_freebayes$POS)

pat_8_plasma_mutect$location <- paste0(pat_8_plasma_mutect$CHROM, ':', pat_8_plasma_mutect$POS)

#pairwise intersections for the 3
pat_8_plasma_ug_fb <- intersect(pat_8_plasma_ug$location, pat_8_plasma_freebayes$location) #1343
pat_8_plasma_ug_vs <- intersect(pat_8_plasma_ug$location, pat_8_plasma_varscan$location) #437
pat_8_plasma_vs_fb <- intersect(pat_8_plasma_varscan$location, pat_8_plasma_freebayes$location) #464

#intersect b/w mutect calls and two out of the 3 others
pat_8_plasma_mutect_ug_fb <- intersect(pat_8_plasma_mutect$location, pat_8_plasma_ug_fb) #3386
pat_8_plasma_mutect_ug_vs <- intersect(pat_8_plasma_mutect$location, pat_8_plasma_ug_vs) #274
pat_8_plasma_mutect_vs_fb <- intersect(pat_8_plasma_mutect$location, pat_8_plasma_vs_fb) #258

#combine them
pat_8_plasma_all <- c(pat_8_plasma_mutect_ug_fb, pat_8_plasma_mutect_ug_vs, pat_8_plasma_mutect_vs_fb)

#subset mutect calls
pat_8_plasma_mutect_short <- pat_8_plasma_mutect[pat_8_plasma_mutect$location %in% pat_8_plasma_all, ] #1367

#trim to trusight gene calls
pat_8_plasma_mutect_trusight <- trusight_from_mutect(pat_8_plasma_mutect_short) #13

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 mets
pat_8_met_pool <- c(pat_8_axillary_mutect_trusight$location, pat_8_breast_1_mutect_trusight$location, pat_8_breast_2_mutect_trusight$location)
pat_8_met_pool <- unique(pat_8_met_pool) #65

pat_8_plasma_link <- pat_8_plasma_mutect_trusight[pat_8_plasma_mutect_trusight$location %in% pat_8_met_pool, ] #10 mutations in 8 genes



## PATIENT 2 ----
## breast 1----
# import data
pat_2_breast_1_ug <- read.delim('Data/pat_2_breast_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #8570
pat_2_breast_1_varscan <- read.delim('Data/pat_2_breast_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3284
pat_2_breast_1_freebayes <- read.delim('Data/pat_2_breast_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #82064

pat_2_breast_1_mutect <- read.delim('Data/pat_2_breast_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #110691

# dummy variable for matching
pat_2_breast_1_ug$location <- paste0(pat_2_breast_1_ug$CHROM, ':', pat_2_breast_1_ug$POS)
pat_2_breast_1_varscan$location <- paste0(pat_2_breast_1_varscan$CHROM, ':', pat_2_breast_1_varscan$POS)
pat_2_breast_1_freebayes$location <- paste0(pat_2_breast_1_freebayes$CHROM, ':', pat_2_breast_1_freebayes$POS)

pat_2_breast_1_mutect$location <- paste0(pat_2_breast_1_mutect$CHROM, ':', pat_2_breast_1_mutect$POS)

#pairwise intersections for the 3
pat_2_breast_1_ug_fb <- intersect(pat_2_breast_1_ug$location, pat_2_breast_1_freebayes$location) #7342
pat_2_breast_1_ug_vs <- intersect(pat_2_breast_1_ug$location, pat_2_breast_1_varscan$location) #1891
pat_2_breast_1_vs_fb <- intersect(pat_2_breast_1_varscan$location, pat_2_breast_1_freebayes$location) #1772

#intersect b/w mutect calls and two out of the 3 others
pat_2_breast_1_mutect_ug_fb <- intersect(pat_2_breast_1_mutect$location, pat_2_breast_1_ug_fb) #3401
pat_2_breast_1_mutect_ug_vs <- intersect(pat_2_breast_1_mutect$location, pat_2_breast_1_ug_vs) #331
pat_2_breast_1_mutect_vs_fb <- intersect(pat_2_breast_1_mutect$location, pat_2_breast_1_vs_fb) #283

#combine them
pat_2_breast_1_all <- c(pat_2_breast_1_mutect_ug_fb, pat_2_breast_1_mutect_ug_vs, pat_2_breast_1_mutect_vs_fb)

#subset mutect calls
pat_2_breast_1_mutect_short <- pat_2_breast_1_mutect[pat_2_breast_1_mutect$location %in% pat_2_breast_1_all, ] #3467

#trim to trusight gene calls
pat_2_breast_1_mutect_trusight <- trusight_from_mutect(pat_2_breast_1_mutect_short) #11

## breast 2----
# import data
pat_2_breast_2_ug <- read.delim('Data/pat_2_breast_2_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #7496
pat_2_breast_2_varscan <- read.delim('Data/pat_2_breast_2_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3369
pat_2_breast_2_freebayes <- read.delim('Data/pat_2_breast_2_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #54925

pat_2_breast_2_mutect <- read.delim('Data/pat_2_breast_2_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #77618

# dummy variable for matching
pat_2_breast_2_ug$location <- paste0(pat_2_breast_2_ug$CHROM, ':', pat_2_breast_2_ug$POS)
pat_2_breast_2_varscan$location <- paste0(pat_2_breast_2_varscan$CHROM, ':', pat_2_breast_2_varscan$POS)
pat_2_breast_2_freebayes$location <- paste0(pat_2_breast_2_freebayes$CHROM, ':', pat_2_breast_2_freebayes$POS)

pat_2_breast_2_mutect$location <- paste0(pat_2_breast_2_mutect$CHROM, ':', pat_2_breast_2_mutect$POS)

#pairwise intersections for the 3
pat_2_breast_2_ug_fb <- intersect(pat_2_breast_2_ug$location, pat_2_breast_2_freebayes$location) #6421
pat_2_breast_2_ug_vs <- intersect(pat_2_breast_2_ug$location, pat_2_breast_2_varscan$location) #1901
pat_2_breast_2_vs_fb <- intersect(pat_2_breast_2_varscan$location, pat_2_breast_2_freebayes$location) #1801

#intersect b/w mutect calls and two out of the 3 others
pat_2_breast_2_mutect_ug_fb <- intersect(pat_2_breast_2_mutect$location, pat_2_breast_2_ug_fb) #2613
pat_2_breast_2_mutect_ug_vs <- intersect(pat_2_breast_2_mutect$location, pat_2_breast_2_ug_vs) #300
pat_2_breast_2_mutect_vs_fb <- intersect(pat_2_breast_2_mutect$location, pat_2_breast_2_vs_fb) #257

#combine them
pat_2_breast_2_all <- c(pat_2_breast_2_mutect_ug_fb, pat_2_breast_2_mutect_ug_vs, pat_2_breast_2_mutect_vs_fb)

#subset mutect calls
pat_2_breast_2_mutect_short <- pat_2_breast_2_mutect[pat_2_breast_2_mutect$location %in% pat_2_breast_2_all, ] #2694

#trim to trusight gene calls
pat_2_breast_2_mutect_trusight <- trusight_from_mutect(pat_2_breast_2_mutect_short) #14

## liver 1----
# import data
pat_2_liver_1_ug <- read.delim('Data/pat_2_liver_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #8201
pat_2_liver_1_varscan <- read.delim('Data/pat_2_liver_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3678
pat_2_liver_1_freebayes <- read.delim('Data/pat_2_liver_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #71556

pat_2_liver_1_mutect <- read.delim('Data/pat_2_liver_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #103999

# dummy variable for matching
pat_2_liver_1_ug$location <- paste0(pat_2_liver_1_ug$CHROM, ':', pat_2_liver_1_ug$POS)
pat_2_liver_1_varscan$location <- paste0(pat_2_liver_1_varscan$CHROM, ':', pat_2_liver_1_varscan$POS)
pat_2_liver_1_freebayes$location <- paste0(pat_2_liver_1_freebayes$CHROM, ':', pat_2_liver_1_freebayes$POS)

pat_2_liver_1_mutect$location <- paste0(pat_2_liver_1_mutect$CHROM, ':', pat_2_liver_1_mutect$POS)

#pairwise intersections for the 3
pat_2_liver_1_ug_fb <- intersect(pat_2_liver_1_ug$location, pat_2_liver_1_freebayes$location) #6822
pat_2_liver_1_ug_vs <- intersect(pat_2_liver_1_ug$location, pat_2_liver_1_varscan$location) #1868
pat_2_liver_1_vs_fb <- intersect(pat_2_liver_1_varscan$location, pat_2_liver_1_freebayes$location) #1812

#intersect b/w mutect calls and two out of the 3 others
pat_2_liver_1_mutect_ug_fb <- intersect(pat_2_liver_1_mutect$location, pat_2_liver_1_ug_fb) #2942
pat_2_liver_1_mutect_ug_vs <- intersect(pat_2_liver_1_mutect$location, pat_2_liver_1_ug_vs) #337
pat_2_liver_1_mutect_vs_fb <- intersect(pat_2_liver_1_mutect$location, pat_2_liver_1_vs_fb) #297

#combine them
pat_2_liver_1_all <- c(pat_2_liver_1_mutect_ug_fb, pat_2_liver_1_mutect_ug_vs, pat_2_liver_1_mutect_vs_fb)

#subset mutect calls
pat_2_liver_1_mutect_short <- pat_2_liver_1_mutect[pat_2_liver_1_mutect$location %in% pat_2_liver_1_all, ] #3016

#trim to trusight gene calls
pat_2_liver_1_mutect_trusight <- trusight_from_mutect(pat_2_liver_1_mutect_short) #20

## liver 2----
# import data
pat_2_liver_2_ug <- read.delim('Data/pat_2_liver_2_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #8296
pat_2_liver_2_varscan <- read.delim('Data/pat_2_liver_2_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3747
pat_2_liver_2_freebayes <- read.delim('Data/pat_2_liver_2_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #57034

pat_2_liver_2_mutect <- read.delim('Data/pat_2_liver_2_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #80398

# dummy variable for matching
pat_2_liver_2_ug$location <- paste0(pat_2_liver_2_ug$CHROM, ':', pat_2_liver_2_ug$POS)
pat_2_liver_2_varscan$location <- paste0(pat_2_liver_2_varscan$CHROM, ':', pat_2_liver_2_varscan$POS)
pat_2_liver_2_freebayes$location <- paste0(pat_2_liver_2_freebayes$CHROM, ':', pat_2_liver_2_freebayes$POS)

pat_2_liver_2_mutect$location <- paste0(pat_2_liver_2_mutect$CHROM, ':', pat_2_liver_2_mutect$POS)

#pairwise intersections for the 3
pat_2_liver_2_ug_fb <- intersect(pat_2_liver_2_ug$location, pat_2_liver_2_freebayes$location) #6502
pat_2_liver_2_ug_vs <- intersect(pat_2_liver_2_ug$location, pat_2_liver_2_varscan$location) #1797
pat_2_liver_2_vs_fb <- intersect(pat_2_liver_2_varscan$location, pat_2_liver_2_freebayes$location) #1710

#intersect b/w mutect calls and two out of the 3 others
pat_2_liver_2_mutect_ug_fb <- intersect(pat_2_liver_2_mutect$location, pat_2_liver_2_ug_fb) #2515
pat_2_liver_2_mutect_ug_vs <- intersect(pat_2_liver_2_mutect$location, pat_2_liver_2_ug_vs) #345
pat_2_liver_2_mutect_vs_fb <- intersect(pat_2_liver_2_mutect$location, pat_2_liver_2_vs_fb) #314

#combine them
pat_2_liver_2_all <- c(pat_2_liver_2_mutect_ug_fb, pat_2_liver_2_mutect_ug_vs, pat_2_liver_2_mutect_vs_fb)

#subset mutect calls
pat_2_liver_2_mutect_short <- pat_2_liver_2_mutect[pat_2_liver_2_mutect$location %in% pat_2_liver_2_all, ] #2584

#trim to trusight gene calls
pat_2_liver_2_mutect_trusight <- trusight_from_mutect(pat_2_liver_2_mutect_short) #17

## plasma----
# import data
pat_2_plasma_ug <- read.delim('Data/pat_2_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #5336
pat_2_plasma_varscan <- read.delim('Data/pat_2_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #5253
pat_2_plasma_freebayes <- read.delim('Data/pat_2_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #12304

pat_2_plasma_mutect <- read.delim('Data/pat_2_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #20973

# dummy variable for matching
pat_2_plasma_ug$location <- paste0(pat_2_plasma_ug$CHROM, ':', pat_2_plasma_ug$POS)
pat_2_plasma_varscan$location <- paste0(pat_2_plasma_varscan$CHROM, ':', pat_2_plasma_varscan$POS)
pat_2_plasma_freebayes$location <- paste0(pat_2_plasma_freebayes$CHROM, ':', pat_2_plasma_freebayes$POS)

pat_2_plasma_mutect$location <- paste0(pat_2_plasma_mutect$CHROM, ':', pat_2_plasma_mutect$POS)

#pairwise intersections for the 3
pat_2_plasma_ug_fb <- intersect(pat_2_plasma_ug$location, pat_2_plasma_freebayes$location) #4133
pat_2_plasma_ug_vs <- intersect(pat_2_plasma_ug$location, pat_2_plasma_varscan$location) #901
pat_2_plasma_vs_fb <- intersect(pat_2_plasma_varscan$location, pat_2_plasma_freebayes$location) #942

#intersect b/w mutect calls and two out of the 3 others
pat_2_plasma_mutect_ug_fb <- intersect(pat_2_plasma_mutect$location, pat_2_plasma_ug_fb) #1944
pat_2_plasma_mutect_ug_vs <- intersect(pat_2_plasma_mutect$location, pat_2_plasma_ug_vs) #563
pat_2_plasma_mutect_vs_fb <- intersect(pat_2_plasma_mutect$location, pat_2_plasma_vs_fb) #535

#combine them
pat_2_plasma_all <- c(pat_2_plasma_mutect_ug_fb, pat_2_plasma_mutect_ug_vs, pat_2_plasma_mutect_vs_fb)

#subset mutect calls
pat_2_plasma_mutect_short <- pat_2_plasma_mutect[pat_2_plasma_mutect$location %in% pat_2_plasma_all, ] #1996

#trim to trusight gene calls
pat_2_plasma_mutect_trusight <- trusight_from_mutect(pat_2_plasma_mutect_short) #45

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 mets
pat_2_met_pool <- c(pat_2_breast_1_mutect_trusight$location, pat_2_breast_2_mutect_trusight$location, pat_2_liver_1_mutect_trusight$location, pat_2_liver_2_mutect_trusight$location)
pat_2_met_pool <- unique(pat_2_met_pool) #35

pat_2_plasma_link <- pat_2_plasma_mutect_trusight[pat_2_plasma_mutect_trusight$location %in% pat_2_met_pool, ] #9 mutations in 7 genes




## PATIENT EMA ----
## heart----
# import data
pat_ema_heart_ug <- read.delim('Data/pat_ema_heart_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #8352
pat_ema_heart_varscan <- read.delim('Data/pat_ema_heart_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3585
pat_ema_heart_freebayes <- read.delim('Data/pat_ema_heart_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #69805

pat_ema_heart_mutect <- read.delim('Data/pat_ema_heart_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #100948

# dummy variable for matching
pat_ema_heart_ug$location <- paste0(pat_ema_heart_ug$CHROM, ':', pat_ema_heart_ug$POS)
pat_ema_heart_varscan$location <- paste0(pat_ema_heart_varscan$CHROM, ':', pat_ema_heart_varscan$POS)
pat_ema_heart_freebayes$location <- paste0(pat_ema_heart_freebayes$CHROM, ':', pat_ema_heart_freebayes$POS)

pat_ema_heart_mutect$location <- paste0(pat_ema_heart_mutect$CHROM, ':', pat_ema_heart_mutect$POS)

#pairwise intersections for the 3
pat_ema_heart_ug_fb <- intersect(pat_ema_heart_ug$location, pat_ema_heart_freebayes$location) #7371
pat_ema_heart_ug_vs <- intersect(pat_ema_heart_ug$location, pat_ema_heart_varscan$location) #2070
pat_ema_heart_vs_fb <- intersect(pat_ema_heart_varscan$location, pat_ema_heart_freebayes$location) #1997

#intersect b/w mutect calls and two out of the 3 others
pat_ema_heart_mutect_ug_fb <- intersect(pat_ema_heart_mutect$location, pat_ema_heart_ug_fb) #2964
pat_ema_heart_mutect_ug_vs <- intersect(pat_ema_heart_mutect$location, pat_ema_heart_ug_vs) #365
pat_ema_heart_mutect_vs_fb <- intersect(pat_ema_heart_mutect$location, pat_ema_heart_vs_fb) #324

#combine them
pat_ema_heart_all <- c(pat_ema_heart_mutect_ug_fb, pat_ema_heart_mutect_ug_vs, pat_ema_heart_mutect_vs_fb)

#subset mutect calls
pat_ema_heart_mutect_short <- pat_ema_heart_mutect[pat_ema_heart_mutect$location %in% pat_ema_heart_all, ] #2039

#trim to trusight gene calls
pat_ema_heart_mutect_trusight <- trusight_from_mutect(pat_ema_heart_mutect_short) #15

## oment 1----
# import data
pat_ema_oment_1_ug <- read.delim('Data/pat_ema_oment_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #13204
pat_ema_oment_1_varscan <- read.delim('Data/pat_ema_oment_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4602
pat_ema_oment_1_freebayes <- read.delim('Data/pat_ema_oment_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #143272

pat_ema_oment_1_mutect <- read.delim('Data/pat_ema_oment_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #181967

# dummy variable for matching
pat_ema_oment_1_ug$location <- paste0(pat_ema_oment_1_ug$CHROM, ':', pat_ema_oment_1_ug$POS)
pat_ema_oment_1_varscan$location <- paste0(pat_ema_oment_1_varscan$CHROM, ':', pat_ema_oment_1_varscan$POS)
pat_ema_oment_1_freebayes$location <- paste0(pat_ema_oment_1_freebayes$CHROM, ':', pat_ema_oment_1_freebayes$POS)

pat_ema_oment_1_mutect$location <- paste0(pat_ema_oment_1_mutect$CHROM, ':', pat_ema_oment_1_mutect$POS)

#pairwise intersections for the 3
pat_ema_oment_1_ug_fb <- intersect(pat_ema_oment_1_ug$location, pat_ema_oment_1_freebayes$location) #11915
pat_ema_oment_1_ug_vs <- intersect(pat_ema_oment_1_ug$location, pat_ema_oment_1_varscan$location) #2318
pat_ema_oment_1_vs_fb <- intersect(pat_ema_oment_1_varscan$location, pat_ema_oment_1_freebayes$location) #2191

#intersect b/w mutect calls and two out of the 3 others
pat_ema_oment_1_mutect_ug_fb <- intersect(pat_ema_oment_1_mutect$location, pat_ema_oment_1_ug_fb) #6734
pat_ema_oment_1_mutect_ug_vs <- intersect(pat_ema_oment_1_mutect$location, pat_ema_oment_1_ug_vs) #539
pat_ema_oment_1_mutect_vs_fb <- intersect(pat_ema_oment_1_mutect$location, pat_ema_oment_1_vs_fb) #484

#combine them
pat_ema_oment_1_all <- c(pat_ema_oment_1_mutect_ug_fb, pat_ema_oment_1_mutect_ug_vs, pat_ema_oment_1_mutect_vs_fb)

#subset mutect calls
pat_ema_oment_1_mutect_short <- pat_ema_oment_1_mutect[pat_ema_oment_1_mutect$location %in% pat_ema_oment_1_all, ] #6863

#trim to trusight gene calls
pat_ema_oment_1_mutect_trusight <- trusight_from_mutect(pat_ema_oment_1_mutect_short) #13

## oment 2----
# import data
pat_ema_oment_2_ug <- read.delim('Data/pat_ema_oment_2_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11632
pat_ema_oment_2_varscan <- read.delim('Data/pat_ema_oment_2_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #5282
pat_ema_oment_2_freebayes <- read.delim('Data/pat_ema_oment_2_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #106115

pat_ema_oment_2_mutect <- read.delim('Data/pat_ema_oment_2_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #174243

# dummy variable for matching
pat_ema_oment_2_ug$location <- paste0(pat_ema_oment_2_ug$CHROM, ':', pat_ema_oment_2_ug$POS)
pat_ema_oment_2_varscan$location <- paste0(pat_ema_oment_2_varscan$CHROM, ':', pat_ema_oment_2_varscan$POS)
pat_ema_oment_2_freebayes$location <- paste0(pat_ema_oment_2_freebayes$CHROM, ':', pat_ema_oment_2_freebayes$POS)

pat_ema_oment_2_mutect$location <- paste0(pat_ema_oment_2_mutect$CHROM, ':', pat_ema_oment_2_mutect$POS)

#pairwise intersections for the 3
pat_ema_oment_2_ug_fb <- intersect(pat_ema_oment_2_ug$location, pat_ema_oment_2_freebayes$location) #10483
pat_ema_oment_2_ug_vs <- intersect(pat_ema_oment_2_ug$location, pat_ema_oment_2_varscan$location) #2554
pat_ema_oment_2_vs_fb <- intersect(pat_ema_oment_2_varscan$location, pat_ema_oment_2_freebayes$location) #2536

#intersect b/w mutect calls and two out of the 3 others
pat_ema_oment_2_mutect_ug_fb <- intersect(pat_ema_oment_2_mutect$location, pat_ema_oment_2_ug_fb) #5333
pat_ema_oment_2_mutect_ug_vs <- intersect(pat_ema_oment_2_mutect$location, pat_ema_oment_2_ug_vs) #568
pat_ema_oment_2_mutect_vs_fb <- intersect(pat_ema_oment_2_mutect$location, pat_ema_oment_2_vs_fb) #500

#combine them
pat_ema_oment_2_all <- c(pat_ema_oment_2_mutect_ug_fb, pat_ema_oment_2_mutect_ug_vs, pat_ema_oment_2_mutect_vs_fb)

#subset mutect calls
pat_ema_oment_2_mutect_short <- pat_ema_oment_2_mutect[pat_ema_oment_2_mutect$location %in% pat_ema_oment_2_all, ] #5473

#trim to trusight gene calls
pat_ema_oment_2_mutect_trusight <- trusight_from_mutect(pat_ema_oment_2_mutect_short) #16

## liver 1----
# import data
pat_ema_liver_1_ug <- read.delim('Data/pat_ema_liver_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11065
pat_ema_liver_1_varscan <- read.delim('Data/pat_ema_liver_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4095
pat_ema_liver_1_freebayes <- read.delim('Data/pat_ema_liver_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #147774

pat_ema_liver_1_mutect <- read.delim('Data/pat_ema_liver_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #130183

# dummy variable for matching
pat_ema_liver_1_ug$location <- paste0(pat_ema_liver_1_ug$CHROM, ':', pat_ema_liver_1_ug$POS)
pat_ema_liver_1_varscan$location <- paste0(pat_ema_liver_1_varscan$CHROM, ':', pat_ema_liver_1_varscan$POS)
pat_ema_liver_1_freebayes$location <- paste0(pat_ema_liver_1_freebayes$CHROM, ':', pat_ema_liver_1_freebayes$POS)

pat_ema_liver_1_mutect$location <- paste0(pat_ema_liver_1_mutect$CHROM, ':', pat_ema_liver_1_mutect$POS)

#pairwise intersections for the 3
pat_ema_liver_1_ug_fb <- intersect(pat_ema_liver_1_ug$location, pat_ema_liver_1_freebayes$location) #9793
pat_ema_liver_1_ug_vs <- intersect(pat_ema_liver_1_ug$location, pat_ema_liver_1_varscan$location) #1549
pat_ema_liver_1_vs_fb <- intersect(pat_ema_liver_1_varscan$location, pat_ema_liver_1_freebayes$location) #1489

#intersect b/w mutect calls and two out of the 3 others
pat_ema_liver_1_mutect_ug_fb <- intersect(pat_ema_liver_1_mutect$location, pat_ema_liver_1_ug_fb) #5489
pat_ema_liver_1_mutect_ug_vs <- intersect(pat_ema_liver_1_mutect$location, pat_ema_liver_1_ug_vs) #382
pat_ema_liver_1_mutect_vs_fb <- intersect(pat_ema_liver_1_mutect$location, pat_ema_liver_1_vs_fb) #312

#combine them
pat_ema_liver_1_all <- c(pat_ema_liver_1_mutect_ug_fb, pat_ema_liver_1_mutect_ug_vs, pat_ema_liver_1_mutect_vs_fb)

#subset mutect calls
pat_ema_liver_1_mutect_short <- pat_ema_liver_1_mutect[pat_ema_liver_1_mutect$location %in% pat_ema_liver_1_all, ] #5585

#trim to trusight gene calls
pat_ema_liver_1_mutect_trusight <- trusight_from_mutect(pat_ema_liver_1_mutect_short) #5

## liver 2----
# import data
pat_ema_liver_2_ug <- read.delim('Data/pat_ema_liver_2_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #10404
pat_ema_liver_2_varscan <- read.delim('Data/pat_ema_liver_2_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #5252
pat_ema_liver_2_freebayes <- read.delim('Data/pat_ema_liver_2_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #114670

pat_ema_liver_2_mutect <- read.delim('Data/pat_ema_liver_2_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #137556

# dummy variable for matching
pat_ema_liver_2_ug$location <- paste0(pat_ema_liver_2_ug$CHROM, ':', pat_ema_liver_2_ug$POS)
pat_ema_liver_2_varscan$location <- paste0(pat_ema_liver_2_varscan$CHROM, ':', pat_ema_liver_2_varscan$POS)
pat_ema_liver_2_freebayes$location <- paste0(pat_ema_liver_2_freebayes$CHROM, ':', pat_ema_liver_2_freebayes$POS)

pat_ema_liver_2_mutect$location <- paste0(pat_ema_liver_2_mutect$CHROM, ':', pat_ema_liver_2_mutect$POS)

#pairwise intersections for the 3
pat_ema_liver_2_ug_fb <- intersect(pat_ema_liver_2_ug$location, pat_ema_liver_2_freebayes$location) #9329
pat_ema_liver_2_ug_vs <- intersect(pat_ema_liver_2_ug$location, pat_ema_liver_2_varscan$location) #2300
pat_ema_liver_2_vs_fb <- intersect(pat_ema_liver_2_varscan$location, pat_ema_liver_2_freebayes$location) #2247

#intersect b/w mutect calls and two out of the 3 others
pat_ema_liver_2_mutect_ug_fb <- intersect(pat_ema_liver_2_mutect$location, pat_ema_liver_2_ug_fb) #4602
pat_ema_liver_2_mutect_ug_vs <- intersect(pat_ema_liver_2_mutect$location, pat_ema_liver_2_ug_vs) #557
pat_ema_liver_2_mutect_vs_fb <- intersect(pat_ema_liver_2_mutect$location, pat_ema_liver_2_vs_fb) #312

#combine them
pat_ema_liver_2_all <- c(pat_ema_liver_2_mutect_ug_fb, pat_ema_liver_2_mutect_ug_vs, pat_ema_liver_2_mutect_vs_fb)

#subset mutect calls
pat_ema_liver_2_mutect_short <- pat_ema_liver_2_mutect[pat_ema_liver_2_mutect$location %in% pat_ema_liver_2_all, ] #4729

#trim to trusight gene calls
pat_ema_liver_2_mutect_trusight <- trusight_from_mutect(pat_ema_liver_2_mutect_short) #13

## left kidney----
# import data
pat_ema_l_kidney_ug <- read.delim('Data/pat_ema_l_kidney_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11508
pat_ema_l_kidney_varscan <- read.delim('Data/pat_ema_l_kidney_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4411
pat_ema_l_kidney_freebayes <- read.delim('Data/pat_ema_l_kidney_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #122035

pat_ema_l_kidney_mutect <- read.delim('Data/pat_ema_l_kidney_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #151398

# dummy variable for matching
pat_ema_l_kidney_ug$location <- paste0(pat_ema_l_kidney_ug$CHROM, ':', pat_ema_l_kidney_ug$POS)
pat_ema_l_kidney_varscan$location <- paste0(pat_ema_l_kidney_varscan$CHROM, ':', pat_ema_l_kidney_varscan$POS)
pat_ema_l_kidney_freebayes$location <- paste0(pat_ema_l_kidney_freebayes$CHROM, ':', pat_ema_l_kidney_freebayes$POS)

pat_ema_l_kidney_mutect$location <- paste0(pat_ema_l_kidney_mutect$CHROM, ':', pat_ema_l_kidney_mutect$POS)

#pairwise intersections for the 3
pat_ema_l_kidney_ug_fb <- intersect(pat_ema_l_kidney_ug$location, pat_ema_l_kidney_freebayes$location) #10483
pat_ema_l_kidney_ug_vs <- intersect(pat_ema_l_kidney_ug$location, pat_ema_l_kidney_varscan$location) #2230
pat_ema_l_kidney_vs_fb <- intersect(pat_ema_l_kidney_varscan$location, pat_ema_l_kidney_freebayes$location) #2128

#intersect b/w mutect calls and two out of the 3 others
pat_ema_l_kidney_mutect_ug_fb <- intersect(pat_ema_l_kidney_mutect$location, pat_ema_l_kidney_ug_fb) #5365
pat_ema_l_kidney_mutect_ug_vs <- intersect(pat_ema_l_kidney_mutect$location, pat_ema_l_kidney_ug_vs) #512
pat_ema_l_kidney_mutect_vs_fb <- intersect(pat_ema_l_kidney_mutect$location, pat_ema_l_kidney_vs_fb) #436

#combine them
pat_ema_l_kidney_all <- c(pat_ema_l_kidney_mutect_ug_fb, pat_ema_l_kidney_mutect_ug_vs, pat_ema_l_kidney_mutect_vs_fb)

#subset mutect calls
pat_ema_l_kidney_mutect_short <- pat_ema_l_kidney_mutect[pat_ema_l_kidney_mutect$location %in% pat_ema_l_kidney_all, ] #5495

#trim to trusight gene calls
pat_ema_l_kidney_mutect_trusight <- trusight_from_mutect(pat_ema_l_kidney_mutect_short) #8

## right kidney----
# import data
pat_ema_r_kidney_ug <- read.delim('Data/pat_ema_r_kidney_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #6608
pat_ema_r_kidney_varscan <- read.delim('Data/pat_ema_r_kidney_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #2703
pat_ema_r_kidney_freebayes <- read.delim('Data/pat_ema_r_kidney_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #69009

pat_ema_r_kidney_mutect <- read.delim('Data/pat_ema_r_kidney_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #65907

# dummy variable for matching
pat_ema_r_kidney_ug$location <- paste0(pat_ema_r_kidney_ug$CHROM, ':', pat_ema_r_kidney_ug$POS)
pat_ema_r_kidney_varscan$location <- paste0(pat_ema_r_kidney_varscan$CHROM, ':', pat_ema_r_kidney_varscan$POS)
pat_ema_r_kidney_freebayes$location <- paste0(pat_ema_r_kidney_freebayes$CHROM, ':', pat_ema_r_kidney_freebayes$POS)

pat_ema_r_kidney_mutect$location <- paste0(pat_ema_r_kidney_mutect$CHROM, ':', pat_ema_r_kidney_mutect$POS)

#pairwise intersections for the 3
pat_ema_r_kidney_ug_fb <- intersect(pat_ema_r_kidney_ug$location, pat_ema_r_kidney_freebayes$location) #5874
pat_ema_r_kidney_ug_vs <- intersect(pat_ema_r_kidney_ug$location, pat_ema_r_kidney_varscan$location) #1569
pat_ema_r_kidney_vs_fb <- intersect(pat_ema_r_kidney_varscan$location, pat_ema_r_kidney_freebayes$location) #1508

#intersect b/w mutect calls and two out of the 3 others
pat_ema_r_kidney_mutect_ug_fb <- intersect(pat_ema_r_kidney_mutect$location, pat_ema_r_kidney_ug_fb) #2163
pat_ema_r_kidney_mutect_ug_vs <- intersect(pat_ema_r_kidney_mutect$location, pat_ema_r_kidney_ug_vs) #259
pat_ema_r_kidney_mutect_vs_fb <- intersect(pat_ema_r_kidney_mutect$location, pat_ema_r_kidney_vs_fb) #226

#combine them
pat_ema_r_kidney_all <- c(pat_ema_r_kidney_mutect_ug_fb, pat_ema_r_kidney_mutect_ug_vs, pat_ema_r_kidney_mutect_vs_fb)

#subset mutect calls
pat_ema_r_kidney_mutect_short <- pat_ema_r_kidney_mutect[pat_ema_r_kidney_mutect$location %in% pat_ema_r_kidney_all, ] #2210

#trim to trusight gene calls
pat_ema_r_kidney_mutect_trusight <- trusight_from_mutect(pat_ema_r_kidney_mutect_short) #1

## plasma----
# import data
pat_ema_plasma_ug <- read.delim('Data/pat_ema_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #52949
pat_ema_plasma_varscan <- read.delim('Data/pat_ema_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #25333
pat_ema_plasma_freebayes <- read.delim('Data/pat_ema_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #391305

pat_ema_plasma_mutect <- read.delim('Data/pat_ema_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #693029

# dummy variable for matching
pat_ema_plasma_ug$location <- paste0(pat_ema_plasma_ug$CHROM, ':', pat_ema_plasma_ug$POS)
pat_ema_plasma_varscan$location <- paste0(pat_ema_plasma_varscan$CHROM, ':', pat_ema_plasma_varscan$POS)
pat_ema_plasma_freebayes$location <- paste0(pat_ema_plasma_freebayes$CHROM, ':', pat_ema_plasma_freebayes$POS)

pat_ema_plasma_mutect$location <- paste0(pat_ema_plasma_mutect$CHROM, ':', pat_ema_plasma_mutect$POS)

#pairwise intersections for the 3
pat_ema_plasma_ug_fb <- intersect(pat_ema_plasma_ug$location, pat_ema_plasma_freebayes$location) #44276
pat_ema_plasma_ug_vs <- intersect(pat_ema_plasma_ug$location, pat_ema_plasma_varscan$location) #6107
pat_ema_plasma_vs_fb <- intersect(pat_ema_plasma_varscan$location, pat_ema_plasma_freebayes$location) #5950

#intersect b/w mutect calls and two out of the 3 others
pat_ema_plasma_mutect_ug_fb <- intersect(pat_ema_plasma_mutect$location, pat_ema_plasma_ug_fb) #32579
pat_ema_plasma_mutect_ug_vs <- intersect(pat_ema_plasma_mutect$location, pat_ema_plasma_ug_vs) #2997
pat_ema_plasma_mutect_vs_fb <- intersect(pat_ema_plasma_mutect$location, pat_ema_plasma_vs_fb) #2763

#combine them
pat_ema_plasma_all <- c(pat_ema_plasma_mutect_ug_fb, pat_ema_plasma_mutect_ug_vs, pat_ema_plasma_mutect_vs_fb)

#subset mutect calls
pat_ema_plasma_mutect_short <- pat_ema_plasma_mutect[pat_ema_plasma_mutect$location %in% pat_ema_plasma_all, ] #33023

#trim to trusight gene calls
pat_ema_plasma_mutect_trusight <- trusight_from_mutect(pat_ema_plasma_mutect_short) #78

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 mets
pat_ema_met_pool <- c(pat_ema_heart_mutect_trusight$location, pat_ema_oment_1_mutect_trusight$location, pat_ema_oment_2_mutect_trusight$location, 
                      pat_ema_liver_1_mutect_trusight$location, pat_ema_liver_2_mutect_trusight$location, 
                      pat_ema_l_kidney_mutect_trusight$location, pat_ema_r_kidney_mutect_trusight$location)
pat_ema_met_pool <- unique(pat_ema_met_pool) #31

pat_ema_plasma_link <- pat_ema_plasma_mutect_trusight[pat_ema_plasma_mutect_trusight$location %in% pat_ema_met_pool, ] #18 mutations in 7 genes





## JUST USING MUTECT WITH ONE INTERSECTION -----
## trusight only ----
## patient 10 ----
## liver 1 ----
#how many trusight calls in each method?
pat_10_liver_1_fb_ts <- trusight_from_mutect(pat_10_liver_1_freebayes) #126
pat_10_liver_1_mu_ts <- trusight_from_mutect(pat_10_liver_1_mutect) #107
pat_10_liver_1_ug_ts <- trusight_from_mutect(pat_10_liver_1_ug) #86
pat_10_liver_1_vs_ts <- trusight_from_mutect(pat_10_liver_1_varscan) #103

#how many in common?
pat_10_liver_1_single_intersect <- Reduce(intersect, list(pat_10_liver_1_fb_ts$location, pat_10_liver_1_mu_ts$location, pat_10_liver_1_ug_ts$location, pat_10_liver_1_vs_ts$location)) #4

## intersect with mutect ----
#how many in freebayes and mutect?
pat_10_liver_1_mu_fb <- intersect(pat_10_liver_1_freebayes$location, pat_10_liver_1_mutect$location) #78343 
pat_10_liver_1_mutect_fb <- pat_10_liver_1_mutect[pat_10_liver_1_mutect$location %in% pat_10_liver_1_mu_fb, ]

#how many trusight calls?
pat_10_liver_1_mutect_fb_ts <- trusight_from_mutect(pat_10_liver_1_mutect_fb) #42

#how many in ug and mutect?
pat_10_liver_1_mu_ug <- intersect(pat_10_liver_1_ug$location, pat_10_liver_1_mutect$location) #6946 
pat_10_liver_1_mutect_ug <- pat_10_liver_1_mutect[pat_10_liver_1_mutect$location %in% pat_10_liver_1_mu_ug, ]

#how many trusight calls?
pat_10_liver_1_mutect_ug_ts <- trusight_from_mutect(pat_10_liver_1_mutect_ug) #19

#how many in vs and mutect?
pat_10_liver_1_mu_vs <- intersect(pat_10_liver_1_varscan$location, pat_10_liver_1_mutect$location) #961 
pat_10_liver_1_mutect_vs <- pat_10_liver_1_mutect[pat_10_liver_1_mutect$location %in% pat_10_liver_1_mu_vs, ]

#how many trusight calls?
pat_10_liver_1_mutect_vs_ts <- trusight_from_mutect(pat_10_liver_1_mutect_vs) #7


#how many mutect hit at least one other method?
pat_10_liver_1_all_calls_mutect <- unique(c(pat_10_liver_1_mutect_fb$location, pat_10_liver_1_mutect_ug$location, pat_10_liver_1_mutect_vs$location)) #79017
pat_10_liver_1_mutect_all <- pat_10_liver_1_mutect[pat_10_liver_1_mutect$location %in% pat_10_liver_1_all_calls_mutect, ]

#how many trusight calls?
pat_10_liver_1_mutect_all_ts <- trusight_from_mutect(pat_10_liver_1_mutect_all) #47

## liver 2a ----
#how many trusight calls in each method?
pat_10_liver_2a_fb_ts <- trusight_from_mutect(pat_10_liver_2a_freebayes) #106
pat_10_liver_2a_mu_ts <- trusight_from_mutect(pat_10_liver_2a_mutect) #71
pat_10_liver_2a_ug_ts <- trusight_from_mutect(pat_10_liver_2a_ug) #79
pat_10_liver_2a_vs_ts <- trusight_from_mutect(pat_10_liver_2a_varscan) #187

#how many in common?
pat_10_liver_2a_single_intersect <- Reduce(intersect, list(pat_10_liver_2a_fb_ts$location, pat_10_liver_2a_mu_ts$location, pat_10_liver_2a_ug_ts$location, pat_10_liver_2a_vs_ts$location)) #2

## intersect with mutect ----
#how many in freebayes and mutect?
pat_10_liver_2a_mu_fb <- intersect(pat_10_liver_2a_freebayes$location, pat_10_liver_2a_mutect$location) #59422 
pat_10_liver_2a_mutect_fb <- pat_10_liver_2a_mutect[pat_10_liver_2a_mutect$location %in% pat_10_liver_2a_mu_fb, ]

#how many trusight calls?
pat_10_liver_2a_mutect_fb_ts <- trusight_from_mutect(pat_10_liver_2a_mutect_fb) #31

#how many in ug and mutect?
pat_10_liver_2a_mu_ug <- intersect(pat_10_liver_2a_ug$location, pat_10_liver_2a_mutect$location) #3951 
pat_10_liver_2a_mutect_ug <- pat_10_liver_2a_mutect[pat_10_liver_2a_mutect$location %in% pat_10_liver_2a_mu_ug, ]

#how many trusight calls?
pat_10_liver_2a_mutect_ug_ts <- trusight_from_mutect(pat_10_liver_2a_mutect_ug) #13

#how many in vs and mutect?
pat_10_liver_2a_mu_vs <- intersect(pat_10_liver_2a_varscan$location, pat_10_liver_2a_mutect$location) #551 
pat_10_liver_2a_mutect_vs <- pat_10_liver_2a_mutect[pat_10_liver_2a_mutect$location %in% pat_10_liver_2a_mu_vs, ]

#how many trusight calls?
pat_10_liver_2a_mutect_vs_ts <- trusight_from_mutect(pat_10_liver_2a_mutect_vs) #4


#how many mutect hit at least one other method?
pat_10_liver_2a_all_calls_mutect <- unique(c(pat_10_liver_2a_mutect_fb$location, pat_10_liver_2a_mutect_ug$location, pat_10_liver_2a_mutect_vs$location)) #59867
pat_10_liver_2a_mutect_all <- pat_10_liver_2a_mutect[pat_10_liver_2a_mutect$location %in% pat_10_liver_2a_all_calls_mutect, ]

#how many trusight calls?
pat_10_liver_2a_mutect_all_ts <- trusight_from_mutect(pat_10_liver_2a_mutect_all) #36

## liver 5 ----
#how many trusight calls in each method?
pat_10_liver_5_fb_ts <- trusight_from_mutect(pat_10_liver_5_freebayes) #120
pat_10_liver_5_mu_ts <- trusight_from_mutect(pat_10_liver_5_mutect) #91
pat_10_liver_5_ug_ts <- trusight_from_mutect(pat_10_liver_5_ug) #84
pat_10_liver_5_vs_ts <- trusight_from_mutect(pat_10_liver_5_varscan) #85

#how many in common?
pat_10_liver_5_single_intersect <- Reduce(intersect, list(pat_10_liver_5_fb_ts$location, pat_10_liver_5_mu_ts$location, pat_10_liver_5_ug_ts$location, pat_10_liver_5_vs_ts$location)) #5

## intersect with mutect ----
#how many in freebayes and mutect?
pat_10_liver_5_mu_fb <- intersect(pat_10_liver_5_freebayes$location, pat_10_liver_5_mutect$location) #85631 
pat_10_liver_5_mutect_fb <- pat_10_liver_5_mutect[pat_10_liver_5_mutect$location %in% pat_10_liver_5_mu_fb, ]

#how many trusight calls?
pat_10_liver_5_mutect_fb_ts <- trusight_from_mutect(pat_10_liver_5_mutect_fb) #46

#how many in ug and mutect?
pat_10_liver_5_mu_ug <- intersect(pat_10_liver_5_ug$location, pat_10_liver_5_mutect$location) #7081 
pat_10_liver_5_mutect_ug <- pat_10_liver_5_mutect[pat_10_liver_5_mutect$location %in% pat_10_liver_5_mu_ug, ]

#how many trusight calls?
pat_10_liver_5_mutect_ug_ts <- trusight_from_mutect(pat_10_liver_5_mutect_ug) #20

#how many in vs and mutect?
pat_10_liver_5_mu_vs <- intersect(pat_10_liver_5_varscan$location, pat_10_liver_5_mutect$location) #792 
pat_10_liver_5_mutect_vs <- pat_10_liver_5_mutect[pat_10_liver_5_mutect$location %in% pat_10_liver_5_mu_vs, ]

#how many trusight calls?
pat_10_liver_5_mutect_vs_ts <- trusight_from_mutect(pat_10_liver_5_mutect_vs) #7


#how many mutect hit at least one other method?
pat_10_liver_5_all_calls_mutect <- unique(c(pat_10_liver_5_mutect_fb$location, pat_10_liver_5_mutect_ug$location, pat_10_liver_5_mutect_vs$location)) #86211
pat_10_liver_5_mutect_all <- pat_10_liver_5_mutect[pat_10_liver_5_mutect$location %in% pat_10_liver_5_all_calls_mutect, ]

#how many trusight calls?
pat_10_liver_5_mutect_all_ts <- trusight_from_mutect(pat_10_liver_5_mutect_all) #51

## plasma ----
#how many trusight calls in each method?
pat_10_plasma_fb_ts <- trusight_from_mutect(pat_10_plasma_freebayes) #197
pat_10_plasma_mu_ts <- trusight_from_mutect(pat_10_plasma_mutect) #1434
pat_10_plasma_ug_ts <- trusight_from_mutect(pat_10_plasma_ug) #383
pat_10_plasma_vs_ts <- trusight_from_mutect(pat_10_plasma_varscan) #1788

#how many in common?
pat_10_plasma_single_intersect <- Reduce(intersect, list(pat_10_plasma_fb_ts$location, pat_10_plasma_mu_ts$location, pat_10_plasma_ug_ts$location, pat_10_plasma_vs_ts$location)) #50

## intersect with mutect ----
#how many in freebayes and mutect?
pat_10_plasma_mu_fb <- intersect(pat_10_plasma_freebayes$location, pat_10_plasma_mutect$location) #259731 
pat_10_plasma_mutect_fb <- pat_10_plasma_mutect[pat_10_plasma_mutect$location %in% pat_10_plasma_mu_fb, ]

#how many trusight calls?
pat_10_plasma_mutect_fb_ts <- trusight_from_mutect(pat_10_plasma_mutect_fb) #152

#how many in ug and mutect?
pat_10_plasma_mu_ug <- intersect(pat_10_plasma_ug$location, pat_10_plasma_mutect$location) #27469
pat_10_plasma_mutect_ug <- pat_10_plasma_mutect[pat_10_plasma_mutect$location %in% pat_10_plasma_mu_ug, ]

#how many trusight calls?
pat_10_plasma_mutect_ug_ts <- trusight_from_mutect(pat_10_plasma_mutect_ug) #114

#how many in vs and mutect?
pat_10_plasma_mu_vs <- intersect(pat_10_plasma_varscan$location, pat_10_plasma_mutect$location) #13843
pat_10_plasma_mutect_vs <- pat_10_plasma_mutect[pat_10_plasma_mutect$location %in% pat_10_plasma_mu_vs, ]

#how many trusight calls?
pat_10_plasma_mutect_vs_ts <- trusight_from_mutect(pat_10_plasma_mutect_vs) #513


#how many mutect hit at least one other method?
pat_10_plasma_all_calls_mutect <- unique(c(pat_10_plasma_mutect_fb$location, pat_10_plasma_mutect_ug$location, pat_10_plasma_mutect_vs$location)) #272333
pat_10_plasma_mutect_all <- pat_10_plasma_mutect[pat_10_plasma_mutect$location %in% pat_10_plasma_all_calls_mutect, ]

#how many trusight calls?
pat_10_plasma_mutect_all_ts <- trusight_from_mutect(pat_10_plasma_mutect_all) #644


pat_10_met_pool_mutect <- c(pat_10_liver_1_mutect_all_ts$location, pat_10_liver_2a_mutect_all_ts$location, pat_10_liver_5_mutect_all_ts$location)
pat_10_met_pool_mutect <- unique(pat_10_met_pool_mutect) #103

pat_10_plasma_link_mutect <- pat_10_plasma_mutect_all_ts[pat_10_plasma_mutect_all_ts$location %in% pat_10_met_pool_mutect, ] #21 mutations in 9 genes


## trusight only ----
## patient 9 ----
## lymph node ----
#how many trusight calls in each method?
pat_9_ln_fb_ts <- trusight_from_mutect(pat_9_ln_freebayes) #85
pat_9_ln_mu_ts <- trusight_from_mutect(pat_9_ln_mutect) #63
pat_9_ln_ug_ts <- trusight_from_mutect(pat_9_ln_ug) #185
pat_9_ln_vs_ts <- trusight_from_mutect(pat_9_ln_varscan) #2

#how many in common?
pat_9_ln_single_intersect <- Reduce(intersect, list(pat_9_ln_fb_ts$location, pat_9_ln_mu_ts$location, pat_9_ln_ug_ts$location, pat_9_ln_vs_ts$location)) #NONE

## intersect with mutect ----
#how many in freebayes and mutect?
pat_9_ln_mu_fb <- intersect(pat_9_ln_freebayes$location, pat_9_ln_mutect$location) #1029 
pat_9_ln_mutect_fb <- pat_9_ln_mutect[pat_9_ln_mutect$location %in% pat_9_ln_mu_fb, ]

#how many trusight calls?
pat_9_ln_mutect_fb_ts <- trusight_from_mutect(pat_9_ln_mutect_fb) #15

#how many in ug and mutect?
pat_9_ln_mu_ug <- intersect(pat_9_ln_ug$location, pat_9_ln_mutect$location) #715
pat_9_ln_mutect_ug <- pat_9_ln_mutect[pat_9_ln_mutect$location %in% pat_9_ln_mu_ug, ]

#how many trusight calls?
pat_9_ln_mutect_ug_ts <- trusight_from_mutect(pat_9_ln_mutect_ug) #18

#how many in vs and mutect?
pat_9_ln_mu_vs <- intersect(pat_9_ln_varscan$location, pat_9_ln_mutect$location) #26 
pat_9_ln_mutect_vs <- pat_9_ln_mutect[pat_9_ln_mutect$location %in% pat_9_ln_mu_vs, ]

#how many trusight calls?
pat_9_ln_mutect_vs_ts <- trusight_from_mutect(pat_9_ln_mutect_vs) #0


#how many mutect hit at least one other method?
pat_9_ln_all_calls_mutect <- unique(c(pat_9_ln_mutect_fb$location, pat_9_ln_mutect_ug$location, pat_9_ln_mutect_vs$location)) #1273
pat_9_ln_mutect_all <- pat_9_ln_mutect[pat_9_ln_mutect$location %in% pat_9_ln_all_calls_mutect, ]

#how many trusight calls?
pat_9_ln_mutect_all_ts <- trusight_from_mutect(pat_9_ln_mutect_all) #23

## liver 2a ----
#how many trusight calls in each method?
pat_9_oment_fb_ts <- trusight_from_mutect(pat_9_oment_freebayes) #134
pat_9_oment_mu_ts <- trusight_from_mutect(pat_9_oment_mutect) #173
pat_9_oment_ug_ts <- trusight_from_mutect(pat_9_oment_ug) #119
pat_9_oment_vs_ts <- trusight_from_mutect(pat_9_oment_varscan) #179

#how many in common?
pat_9_oment_single_intersect <- Reduce(intersect, list(pat_9_oment_fb_ts$location, pat_9_oment_mu_ts$location, pat_9_oment_ug_ts$location, pat_9_oment_vs_ts$location)) #24

## intersect with mutect ----
#how many in freebayes and mutect?
pat_9_oment_mu_fb <- intersect(pat_9_oment_freebayes$location, pat_9_oment_mutect$location) #45393 
pat_9_oment_mutect_fb <- pat_9_oment_mutect[pat_9_oment_mutect$location %in% pat_9_oment_mu_fb, ]

#how many trusight calls?
pat_9_oment_mutect_fb_ts <- trusight_from_mutect(pat_9_oment_mutect_fb) #52

#how many in ug and mutect?
pat_9_oment_mu_ug <- intersect(pat_9_oment_ug$location, pat_9_oment_mutect$location) #4031 
pat_9_oment_mutect_ug <- pat_9_oment_mutect[pat_9_oment_mutect$location %in% pat_9_oment_mu_ug, ]

#how many trusight calls?
pat_9_oment_mutect_ug_ts <- trusight_from_mutect(pat_9_oment_mutect_ug) #44

#how many in vs and mutect?
pat_9_oment_mu_vs <- intersect(pat_9_oment_varscan$location, pat_9_oment_mutect$location) #771 
pat_9_oment_mutect_vs <- pat_9_oment_mutect[pat_9_oment_mutect$location %in% pat_9_oment_mu_vs, ]

#how many trusight calls?
pat_9_oment_mutect_vs_ts <- trusight_from_mutect(pat_9_oment_mutect_vs) #32


#how many mutect hit at least one other method?
pat_9_oment_all_calls_mutect <- unique(c(pat_9_oment_mutect_fb$location, pat_9_oment_mutect_ug$location, pat_9_oment_mutect_vs$location)) #46170
pat_9_oment_mutect_all <- pat_9_oment_mutect[pat_9_oment_mutect$location %in% pat_9_oment_all_calls_mutect, ]

#how many trusight calls?
pat_9_oment_mutect_all_ts <- trusight_from_mutect(pat_9_oment_mutect_all) #64

## ovary ----
#how many trusight calls in each method?
pat_9_ovary_fb_ts <- trusight_from_mutect(pat_9_ovary_freebayes) #101
pat_9_ovary_mu_ts <- trusight_from_mutect(pat_9_ovary_mutect) #102
pat_9_ovary_ug_ts <- trusight_from_mutect(pat_9_ovary_ug) #74
pat_9_ovary_vs_ts <- trusight_from_mutect(pat_9_ovary_varscan) #262

#how many in common?
pat_9_ovary_single_intersect <- Reduce(intersect, list(pat_9_ovary_fb_ts$location, pat_9_ovary_mu_ts$location, pat_9_ovary_ug_ts$location, pat_9_ovary_vs_ts$location)) #11

## intersect with mutect ----
#how many in freebayes and mutect?
pat_9_ovary_mu_fb <- intersect(pat_9_ovary_freebayes$location, pat_9_ovary_mutect$location) #57473 
pat_9_ovary_mutect_fb <- pat_9_ovary_mutect[pat_9_ovary_mutect$location %in% pat_9_ovary_mu_fb, ]

#how many trusight calls?
pat_9_ovary_mutect_fb_ts <- trusight_from_mutect(pat_9_ovary_mutect_fb) #31

#how many in ug and mutect?
pat_9_ovary_mu_ug <- intersect(pat_9_ovary_ug$location, pat_9_ovary_mutect$location) #3354 
pat_9_ovary_mutect_ug <- pat_9_ovary_mutect[pat_9_ovary_mutect$location %in% pat_9_ovary_mu_ug, ]

#how many trusight calls?
pat_9_ovary_mutect_ug_ts <- trusight_from_mutect(pat_9_ovary_mutect_ug) #21

#how many in vs and mutect?
pat_9_ovary_mu_vs <- intersect(pat_9_ovary_varscan$location, pat_9_ovary_mutect$location) #817 
pat_9_ovary_mutect_vs <- pat_9_ovary_mutect[pat_9_ovary_mutect$location %in% pat_9_ovary_mu_vs, ]

#how many trusight calls?
pat_9_ovary_mutect_vs_ts <- trusight_from_mutect(pat_9_ovary_mutect_vs) #29


#how many mutect hit at least one other method?
pat_9_ovary_all_calls_mutect <- unique(c(pat_9_ovary_mutect_fb$location, pat_9_ovary_mutect_ug$location, pat_9_ovary_mutect_vs$location)) #58134
pat_9_ovary_mutect_all <- pat_9_ovary_mutect[pat_9_ovary_mutect$location %in% pat_9_ovary_all_calls_mutect, ]

#how many trusight calls?
pat_9_ovary_mutect_all_ts <- trusight_from_mutect(pat_9_ovary_mutect_all) #50

## plasma ----
#how many trusight calls in each method?
pat_9_plasma_fb_ts <- trusight_from_mutect(pat_9_plasma_freebayes) #244
pat_9_plasma_mu_ts <- trusight_from_mutect(pat_9_plasma_mutect) #579
pat_9_plasma_ug_ts <- trusight_from_mutect(pat_9_plasma_ug) #173
pat_9_plasma_vs_ts <- trusight_from_mutect(pat_9_plasma_varscan) #390

#how many in common?
pat_9_plasma_single_intersect <- Reduce(intersect, list(pat_9_plasma_fb_ts$location, pat_9_plasma_mu_ts$location, pat_9_plasma_ug_ts$location, pat_9_plasma_vs_ts$location)) #35

## intersect with mutect ----
#how many in freebayes and mutect?
pat_9_plasma_mu_fb <- intersect(pat_9_plasma_freebayes$location, pat_9_plasma_mutect$location) #318281 
pat_9_plasma_mutect_fb <- pat_9_plasma_mutect[pat_9_plasma_mutect$location %in% pat_9_plasma_mu_fb, ]

#how many trusight calls?
pat_9_plasma_mutect_fb_ts <- trusight_from_mutect(pat_9_plasma_mutect_fb) #164

#how many in ug and mutect?
pat_9_plasma_mu_ug <- intersect(pat_9_plasma_ug$location, pat_9_plasma_mutect$location) #42290
pat_9_plasma_mutect_ug <- pat_9_plasma_mutect[pat_9_plasma_mutect$location %in% pat_9_plasma_mu_ug, ]

#how many trusight calls?
pat_9_plasma_mutect_ug_ts <- trusight_from_mutect(pat_9_plasma_mutect_ug) #65

#how many in vs and mutect?
pat_9_plasma_mu_vs <- intersect(pat_9_plasma_varscan$location, pat_9_plasma_mutect$location) #4267
pat_9_plasma_mutect_vs <- pat_9_plasma_mutect[pat_9_plasma_mutect$location %in% pat_9_plasma_mu_vs, ]

#how many trusight calls?
pat_9_plasma_mutect_vs_ts <- trusight_from_mutect(pat_9_plasma_mutect_vs) #60


#how many mutect hit at least one other method?
pat_9_plasma_all_calls_mutect <- unique(c(pat_9_plasma_mutect_fb$location, pat_9_plasma_mutect_ug$location, pat_9_plasma_mutect_vs$location)) #321199
pat_9_plasma_mutect_all <- pat_9_plasma_mutect[pat_9_plasma_mutect$location %in% pat_9_plasma_all_calls_mutect, ]

#how many trusight calls?
pat_9_plasma_mutect_all_ts <- trusight_from_mutect(pat_9_plasma_mutect_all) #195


pat_9_met_pool_mutect <- c(pat_9_ln_mutect_all_ts$location, pat_9_oment_mutect_all_ts$location, pat_9_ovary_mutect_all_ts$location)
pat_9_met_pool_mutect <- unique(pat_9_met_pool_mutect) #118

pat_9_plasma_link_mutect <- pat_9_plasma_mutect_all_ts[pat_9_plasma_mutect_all_ts$location %in% pat_9_met_pool_mutect, ] #40 mutations in 11 genes



## CHECKING HOW TO DRAW EVOL CLADE ----
## mutect trusight only
pat_ema_heart_mu_ts <- trusight_from_mutect(pat_ema_heart_mutect) #61
pat_ema_l_kidney_mu_ts <- trusight_from_mutect(pat_ema_l_kidney_mutect) #62
pat_ema_r_kidney_mu_ts <- trusight_from_mutect(pat_ema_r_kidney_mutect) #30
pat_ema_liver_1_mu_ts <- trusight_from_mutect(pat_ema_liver_1_mutect) #54
pat_ema_liver_2_mu_ts <- trusight_from_mutect(pat_ema_liver_2_mutect) #69
pat_ema_oment_1_mu_ts <- trusight_from_mutect(pat_ema_oment_1_mutect) #87
pat_ema_oment_2_mu_ts <- trusight_from_mutect(pat_ema_oment_2_mutect) #100





#how many in all? this is the length of the main branch
pat_ema_all_common <- Reduce(intersect, list(pat_9_ln_mu_ts$location, pat_9_oment_mu_ts$location, pat_9_ovary_mu_ts$location, pat_ema_heart_mu_ts$location, 
                                           pat_ema_l_kidney_mu_ts$location, pat_ema_r_kidney_mu_ts$location, pat_ema_liver_1_mu_ts$location, 
                                           pat_ema_liver_2_mu_ts$location, pat_ema_oment_1_mu_ts$location, pat_ema_oment_2_mu_ts$location)) #1
#which ones?
pat_9_mu_all <- pat_9_ln_mu_ts[pat_9_ln_mu_ts$location %in% pat_9_all_common, ]
#which in plasma?
length(intersect(pat_9_all_common, pat_9_plasma_mu_ts$location)) #2

#exclude those, these can be the lengths of the next branches if split solo
pat_9_ln_mu_not_all <- pat_9_ln_mu_ts[pat_9_ln_mu_ts$location %!in% pat_9_all_common, ] #59
pat_9_oment_mu_not_all <- pat_9_oment_mu_ts[pat_9_oment_mu_ts$location %!in% pat_9_all_common, ] #169
pat_9_ovary_mu_not_all <- pat_9_ovary_mu_ts[pat_9_ovary_mu_ts$location %!in% pat_9_all_common, ] #98

#gather all locations
pat_9_2nd_locations <- c(pat_9_ln_mu_not_all$location, pat_9_oment_mu_not_all$location, pat_9_ovary_mu_not_all$location)
pat_9_2nd_locations <- table(pat_9_2nd_locations)
pat_9_2nd_locations <- pat_9_2nd_locations[order(pat_9_2nd_locations, decreasing = TRUE)]
#only keep those that hit 2 out of 3
pat_9_2nd_locations_2 <- pat_9_2nd_locations[pat_9_2nd_locations == 2]
pat_9_2nd_locations_2 <- names(pat_9_2nd_locations_2)

#subset to those
pat_9_ln_2_of_3 <- pat_9_ln_mu_not_all[pat_9_ln_mu_not_all$location %in% pat_9_2nd_locations_2, ] #3
#ln excluding commons, length of ln solo branch
pat_9_ln_only <- pat_9_ln_mu_not_all[pat_9_ln_mu_not_all$location %!in% pat_9_2nd_locations_2, ] #56
# in plasma?
length(intersect(pat_9_ln_only$location, pat_9_plasma_mu_ts$location))

pat_9_oment_2_of_3 <- pat_9_oment_mu_not_all[pat_9_oment_mu_not_all$location %in% pat_9_2nd_locations_2, ] #13
pat_9_ovary_2_of_3 <- pat_9_ovary_mu_not_all[pat_9_ovary_mu_not_all$location %in% pat_9_2nd_locations_2, ] #12
#ln is most different

#other 2 but not ln
pat_9_oment_not_ln <- pat_9_oment_2_of_3[pat_9_oment_2_of_3$location %!in% pat_9_ln_2_of_3$location, ] #11
pat_9_ovary_not_ln <- pat_9_ovary_2_of_3[pat_9_ovary_2_of_3$location %!in% pat_9_ln_2_of_3$location, ] #11
#in common?
length(intersect(pat_9_oment_not_ln$location, pat_9_ovary_not_ln$location)) #all 11
#in plasma?
length(Reduce(intersect, list(pat_9_oment_not_ln$location, pat_9_ovary_not_ln$location, pat_9_plasma_mu_ts$location))) #all 11

#subset out ln 3
pat_9_oment_not_ln_long <- pat_9_oment_mu_not_all[pat_9_oment_mu_not_all$location %!in% pat_9_ln_mu_not_all$location, ] #167
pat_9_ovary_not_ln_long <- pat_9_ovary_mu_not_all[pat_9_ovary_mu_not_all$location %!in% pat_9_ln_mu_not_all$location, ] #97

#length of common branch
pat_9_ovary_oment <- intersect(pat_9_oment_not_ln_long$location, pat_9_ovary_not_ln_long$location) #11

#subset out those, lengths of solo branches
pat_9_oment_final_length <- pat_9_oment_not_ln_long[pat_9_oment_not_ln_long$location %!in% pat_9_ovary_oment, ] #156
pat_9_ovary_final_length <- pat_9_ovary_not_ln_long[pat_9_ovary_not_ln_long$location %!in% pat_9_ovary_oment, ] #86
#in plasma?
length(intersect(pat_9_oment_final_length$location, pat_9_plasma_mu_ts$location)) #19
length(intersect(pat_9_ovary_final_length$location, pat_9_plasma_mu_ts$location)) #9




## CHECKING HOW TO DRAW EVOL CLADE WITH EMA SAMPLES----
## mutect trusight only
#how many in all? this is the length of the main branch
pat_9_all_common <- Reduce(intersect, list(pat_9_ln_mu_ts$location, pat_9_oment_mu_ts$location, pat_9_ovary_mu_ts$location)) #4
#which ones?
pat_9_mu_all <- pat_9_ln_mu_ts[pat_9_ln_mu_ts$location %in% pat_9_all_common, ]
#which in plasma?
length(intersect(pat_9_all_common, pat_9_plasma_mu_ts$location)) #2

#exclude those, these can be the lengths of the next branches if split solo
pat_9_ln_mu_not_all <- pat_9_ln_mu_ts[pat_9_ln_mu_ts$location %!in% pat_9_all_common, ] #59
pat_9_oment_mu_not_all <- pat_9_oment_mu_ts[pat_9_oment_mu_ts$location %!in% pat_9_all_common, ] #169
pat_9_ovary_mu_not_all <- pat_9_ovary_mu_ts[pat_9_ovary_mu_ts$location %!in% pat_9_all_common, ] #98

#gather all locations
pat_9_2nd_locations <- c(pat_9_ln_mu_not_all$location, pat_9_oment_mu_not_all$location, pat_9_ovary_mu_not_all$location)
pat_9_2nd_locations <- table(pat_9_2nd_locations)
pat_9_2nd_locations <- pat_9_2nd_locations[order(pat_9_2nd_locations, decreasing = TRUE)]
#only keep those that hit 2 out of 3
pat_9_2nd_locations_2 <- pat_9_2nd_locations[pat_9_2nd_locations == 2]
pat_9_2nd_locations_2 <- names(pat_9_2nd_locations_2)

#subset to those
pat_9_ln_2_of_3 <- pat_9_ln_mu_not_all[pat_9_ln_mu_not_all$location %in% pat_9_2nd_locations_2, ] #3
#ln excluding commons, length of ln solo branch
pat_9_ln_only <- pat_9_ln_mu_not_all[pat_9_ln_mu_not_all$location %!in% pat_9_2nd_locations_2, ] #56
# in plasma?
length(intersect(pat_9_ln_only$location, pat_9_plasma_mu_ts$location))

pat_9_oment_2_of_3 <- pat_9_oment_mu_not_all[pat_9_oment_mu_not_all$location %in% pat_9_2nd_locations_2, ] #13
pat_9_ovary_2_of_3 <- pat_9_ovary_mu_not_all[pat_9_ovary_mu_not_all$location %in% pat_9_2nd_locations_2, ] #12
#ln is most different

#other 2 but not ln
pat_9_oment_not_ln <- pat_9_oment_2_of_3[pat_9_oment_2_of_3$location %!in% pat_9_ln_2_of_3$location, ] #11
pat_9_ovary_not_ln <- pat_9_ovary_2_of_3[pat_9_ovary_2_of_3$location %!in% pat_9_ln_2_of_3$location, ] #11
#in common?
length(intersect(pat_9_oment_not_ln$location, pat_9_ovary_not_ln$location)) #all 11
#in plasma?
length(Reduce(intersect, list(pat_9_oment_not_ln$location, pat_9_ovary_not_ln$location, pat_9_plasma_mu_ts$location))) #all 11

#subset out ln 3
pat_9_oment_not_ln_long <- pat_9_oment_mu_not_all[pat_9_oment_mu_not_all$location %!in% pat_9_ln_mu_not_all$location, ] #167
pat_9_ovary_not_ln_long <- pat_9_ovary_mu_not_all[pat_9_ovary_mu_not_all$location %!in% pat_9_ln_mu_not_all$location, ] #97

#length of common branch
pat_9_ovary_oment <- intersect(pat_9_oment_not_ln_long$location, pat_9_ovary_not_ln_long$location) #11

#subset out those, lengths of solo branches
pat_9_oment_final_length <- pat_9_oment_not_ln_long[pat_9_oment_not_ln_long$location %!in% pat_9_ovary_oment, ] #156
pat_9_ovary_final_length <- pat_9_ovary_not_ln_long[pat_9_ovary_not_ln_long$location %!in% pat_9_ovary_oment, ] #86
#in plasma?
length(intersect(pat_9_oment_final_length$location, pat_9_plasma_mu_ts$location)) #19
length(intersect(pat_9_ovary_final_length$location, pat_9_plasma_mu_ts$location)) #9
