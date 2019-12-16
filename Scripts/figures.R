## this is the script for all figures that involve multiple patients
library(dplyr)
library(gplots)
library(ggplot2)
library(VennDiagram)
library(limma)
## add all scripts to load data

## big heatmap figure ----
# get names of genes in each sample
# pat 9
pat_9_ln_names <- names(table(pat_9_ln$gene_name))
pat_9_oment_names <- names(table(pat_9_oment$gene_name))
pat_9_ovary_names <- names(table(pat_9_ovary$gene_name))
pat_9_plasma_names <- names(table(pat_9_plasma$gene_name))

all_pat_9_names <- c(pat_9_ln_names, pat_9_oment_names, pat_9_ovary_names, pat_9_plasma_names)

#pat ema
pat_ema_heart_names <- names(table(pat_ema_heart$gene_name))
pat_ema_l_kidney_names <- names(table(pat_ema_l_kidney$gene_name))
pat_ema_r_kidney_names <- names(table(pat_ema_r_kidney$gene_name))
pat_ema_liver_1_names <- names(table(pat_ema_liver_1$gene_name))
pat_ema_liver_2_names <- names(table(pat_ema_liver_2$gene_name))
pat_ema_oment_1_names <- names(table(pat_ema_oment_1$gene_name))
pat_ema_oment_2_names <- names(table(pat_ema_oment_2$gene_name))
pat_ema_plasma_names <- names(table(pat_ema_plasma$gene_name))

all_pat_ema_names <- c(pat_ema_heart_names, pat_ema_l_kidney_names, pat_ema_r_kidney_names, 
                       pat_ema_liver_1_names, pat_ema_liver_2_names, pat_ema_oment_1_names, 
                       pat_ema_oment_2_names, pat_ema_plasma_names)

#pat 10
pat_10_liver_1_names <- names(table(pat_10_liver_1$gene_name))
pat_10_liver_2a_names <- names(table(pat_10_liver_2a$gene_name))
pat_10_liver_5_names <- names(table(pat_10_liver_5$gene_name))
pat_10_plasma_names <- names(table(pat_10_plasma$gene_name))

all_pat_10_names <- c(pat_10_liver_1_names, pat_10_liver_2a_names, pat_10_liver_5_names, pat_10_plasma_names)

#pat 8
pat_8_axillary_names <- names(table(pat_8_axillary$gene_name))
pat_8_breast_1_names <- names(table(pat_8_breast_1$gene_name))
pat_8_breast_2_names <- names(table(pat_8_breast_2$gene_name))
pat_8_plasma_names <- names(table(pat_8_plasma$gene_name))

all_pat_8_names <- c(pat_8_axillary_names, pat_8_breast_1_names, pat_8_breast_2_names, pat_8_plasma_names)

#pat 2
pat_2_breast_2_names <- names(table(pat_2_breast_2$gene_name))
pat_2_liver_1_names <- names(table(pat_2_liver_1$gene_name))
pat_2_liver_2_names <- names(pat_2_liver_2$gene_name)
pat_2_plasma_names <- names(pat_2_plasma$gene_name)

all_pat_2_names <- c(pat_2_breast_2_names, pat_2_liver_1_names, pat_2_liver_2_names, pat_2_plasma_names)

## put them together
all_mut_names <- c(all_pat_9_names, all_pat_ema_names, all_pat_10_names, all_pat_8_names, all_pat_2_names)

#plot it
all_freqs <- table(all_mut_names)
all_freqs <- all_freqs/24
all_freqs_ordered <- all_freqs[order(all_freqs, decreasing = TRUE)]
all_freqs_ordered <- round(all_freqs_ordered, digits = 2)
barplot(all_freqs_ordered, las = 2, ylim = c(0,1.0), cex.names = 0.8, col = 'dodgerblue')


tst_genes <- c('AKT1', 'BRIP1', 'CREBBP', 'FANCI', 'FGFR2', 'JAK3', 'MSH3', 'PALB2', 'RAD51D', 'TSC1',
'AKT2', 'BTK', 'CSF1R', 'FANCL', 'FGFR3', 'KDR', 'MSH6', 'PDGFRA', 'RAD54L', 'TSC2',
'AKT3', 'CARD11', 'CTNNB1', 'FBXW7', 'FGFR4', 'KIT', 'MTOR', 'PDGFRB', 'RB1', 'VHL',
'ALK', 'CCND1', 'DDR2', 'FGF1', 'FLT1', 'KMT2A', 'MLL', 'MUTYH', 'PIK3CA', 'RET', 'XRCC2',
'APC', 'CCND2', 'DNMT3A', 'FGF2', 'FLT3', 'KRAS', 'MYC', 'PIK3CB', 'RICTOR',
'AR', 'CCNE1', 'EGFR', 'FGF3', 'FOXL2', 'MAP2K1', 'MYCL1', 'PIK3CD', 'ROS1',
'ARID1A', 'CD79A', 'EP300', 'FGF4', 'GEN1', 'MAP2K2', 'MYCN', 'PIK3CG', 'RPS6KB1',
'ATM', 'CD79B', 'ERBB2', 'FGF5', 'GNA11', 'MCL1', 'MYD88', 'PIK3R1', 'SLX4',
'ATR', 'CDH1', 'ERBB3', 'FGF6', 'GNAQ', 'MDM2', 'NBN', 'PMS2', 'SMAD4',
'BAP1', 'CDK12', 'ERBB4', 'FGF7', 'GNAS', 'MDM4', 'NF1', 'PPP2R2A', 'SMARCB1',
'BARD1', 'CDK4', 'ERCC1', 'FGF8', 'HNF1A', 'MET', 'NOTCH1', 'PTCH1', 'SMO',
'BCL2', 'CDK6', 'ERCC2', 'FGF9', 'HRAS', 'MLH1', 'NOTCH2', 'PTEN', 'SRC',
'BCL6', 'CDKN2A', 'ERG', 'FGF10', 'IDH1', 'MLLT3', 'NOTCH3', 'PTPN11', 'STK11',
'BRAF', 'CEBPA', 'ESR1', 'FGF14', 'IDH2', 'MPL', 'NPM1', 'RAD51', 'TERT',
'BRCA1', 'CHEK1', 'EZH2', 'FGF23', 'INPP4B', 'MRE11A', 'NRAS', 'RAD51B', 'TET2',
'BRCA2', 'CHEK2', 'FAM175A', 'FGFR1', 'JAK2', 'MSH2', 'NRG1', 'RAD51C', 'TP53')

length(intersect(names(all_freqs_ordered), tst_genes))
names(all_freqs_ordered)[names(all_freqs_ordered) %!in% tst_genes]

# [1] "VAT1L"        "FGF19"        "CASC11"       "LAMP1"        "MIR6759"      "MYCL"         "CCND3"        "CD3EAP"       "MYCNOS"       "TFRC"         "RAF1"        
# [12] "TRIM5"        "GRTP1"        "HSPA12A"      "LOC100130075" "MIR641"       "MROH6"        "SDCCAG8"      "SDK1"         "TAF8"         "CSDE1"        "FGF10-AS1"   
# [23] "GNB3"         "KLLN"         "LOC401177"    "NUDT6"        "PTN"          "SYNE1"        "SYNE1-AS1"    "TOE1"         "UBAC2"



## missense/nonsense/silent ratios (from snpEff) ----

pat_9_ln_ratios <- c(51.465, 2.93, 45.604)
pat_9_oment_ratios <- c(66.812, 2.62, 30.568)
pat_9_ovary_ratios <- c(71.366, 1.322, 27.313)
pat_9_plasma_ratios <- c(70.207, 3.368, 26.425)

pat_ema_heart_ratios <- c(69.412, 1.176, 29.412)
pat_ema_l_kidney_ratios <- c(81.69, 0.00, 18.31)
pat_ema_r_kidney_ratios <- c(78.571, 0.714, 20.714)
pat_ema_liver_1_ratios <- c(77.703, 0.676, 21.622)
pat_ema_liver_2_ratios <- c(72.0, 1.0, 27.0)
pat_ema_oment_1_ratios <- c(82.009, 0.00, 17.901)
pat_ema_oment_2_ratios <- c(83.908, 0.00, 16.092)
pat_ema_plasma_ratios <- c(68.696, 5.843, 25.461)

pat_10_liver_1_ratios <- c(64.516, 0.00, 35.484)
pat_10_liver_2a_ratios <- c(70.968, 0.00, 29.032)
pat_10_liver_5_ratios <- c(65.263, 1.053, 33.684)
pat_10_plasma_ratios <- c(76.426, 3.57, 20.004)

pat_8_axillary_ratios <- c(77.841, 1.705, 20.455)
pat_8_breast_1_ratios <- c(70.984, 0.00, 29.016)
pat_8_breast_2_ratios <- c(72.34, 2.464, 25.196)
pat_8_plasma_ratios <- c(50.569, 0.00, 49.431)

pat_2_breast_2_ratios <- c(81.373, 1.961, 16.667)
pat_2_liver_1_ratios <- c(80.925, 1.734, 17.341)
pat_2_liver_2_ratios <- c(81.522, 3.804, 14.674)
pat_2_plasma_ratios <- c(56.164, 1.027, 42.808)



#combine and plot it
df <- rbind(pat_9_ln_ratios, pat_9_oment_ratios, pat_9_ovary_ratios, pat_9_plasma_ratios, pat_ema_heart_ratios, pat_ema_l_kidney_ratios, 
            pat_ema_r_kidney_ratios, pat_ema_liver_1_ratios, pat_ema_liver_2_ratios, pat_ema_oment_1_ratios, pat_ema_oment_2_ratios, 
            pat_ema_plasma_ratios, pat_10_liver_1_ratios, pat_10_liver_2a_ratios, pat_10_liver_5_ratios, pat_10_plasma_ratios, 
            pat_8_axillary_ratios, pat_8_breast_1_ratios, pat_8_breast_2_ratios, pat_8_plasma_ratios, pat_2_breast_2_ratios, pat_2_liver_1_ratios, 
            pat_2_liver_2_ratios, pat_2_plasma_ratios)
colnames(df) <- c('missense', 'nonsense', 'silent')
sample <- rownames(df)
df <- as.data.frame(df)
df$sample <- sample
df1 <- melt(df, id.vars = 'sample')
df1$sample <- factor(df1$sample)
library(ggplot2)
ggplot(df1, aes(x = sample, y = value, fill = variable)) + 
  geom_bar(stat = "identity")


## heatmap style plot ----

# function to get effects of genes
gene_effect <- function(df) {
  df_table <- data.frame(table(df$effect, df$gene_name))
  df_most <- df_table %>% group_by(Var2) %>% top_n(1, Freq) #most common mutation per gene
  df_most$Var2 <- as.character(df_most$Var2)
  df_most$Var1 <- as.character(df_most$Var1)
  df_first <- df_most[match(unique(df_most$Var2), df_most$Var2),] #get first listed if tie
  df_first <- df_first[, -3]
  colnames(df_first) <- c('effect', 'gene')
  return(df_first)
}
#act on pat 9
pat_9_ln_effects <- gene_effect(pat_9_ln)
pat_9_oment_effects <- gene_effect(pat_9_oment)
pat_9_ovary_effects <- gene_effect(pat_9_ovary)
pat_9_plasma_effects <- gene_effect(pat_9_plasma)

pat_ema_heart_effects <- gene_effect(pat_ema_heart)
pat_ema_l_kidney_effects <- gene_effect(pat_ema_l_kidney)
pat_ema_r_kidney_effects <- gene_effect(pat_ema_r_kidney)
pat_ema_liver_1_effects <- gene_effect(pat_ema_liver_1)
pat_ema_liver_2_effects <- gene_effect(pat_ema_liver_2)
pat_ema_oment_1_effects <- gene_effect(pat_ema_oment_1)
pat_ema_oment_2_effects <- gene_effect(pat_ema_oment_2)
pat_ema_plasma_effects <- gene_effect(pat_ema_plasma)

pat_10_liver_1_effects <- gene_effect(pat_10_liver_1)
pat_10_liver_2a_effects <- gene_effect(pat_10_liver_2a)
pat_10_liver_5_effects <- gene_effect(pat_10_liver_5)
pat_10_plasma_effects <- gene_effect(pat_10_plasma)

pat_8_axillary_effects <- gene_effect(pat_8_axillary)
pat_8_breast_1_effects <- gene_effect(pat_8_breast_1)
pat_8_breast_2_effects <- gene_effect(pat_8_breast_2)
pat_8_plasma_effects <- gene_effect(pat_8_plasma)

pat_2_breast_2_effects <- gene_effect(pat_2_breast_2)
pat_2_liver_1_effects <- gene_effect(pat_2_liver_1)
pat_2_liver_2_effects <- gene_effect(pat_2_liver_2)
pat_2_plasma_effects <- gene_effect(pat_2_plasma)

# get all 178 for commensurate plotting
all_freqs_df <- data.frame(all_freqs_ordered)
colnames(all_freqs_df)[1] <- 'gene'

#start merging
all_effects <- merge(all_freqs_df, pat_ema_plasma_effects, by = 'gene', all = TRUE)
colnames(all_effects)[3] <- 'pat_ema_plasma'
all_effects <- merge(all_effects, pat_ema_r_kidney_effects, by = 'gene', all = TRUE)
colnames(all_effects)[4] <- 'pat_ema_r_kidney'
all_effects <- merge(all_effects, pat_ema_oment_2_effects, by = 'gene', all = TRUE)
colnames(all_effects)[5] <- 'pat_ema_oment_2'
all_effects <- merge(all_effects, pat_ema_oment_1_effects, by = 'gene', all = TRUE)
colnames(all_effects)[6] <- 'pat_ema_oment_1'
all_effects <- merge(all_effects, pat_ema_liver_2_effects, by = 'gene', all = TRUE)
colnames(all_effects)[7] <- 'pat_ema_liver_2'
all_effects <- merge(all_effects, pat_ema_liver_1_effects, by = 'gene', all = TRUE)
colnames(all_effects)[8] <- 'pat_ema_liver_1'
all_effects <- merge(all_effects, pat_ema_l_kidney_effects, by = 'gene', all = TRUE)
colnames(all_effects)[9] <- 'pat_ema_l_kidney'
all_effects <- merge(all_effects, pat_ema_heart_effects, by = 'gene', all = TRUE)
colnames(all_effects)[10] <- 'pat_ema_heart'
all_effects <- merge(all_effects, pat_9_plasma_effects, by = 'gene', all = TRUE)
colnames(all_effects)[11] <- 'pat_9_plasma'
all_effects <- merge(all_effects, pat_9_ovary_effects, by = 'gene', all = TRUE)
colnames(all_effects)[12] <- 'pat_9_ovary'
all_effects <- merge(all_effects, pat_9_oment_effects, by = 'gene', all = TRUE)
colnames(all_effects)[13] <- 'pat_9_oment'
all_effects <- merge(all_effects, pat_9_ln_effects, by = 'gene', all = TRUE)
colnames(all_effects)[14] <- 'pat_9_ln'
all_effects <- merge(all_effects, pat_8_plasma_effects, by = 'gene', all = TRUE)
colnames(all_effects)[15] <- 'pat_8_plasma'
all_effects <- merge(all_effects, pat_8_breast_2_effects, by = 'gene', all = TRUE)
colnames(all_effects)[16] <- 'pat_8_breast_2'
all_effects <- merge(all_effects, pat_8_breast_1_effects, by = 'gene', all = TRUE)
colnames(all_effects)[17] <- 'pat_8_breast_1'
all_effects <- merge(all_effects, pat_8_axillary_effects, by = 'gene', all = TRUE)
colnames(all_effects)[18] <- 'pat_8_axillary'
all_effects <- merge(all_effects, pat_2_plasma_effects, by = 'gene', all = TRUE)
colnames(all_effects)[19] <- 'pat_2_plasma'
all_effects <- merge(all_effects, pat_2_liver_2_effects, by = 'gene', all = TRUE)
colnames(all_effects)[20] <- 'pat_2_liver_2'
all_effects <- merge(all_effects, pat_2_liver_1_effects, by = 'gene', all = TRUE)
colnames(all_effects)[21] <- 'pat_2_liver_1'
all_effects <- merge(all_effects, pat_2_breast_2_effects, by = 'gene', all = TRUE)
colnames(all_effects)[22] <- 'pat_2_breast_2'

# a more condensed way to do this?
#Reduce(function(x,y) merge(x = x, y = y, by = "Character"), 
       #list(height, gender, eyeColour))

rownames(all_effects) <- all_effects$gene
all_effects <- all_effects[, -c(1,2)]







heat_colors <- colorRampPalette(c('white', 'red', 'green', 'purple', 'blue', 'black', 'orange'))(n = 7)

# can probably get rid of some of these and lump into other

color_picker <- function(x) {
    if (is.na(x) == TRUE) {
      x <- 0
    }
    else if (x == '3_prime_UTR_variant') {
      x <- 0.1
    }
    else if (x == '5_prime_UTR_variant') {
      x <- 0.2
    }
    else if (x == 'missense_variant') {
      x <- 0.3
    }
    else if (x == 'structural_interaction_variant') {
      x <- 0.4
    }
    else if (x == 'stop_gained') {
      x <- 0.5
    }
    else if (x == 'synonymous_variant') {
      x <- 0.6
    }
    else {
      x <- 0.7
    }
  }

df2_copy <- as.matrix(df2)
df2_copy <- apply(all_effects, c(1,2), FUN = color_picker)
df2_copy <- matrix(unlist(df2_copy), ncol = 5, byrow = FALSE)
par(mar = c(5.1, 4.1, 13.1, 2.1))
heatmap.2(t(df2_copy), col = heat_colors, trace = 'none', Rowv = FALSE, Colv = FALSE, sepwidth=c(0.05,0.05),
          sepcolor="gray", colsep=1:nrow(df2_copy), rowsep=1:ncol(df2_copy))
par(mar = c(5.1, 4.1, 4.1, 2.1))

# pat_9_ln_first <- pat_9_ln[match(unique(pat_9_ln$gene_name), pat_9_ln$gene_name),]
# pat_9_ln_effects <- pat_9_ln_first$effect


## percent found figure ----
samples <- c('pat_9_ln', 'pat_9_oment', 'pat_9_ovary', 'pat_ema_heart', 'pat_ema_l_kidney', 'pat_ema_r_kidney', 'pat_ema_liver_1', 
             'pat_ema_liver_2', 'pat_ema_oment_1', 'pat_ema_oment_2', 'pat_10_liver_1', 'pat_10_liver_2a', 'pat_10_liver_5', 'pat_8_axillary', 
             'pat_8_breast_1', 'pat_8_breast_2', 'pat_2_breast_2', 'pat_2_liver_1', 'pat_2_liver_2')
number_muts <- c(154, 30, 15, 11, 7, 6, 4, 4, 4, 6, 14, 9, 10, 23, 35, 512, 23, 20, 27)
percent_found <- c(11/154, 11/30, 3/15, 4/11, 3/7, 2/6, 1/4, 2/4, 2/4, 3/6, 7/14, 5/9, 5/10, 7/23, 8/35, 42/512, 7/23, 4/20, 4/27)


mean_maf <- c(mean(pat_9_ln$AF), mean(pat_9_oment$AF), mean(pat_9_ovary$AF), mean(pat_ema_heart$AF), mean(pat_ema_l_kidney$AF), 
              mean(pat_ema_r_kidney$AF), mean(pat_ema_liver_1$AF), mean(pat_ema_liver_2$AF), mean(pat_ema_oment_1$AF), mean(pat_ema_oment_2$AF), 
              mean(pat_10_liver_1$AF), mean(pat_10_liver_2a$AF), mean(pat_10_liver_5$AF), mean(pat_8_axillary$AF), mean(pat_8_breast_1$AF), 
              mean(pat_8_breast_2$AF), mean(pat_2_breast_2$AF), mean(pat_2_liver_1$AF), mean(pat_2_liver_2$AF))
patient <- c(rep('patient_9', 3), rep('patient_ema', 7), rep('patient_10', 3), rep('patient_8', 3), rep('patient_2', 3))

found_df <- data.frame(number_muts, percent_found, mean_maf, patient)
found_df <- found_df[-16, ]
ggplot(found_df, aes(x = number_muts, y = percent_found)) +
  geom_point(aes(size = found_df$mean_maf, shape = as.factor(found_df$patient)))



#regression figures
# done witth ema upset df loaded
pat_ema_pl_h <- upset_df[upset_df$Heart == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_h <- rownames(pat_ema_pl_h)

pat_ema_heart_new <- pat_ema_heart[pat_ema_heart$location %in% pat_ema_pl_h, ]
pat_ema_heart_new <- pat_ema_heart_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_h, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]

pat_ema_pl_h_bind <- merge(pat_ema_plasma_new, pat_ema_heart_new, by = 'location', all = TRUE)
pat_ema_pl_h_bind <- pat_ema_pl_h_bind[, -1]
pat_ema_pl_h_bind[, 1] <- as.numeric(pat_ema_pl_h_bind[, 1])
pat_ema_pl_h_bind[, 2] <- as.numeric(pat_ema_pl_h_bind[, 2])
pat_ema_pl_h_lm <- lm(AF.x ~ ., data = pat_ema_pl_h_bind)

pat_ema_pl_lk <- upset_df[upset_df$L_Kidney == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_lk <- rownames(pat_ema_pl_lk)

pat_ema_l_kidney_new <- pat_ema_l_kidney[pat_ema_l_kidney$location %in% pat_ema_pl_lk, ]
pat_ema_l_kidney_new <- pat_ema_l_kidney_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_lk, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]

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

pat_ema_pl_l1_bind <- merge(pat_ema_plasma_new, pat_ema_liver_1_new, by = 'location', all = TRUE)
pat_ema_pl_l1_bind <- pat_ema_pl_l1_bind[, -1]
pat_ema_pl_l1_bind[, 1] <- as.numeric(pat_ema_pl_l1_bind[, 1])
pat_ema_pl_l1_bind[, 2] <- as.numeric(pat_ema_pl_l1_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_l1_bind)

pat_ema_l2 <- upset_df[upset_df$Liver_2 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_l2 <- rownames(pat_ema_l2)

pat_ema_liver_2_new <- pat_ema_liver_2[pat_ema_liver_2$location %in% pat_ema_l2, ]
pat_ema_liver_2_new <- pat_ema_liver_2_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_l2, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]

pat_ema_l2_bind <- merge(pat_ema_plasma_new, pat_ema_liver_2_new, by = 'location', all = TRUE)
pat_ema_l2_bind <- pat_ema_l2_bind[, -1]
pat_ema_l2_bind[, 1] <- as.numeric(pat_ema_l2_bind[, 1])
pat_ema_l2_bind[, 2] <- as.numeric(pat_ema_l2_bind[, 2])
lm(AF.x ~ ., data = pat_ema_l2_bind)

pat_ema_pl_om1 <- upset_df[upset_df$Oment_1 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_om1 <- rownames(pat_ema_pl_om1)

pat_ema_oment_1_new <- pat_ema_oment_1[pat_ema_oment_1$location %in% pat_ema_pl_om1, ]
pat_ema_oment_1_new <- pat_ema_oment_1_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_om1, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]

pat_ema_pl_om1_bind <- merge(pat_ema_plasma_new, pat_ema_oment_1_new, by = 'location', all = TRUE)
pat_ema_pl_om1_bind <- pat_ema_pl_om1_bind[, -1]
pat_ema_pl_om1_bind[, 1] <- as.numeric(pat_ema_pl_om1_bind[, 1])
pat_ema_pl_om1_bind[, 2] <- as.numeric(pat_ema_pl_om1_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_om1_bind)

pat_ema_pl_om2 <- upset_df[upset_df$Oment_2 == 1 & upset_df$Plasma == 1, ] #12
pat_ema_pl_om2 <- rownames(pat_ema_pl_om2)

pat_ema_oment_2_new <- pat_ema_oment_2[pat_ema_oment_2$location %in% pat_ema_pl_om2, ]
pat_ema_oment_2_new <- pat_ema_oment_2_new[, c('location', 'AF')]

pat_ema_plasma_new <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_pl_om2, ]
pat_ema_plasma_new <- pat_ema_plasma_new[, c('location', 'AF')]

pat_ema_pl_om2_bind <- merge(pat_ema_plasma_new, pat_ema_oment_2_new, by = 'location', all = TRUE)
pat_ema_pl_om2_bind <- pat_ema_pl_om2_bind[, -1]
pat_ema_pl_om2_bind[, 1] <- as.numeric(pat_ema_pl_om2_bind[, 1])
pat_ema_pl_om2_bind[, 2] <- as.numeric(pat_ema_pl_om2_bind[, 2])
lm(AF.x ~ ., data = pat_ema_pl_om2_bind)

# all kidney
pat_ema_all_kidney <- rbind(pat_ema_pl_lk_bind, pat_ema_pl_rk_bind)
lm(pat_ema_all_kidney$AF.x ~ pat_ema_all_kidney$AF.y)

# all liver
pat_ema_all_liver <- rbind(pat_ema_pl_l1_bind, pat_ema_l2_bind)
lm(pat_ema_all_liver$AF.x ~ pat_ema_all_liver$AF.y)

#need to add pat 9 to oment plus ovary and lymph THIS CAN ONLY BE DONE WITH PAT 9 UPSET DF IN ENVIRON!!!!!!!!!!!!!!!!!!!!!!!!!!
pat_9_pl_om <- upset_df[upset_df$Omental_Met == 1 & upset_df$Plasma == 1, ] #12
pat_9_pl_om <- rownames(pat_9_pl_om)

pat_9_oment_new <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_om, ]
pat_9_oment_new <- pat_9_oment_new[, c('location', 'AF')]

pat_9_plasma_new <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_om, ]
pat_9_plasma_new <- pat_9_plasma_new[, c('location', 'AF')]

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

pat_9_plasma_new <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov, ]
pat_9_plasma_new <- pat_9_plasma_new[, c('location', 'AF')]

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
# tumor size (cm^2)
ln_size <- 9.92
oment_size <- 10.64
ovary_size <- 7.92
kidney_size <- 2.07
liver_size <- 2.0
tumor_size <- c(ln_size, oment_size, ovary_size, kidney_size, liver_size)
mean_plasma <- c(0.3102, 0.2987, 0.2013, 0.2285, 0.7834)
tumor_graph_df <- data.frame(tumor_size, mean_plasma)
plot(tumor_graph_df$tumor_size, tumor_graph_df$mean_plasma, ylim = c(0, max(tumor_graph_df$mean_plasma)), xlim = c(0, max(tumor_graph_df$tumor_size)), 
     col = c('blue', 'orange', 'dodgerblue', 'purple', 'green'), pch = 16, cex = 1.2, xlab = expression(paste('Tumor Size (cm' ^ 2, ')')), ylab = 'Mean MAF in Plasma')
fit <- lm(tumor_graph_df$mean_plasma ~ tumor_graph_df$tumor_size)
abline(a = 0.54597, b = -0.02789)
summary(fit)
legend(x = 0.0, y = 0.75, legend = c('Kidney', 'Lymph', 'Liver', 'Omental', 'Ovary'), pch = 16, cex = 1.2, col = c('blue', 'green', 'orange', 'purple', 'dodgerblue'), bty = 'n')
legend(x = 6.5, y = 0.6, legend = expression(paste('R'^2, '= 0.08')), bty = 'n')

# perfusion
heart_perf <- 0.62 #260/417 0.62 ml/min/g
ln_perf <- 0.4 #0.4 ml/min/g
oment_perf <- 0.048 #0.048 ml/min/g
ovary_perf <- 0.000491 #STILL NEED
kidney_perf <- 3.03 #1325/438 3.03 ml/min/g
liver_perf <- 0.13 #423/3257 0.130 ml/min/g

heart_mean <- mean(pat_ema_pl_h_bind$AF.x)
ln_mean <- mean(pat_9_pl_ln_bind$AF.x)
oment_mean <- mean(pat_ema_all_oment$AF.x)
ovary_mean <- mean(pat_9_pl_ov_bind$AF.x)
kidney_mean <- mean(pat_ema_all_liver$AF.x)
liver_mean <- mean(pat_ema_all_liver$AF.x)

mean_plasma <- c(heart_mean, ln_mean, oment_mean, ovary_mean, kidney_mean, liver_mean)
perfusion <- c(heart_perf, ln_perf, oment_perf, ovary_perf, kidney_perf, liver_perf)

cor(perfusion, mean_plasma, method = 'spearman')

perf_graph_df <- data.frame(mean_plasma, perfusion)
plot(perf_graph_df$perfusion, perf_graph_df$mean_plasma, ylim = c(0,1.0), ylab = 'Mean MAF in Plasma', xlab = 'Relative Blood Perfusion of Organ')
fit2 <- lm(perf_graph_df$mean_plasma ~ perf_graph_df$perfusion)
fit2
abline(a = 0.8823, b = 0.4102)
summary(fit2)

# comparing pat 9 plasma w ema plasma
length(intersect(pat_9_plasma$location, pat_ema_plasma$location)) # only 16
#which ones?
pat_9_plasma_common <- pat_9_plasma[pat_9_plasma$location %in% intersect(pat_9_plasma$location, pat_ema_plasma$location), ]
table(pat_9_plasma_common$gene_name)
mean(pat_9_plasma_common$AF)

pat_9_plasma_not_common <- pat_9_plasma[pat_9_plasma$location %!in% intersect(pat_9_plasma$location, pat_ema_plasma$location), ]
mean(pat_9_plasma_not_common$AF)

pat_ema_plasma_not_common <- pat_ema_plasma[pat_ema_plasma$location %!in% intersect(pat_9_plasma$location, pat_ema_plasma$location), ]
mean(pat_ema_plasma_not_common$AF)

pat_9_locations <- pat_9_plasma_not_common$location
pat_ema_locations <- pat_ema_plasma_not_common$location



venn.diagram(list(Pat_9 = pat_9_plasma$location, Pat_9A = pat_ema_plasma$location), fill = c("dodgerblue", "purple"),
             alpha = c(0.5, 0.5), cex = 2, cat.fontface = 4, lty =1, fontfamily = 3, 
             filename = '../venn.png')
