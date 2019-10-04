source('Scripts/variant_fxn.R')
library(readxl)
library(ggsci)
library(UpSetR)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(gplots)
library(plotrix)
library(stringr)

## PATIENT EMA ----
pat_ema_heart <- read.delim('Data/Patient_EMA/pat_ema_heart_all_int_clean_hg19_ann.txt', header = TRUE, 
                            stringsAsFactors = FALSE, sep = '\t')
pat_ema_heart <- mutect_process(pat_ema_heart, sample_type = 'plasma') #40

pat_ema_l_kidney <- read.delim('Data/Patient_EMA/pat_ema_l_kidney_all_int_clean_hg19_ann.txt', header = TRUE, 
                               stringsAsFactors = FALSE, sep = '\t')
pat_ema_l_kidney <- mutect_process(pat_ema_l_kidney, sample_type = 'plasma') #36

pat_ema_r_kidney <- read.delim('Data/Patient_EMA/pat_ema_r_kidney_all_int_clean_hg19_ann.txt', header = TRUE, 
                               stringsAsFactors = FALSE, sep = '\t')
pat_ema_r_kidney <- mutect_process(pat_ema_r_kidney, sample_type = 'plasma') #44

pat_ema_liver_1 <- read.delim('Data/Patient_EMA/pat_ema_liver_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_liver_1 <- mutect_process(pat_ema_liver_1, sample_type = 'plasma') #38

pat_ema_liver_2 <- read.delim('Data/Patient_EMA/pat_ema_liver_2_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_liver_2 <- mutect_process(pat_ema_liver_2, sample_type = 'plasma') #32

pat_ema_oment_1 <- read.delim('Data/Patient_EMA/pat_ema_oment_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_oment_1 <- mutect_process(pat_ema_oment_1, sample_type = 'plasma') #41

pat_ema_oment_2 <- read.delim('Data/Patient_EMA/pat_ema_oment_2_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_oment_2 <- mutect_process(pat_ema_oment_2, sample_type = 'plasma') #43

pat_ema_plasma <- read.delim('Data/Patient_EMA/pat_ema_plasma_all_int_clean_hg19_ann.txt', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_ema_plasma <- mutect_process(pat_ema_plasma, sample_type = 'plasma') #1116

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_ema_heart$location, pat_ema_plasma$location)) #7/40
length(intersect(pat_ema_l_kidney$location, pat_ema_plasma$location)) #5/36
length(intersect(pat_ema_r_kidney$location, pat_ema_plasma$location)) #4/44
length(intersect(pat_ema_liver_1$location, pat_ema_plasma$location)) #2/38
length(intersect(pat_ema_liver_2$location, pat_ema_plasma$location)) #4/32
length(intersect(pat_ema_oment_1$location, pat_ema_plasma$location)) #4/41
length(intersect(pat_ema_oment_2$location, pat_ema_plasma$location)) #4/43

#pool mutations from 3 mets
pat_ema_met_pool <- unique(c(pat_ema_heart$location, pat_ema_l_kidney$location, pat_ema_r_kidney$location, 
                             pat_ema_liver_1$location, pat_ema_liver_2$location, pat_ema_oment_1$location, 
                             pat_ema_oment_2$location)) #72 unique mutations


pat_ema_plasma_found <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_met_pool, ] #4 mutations
pat_ema_plasma_found_vars <- (pat_ema_plasma_found$location) 

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

pat_ema_l_kidney_pooled <- pat_ema_l_kidney_pooled[, c('location', 'AF')]
colnames(pat_ema_l_kidney_pooled) <- c('location', 'AF_l_kidney')

pat_ema_r_kidney_pooled <- pat_ema_r_kidney_pooled[, c('location', 'AF')]
colnames(pat_ema_r_kidney_pooled) <- c('location', 'AF_r_kidney')

pat_ema_liver_1_pooled <- pat_ema_liver_1_pooled[, c('location', 'AF')]
colnames(pat_ema_liver_1_pooled) <- c('location', 'AF_liver_1')

pat_ema_liver_2_pooled <- pat_ema_liver_2_pooled[, c('location', 'AF')]
colnames(pat_ema_liver_2_pooled) <- c('location', 'AF_liver_2')

pat_ema_oment_1_pooled <- pat_ema_oment_1_pooled[, c('location', 'AF')]
colnames(pat_ema_oment_1_pooled) <- c('location', 'AF_oment_1')

pat_ema_oment_2_pooled <- pat_ema_oment_2_pooled[, c('location', 'AF')]
colnames(pat_ema_oment_2_pooled) <- c('location', 'AF_oment_2')

pat_ema_plasma_pooled <- pat_ema_plasma_pooled[, c('location', 'AF')]
colnames(pat_ema_plasma_pooled) <- c('location', 'AF_plasma')

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

pat_ema_7 <- c(pat_ema_muts_pooled$AF_heart, pat_ema_muts_pooled$AF_l_kidney, pat_ema_muts_pooled$AF_r_kidney, 
               pat_ema_muts_pooled$AF_liver_1, pat_ema_muts_pooled$AF_liver_2, pat_ema_muts_pooled$AF_oment_1, 
               pat_ema_muts_pooled$AF_oment_2)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_ema_muts_pooled)[1] * dim(pat_ema_muts_pooled)[2])
# plasma reds
cell_cols[50:56] <- color.scale(pat_ema_muts_pooled[, 8], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:49] <- color.scale(pat_ema_7, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 7, byrow = FALSE)
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
color.legend(0.5,-1.7,2.5,-0.9,round(c(min(pat_ema_muts_pooled[, 8], na.rm = TRUE), max(pat_ema_muts_pooled[, 8], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.2, at=0.75, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_ema_muts_pooled[, 1:7], na.rm = TRUE),max(pat_ema_muts_pooled[, 1:7], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(3,-1.7,5,-0.9,round(c(min(pat_ema_muts_pooled[, 1:7], na.rm = TRUE), max(pat_ema_muts_pooled[, 1:7], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.2, at=3.225, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 2.8, cex = 1.1, font = 2)

# add NA legend
color.legend(6.15, -1.7, 6.35, -0.9, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.55, at=5.7, cex = 1.1, font = 2)
legend(x=6.165,y=-0.83,legend='',pch=16,bty="n",xpd = NA)

#plot labels
#mut_col_labels <- rownames(pat_ema_muts_pooled)
mut_col_labels <- c('TP53\nc.-123C>G', 'FLT1\nc.*1999G>A', 'TSC2\np.D1734D', 'NRG1\np.M349T', 'MSH2\np.A305T', 'FLT3\np.D324N', 'FGFR2\nc.*303G>A')
axis(3, at = (1:ncol(pat_ema_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Heart', 'L Kidney', 'R Kidney', 'Liver 1', 'Liver 2', 'Omental 1', 'Omental 2')
axis(2, at = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2, hadj = 0.9)

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










#regression figures
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

