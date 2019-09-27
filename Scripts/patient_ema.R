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
pat_ema_heart <- mutect_process(pat_ema_heart) #50

pat_ema_l_kidney <- read.delim('Data/Patient_EMA/pat_ema_l_kidney_all_int_clean_hg19_ann.txt', header = TRUE, 
                               stringsAsFactors = FALSE, sep = '\t')
pat_ema_l_kidney <- mutect_process(pat_ema_l_kidney) #49

pat_ema_r_kidney <- read.delim('Data/Patient_EMA/pat_ema_r_kidney_all_int_clean_hg19_ann.txt', header = TRUE, 
                               stringsAsFactors = FALSE, sep = '\t')
pat_ema_r_kidney <- mutect_process(pat_ema_r_kidney) #59

pat_ema_liver_1 <- read.delim('Data/Patient_EMA/pat_ema_liver_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_liver_1 <- mutect_process(pat_ema_liver_1) #42

pat_ema_liver_2 <- read.delim('Data/Patient_EMA/pat_ema_liver_2_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_liver_2 <- mutect_process(pat_ema_liver_2) #41

pat_ema_oment_1 <- read.delim('Data/Patient_EMA/pat_ema_oment_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_oment_1 <- mutect_process(pat_ema_oment_1) #46

pat_ema_oment_2 <- read.delim('Data/Patient_EMA/pat_ema_oment_2_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_ema_oment_2 <- mutect_process(pat_ema_oment_2) #55

pat_ema_plasma <- read.delim('Data/Patient_EMA/pat_ema_plasma_all_int_clean_hg19_ann.txt', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_ema_plasma <- mutect_process(pat_ema_plasma) #1276

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_ema_heart$location, pat_ema_plasma$location)) #7/50
length(intersect(pat_ema_l_kidney$location, pat_ema_plasma$location)) #5/49
length(intersect(pat_ema_r_kidney$location, pat_ema_plasma$location)) #4/59
length(intersect(pat_ema_liver_1$location, pat_ema_plasma$location)) #2/42
length(intersect(pat_ema_liver_2$location, pat_ema_plasma$location)) #4/41
length(intersect(pat_ema_oment_1$location, pat_ema_plasma$location)) #4/46
length(intersect(pat_ema_oment_2$location, pat_ema_plasma$location)) #4/55

#pool mutations from 3 mets
pat_ema_met_pool <- unique(c(pat_ema_heart$location, pat_ema_l_kidney$location, pat_ema_r_kidney$location, 
                             pat_ema_liver_1$location, pat_ema_liver_2$location, pat_ema_oment_1$location, 
                             pat_ema_oment_2$location)) #97 unique mutations


pat_ema_plasma_found <- pat_ema_plasma[pat_ema_plasma$location %in% pat_ema_met_pool, ] #7 mutations
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
mut_col_labels <- c('', '', '', 'NRG1\np.M349T', 'MSH2\np.A305T', 'FLT3\np.D324N', '')
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
