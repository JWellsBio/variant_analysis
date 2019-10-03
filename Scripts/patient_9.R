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



## PATIENT 9 ----
## import data ----
pat_9_ln    <- read.delim('Data/Patient_9/pat_9_ln_all_int_clean_hg19_ann.txt', header = TRUE, 
                          stringsAsFactors = FALSE, sep = '\t')
pat_9_ln    <- mutect_process(pat_9_ln) #219

pat_9_oment <- read.delim('Data/Patient_9/pat_9_oment_all_int_clean_hg19_ann.txt', header = TRUE, 
                          stringsAsFactors = FALSE, sep = '\t')
pat_9_oment <- mutect_process(pat_9_oment) #124

pat_9_ovary <- read.delim('Data/Patient_9/pat_9_ovary_all_int_clean_hg19_ann.txt', header = TRUE, 
                          stringsAsFactors = FALSE, sep = '\t')
pat_9_ovary <- mutect_process(pat_9_ovary) #89

# how much in common between samples? ----
length(intersect(pat_9_ln$location, pat_9_oment$location)) #10
length(intersect(pat_9_ln$location, pat_9_ovary$location)) #6
length(intersect(pat_9_oment$location, pat_9_ovary$location)) #21

## plasma ----
pat_9_plasma <- read.delim('Data/Patient_9/pat_9_plasma_all_int_clean_hg19_ann.txt', header = TRUE, 
                           stringsAsFactors = FALSE, sep = '\t')
pat_9_plasma <- mutect_process(pat_9_plasma) #651

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_9_ln$location, pat_9_plasma$location)) #16/219
length(intersect(pat_9_oment$location, pat_9_plasma$location)) #30/124
length(intersect(pat_9_ovary$location, pat_9_plasma$location)) #28/89

#pool mutations from 3 mets
pat_9_met_pool <- unique(c(pat_9_ln$location, pat_9_oment$location, pat_9_ovary$location)) #399 unique mutations


pat_9_plasma_found <- pat_9_plasma[pat_9_plasma$location %in% pat_9_met_pool, ] #49 mutations
pat_9_plasma_found_vars <- (pat_9_plasma_found$location) 

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
pat_9_oment_pooled$AF_oment <- as.numeric(pat_9_oment_pooled$AF_oment)

pat_9_ovary_pooled <- pat_9_ovary_pooled[, c('location', 'AF')]
colnames(pat_9_ovary_pooled) <- c('location', 'AF_ovary')
pat_9_ovary_pooled$AF_ovary <- as.numeric(pat_9_ovary_pooled$AF_ovary)

pat_9_plasma_pooled <- pat_9_plasma_pooled[, c('location', 'AF')]
colnames(pat_9_plasma_pooled) <- c('location', 'AF_plasma')
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
cell_cols[148:196] <- color.scale(pat_9_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:147] <- color.scale(pat_9_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 49, byrow = FALSE)
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
mut_col_labels[1:15] <- c('', '', '', '', 'ERBB2 p.P1170A', '', '', '', '', '', '', '', 'NRG1 p.M349T', 'FGFR4 p.G388R', '')
mut_col_labels[16:25] <- c('', '', '', '', 'MSH2 p.A305T', '', 'BRIP1 p.S919P', 'TET2 p.L1721W', 'CCND3 p.S259A', 'MUTYH p.V22M')
mut_col_labels[26:35] <- c('FLT3 p.D324N', '', '', '', '', '', 'ERCC1 p.F40L', 'ERCC1 p.P38R', 'ERCC1 p.A36V', '')
mut_col_labels[36:45] <- c('FGFR1 p.L703F', '', '', '', 'FGFR1 p.D705Y', '', '', 'EGFR p.P1119H', '', 'EGFR p.P1123L')
mut_col_labels[46:49] <- c('EGFR p.A1118G', 'EGFR p.P1119A', 'EGFR p.S1120I', '')
axis(3, at = (1:ncol(pat_9_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.7, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Lymph\nMet   ', 'Omental\nMet   ', 'Ovary\nMet  ')
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
blue_pal <- blue_pal(8)
orange_pal <- pal_material('orange', n = 8, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(8)
purple_pal <- pal_material('purple', n= 4, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(4)
bar_colors <- c(blue_pal[c(1,3,4)], orange_pal[c(1:3, 5,5,5)], purple_pal[c(1,2,3,3)], 'darkgreen')

upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray60', Lymph_Met = 'white', Omental_Met = 'white', 
                                                                 Ovary_Met = 'white')))),
      intersections = list(list('Lymph_Met'), 
                           list('Ovary_Met'), 
                           list('Omental_Met'), 
                           list('Ovary_Met', 'Plasma'), 
                           list('Omental_Met', 'Plasma'), 
                           list('Lymph_Met', 'Plasma'), 
                           list('Omental_Met', 'Ovary_Met'), 
                           list('Omental_Met', 'Lymph_Met'), 
                           list('Ovary_Met', 'Lymph_Met'), 
                           list('Omental_Met', 'Ovary_Met', 'Plasma'), 
                           list('Omental_Met', 'Lymph_Met', 'Plasma'), 
                           list('Ovary_Met', 'Omental_Met', 'Lymph_Met'), 
                           list('Ovary_Met', 'Lymph_Met', 'Plasma'), 
                           list(colnames(upset_df))),
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray60', 'goldenrod4', 'aquamarine3', 'chocolate3'), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
      text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 2.0))


## looking at AF for each set of mutations ----
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
pat_9_ln_sd <- sd(pat_9_ln_af)
pat_9_ln_af <- median(pat_9_ln_af)

pat_9_om_only <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_om_only <- pat_9_om_only$mutation
pat_9_om_af <- pat_9_oment[pat_9_oment$location %in% pat_9_om_only, ]
pat_9_om_af <- pat_9_om_af$AF
pat_9_om_sd <- sd(pat_9_om_af)
pat_9_om_af <- median(pat_9_om_af)

pat_9_ov_only <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_ov_only <- pat_9_ov_only$mutation
pat_9_ov_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_only, ]
pat_9_ov_af <- pat_9_ov_af$AF
pat_9_ov_sd <- sd(pat_9_ov_af)
pat_9_ov_af <- median(pat_9_ov_af)

pat_9_pl_ov <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov <- pat_9_pl_ov$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_pl_ov_af <- c(pat_9_ovary_met_af, pat_9_plasma_af)
pat_9_pl_ov_sd <- sd(pat_9_pl_ov_af)
pat_9_pl_ov_af <- median(pat_9_pl_ov_af)

pat_9_pl_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_pl_om <- pat_9_pl_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_pl_om_af <- c(pat_9_oment_met_af, pat_9_plasma_af)
pat_9_pl_om_sd <- sd(pat_9_pl_om_af)
pat_9_pl_om_af <- median(pat_9_pl_om_af)

pat_9_pl_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_pl_ly <- pat_9_pl_ly$mutation
pat_9_lymph_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ly, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ly, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_pl_ly_af <- c(pat_9_lymph_met_af, pat_9_plasma_af)
pat_9_pl_ly_sd <- sd(pat_9_pl_ly_af)
pat_9_pl_ly_af <- median(pat_9_pl_ly_af)

pat_9_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_ov_om <- pat_9_ov_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ov_om_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af)
pat_9_ov_om_sd <- sd(pat_9_ov_om_af)
pat_9_ov_om_af <- median(pat_9_ov_om_af)

pat_9_om_ly <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_om_ly <- pat_9_om_ly$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_om_ly, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ln_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_om_ly, ]
pat_9_ln_met_af <- pat_9_ln_met_af$AF
pat_9_om_ly_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af)
pat_9_om_ly_sd <- sd(pat_9_om_ly_af)
pat_9_om_ly_af <- median(pat_9_om_ly_af)

pat_9_ov_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
pat_9_ov_ly <- pat_9_ov_ly$mutation
pat_9_ln_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_ov_ly, ]
pat_9_ln_met_af <- pat_9_ln_met_af$AF
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_ly, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ov_ly_af <- c(pat_9_ln_met_af, pat_9_ovary_met_af)
pat_9_ov_ly_sd <- sd(pat_9_ov_ly_af)
pat_9_ov_ly_af <- median(pat_9_ov_ly_af)

pat_9_pl_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov_om <- pat_9_pl_ov_om$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_pl_ov_om_af <- c(pat_9_ovary_met_af, pat_9_plasma_af, pat_9_oment_met_af)
pat_9_pl_ov_om_sd <- sd(pat_9_pl_ov_om_af)
pat_9_pl_ov_om_af <- median(pat_9_pl_ov_om_af)

pat_9_pl_ly_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_pl_ly_om <- pat_9_pl_ly_om$mutation
pat_9_ln_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ly_om, ]
pat_9_ln_met_af <- pat_9_ln_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ly_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_ly_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_pl_ly_om_af <- c(pat_9_ln_met_af, pat_9_plasma_af, pat_9_oment_met_af)
pat_9_pl_ly_om_sd <- sd(pat_9_pl_ly_om_af)
pat_9_pl_ly_om_af <- median(pat_9_pl_ly_om_af)

pat_9_ly_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
pat_9_ly_ov_om <- pat_9_ly_ov_om$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ly_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_ln_af <- pat_9_ln[pat_9_ln$location %in% pat_9_ly_ov_om, ]
pat_9_ln_af <- pat_9_ln_af$AF
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_ly_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_ly_ov_om_af <- c(pat_9_ovary_met_af, pat_9_ln_af, pat_9_oment_met_af)
pat_9_ly_ov_om_sd <- sd(pat_9_ly_ov_om_af)
pat_9_ly_ov_om_af <- median(pat_9_ly_ov_om_af)

pat_9_pl_ov_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
pat_9_pl_ov_ly <- pat_9_pl_ov_ly$mutation
pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov_ly, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov_ly, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_ln_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ov_ly, ]
pat_9_ln_met_af <- pat_9_ln_met_af$AF
pat_9_pl_ov_ly_af <- c(pat_9_ovary_met_af, pat_9_plasma_af, pat_9_ln_met_af)
pat_9_pl_ov_ly_sd <- sd(pat_9_pl_ov_ly_af)
pat_9_pl_ov_ly_af <- median(pat_9_pl_ov_ly_af)

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
pat_9_all_sd <- sd(pat_9_all_af)
pat_9_all_af <- median(pat_9_all_af)


pat_9_afs <- c(pat_9_ln_af, pat_9_om_af, pat_9_ov_af, pat_9_pl_ov_af, pat_9_pl_om_af, pat_9_pl_ly_af,  
               pat_9_ov_om_af, pat_9_om_ly_af, pat_9_ov_ly_af, pat_9_pl_ov_om_af, pat_9_pl_ly_om_af, pat_9_ly_ov_om_af, 
               pat_9_pl_ov_ly_af, pat_9_all_af)
pat_9_sds <- c(pat_9_ln_sd, pat_9_om_sd, pat_9_ov_sd, pat_9_pl_ov_sd, pat_9_pl_om_sd, pat_9_pl_ly_sd,  
               pat_9_ov_om_sd, pat_9_om_ly_sd, pat_9_ov_ly_sd, pat_9_pl_ov_om_sd, pat_9_pl_ly_om_sd, pat_9_ly_ov_om_sd, 
               pat_9_pl_ov_ly_sd, pat_9_all_sd)
pat_9_labels <- c('Lymph\nMet', 'Omental\nMet', 'Ovary\nMet', 'Ovary Met\n+ Plasma', 'Omental Met\n+ Plasma', 
                  'Lymph Met\n+ Plasma', 'Omental Met\n+ Ovary Met', 'Omental Met\n+ Lymph Met', 
                  'Ovary Met\n+ Lymph Met', 'Omental Met\n+ Ovary Met\n+ Plasma', 'Lymph Met\n+ Omental Met\n+ Plasma', 
                  'Lymph Met\n+ Ovary Met\n+ Omental Met', 'Lymph Met\n+ Ovary Met\n+ Plasma', 'All Samples')
pat_9_afs <- data.frame(pat_9_labels, pat_9_afs, pat_9_sds, bar_colors)
p<-ggplot(data=pat_9_afs, aes(x=pat_9_labels, y=pat_9_afs, fill = bar_colors)) +
  geom_bar(stat="identity", width = 0.5, color = 'black') + scale_x_discrete(limits=pat_9_afs$pat_9_labels) + ylab('Median Mutant Allele Frequency') +
  xlab('Intersections') + theme_bw() + 
  theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) + 
  scale_fill_manual("legend", values = bar_colors) +
  geom_errorbar(aes(ymin=pat_9_afs, ymax=pat_9_afs+pat_9_sds), width=.2,
                position=position_dodge(.9)) + theme(axis.text=element_text(size=8, face = 'bold'),
                                                     axis.title=element_text(size=16,face="bold"))

p
###ADJUST COLORS AND FINALIZE THIS!!!!!!!!!!!!!!!








