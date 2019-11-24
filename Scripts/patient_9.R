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
#lymph node met
pat_9_ln    <- read.delim('Data/Patient_9/pat_9_ln_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                          stringsAsFactors = FALSE, sep = '\t')
pat_9_ln    <- mutect_process(pat_9_ln) #32

#omental met
pat_9_oment <- read.delim('Data/Patient_9/pat_9_oment_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                          stringsAsFactors = FALSE, sep = '\t')
pat_9_oment <- mutect_process(pat_9_oment) #5

#ovarian met
pat_9_ovary <- read.delim('Data/Patient_9/pat_9_ovary_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                          stringsAsFactors = FALSE, sep = '\t')
pat_9_ovary <- mutect_process(pat_9_ovary) #3


## plasma ----
pat_9_plasma <- read.delim('Data/Patient_9/pat_9_plasma_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                           stringsAsFactors = FALSE, sep = '\t')
pat_9_plasma <- mutect_process(pat_9_plasma, sample_type = 'plasma') #98



## looking at how well plasma detects tumor mutations ----
length(intersect(pat_9_ln$location, pat_9_plasma$location)) #9/136
length(intersect(pat_9_oment$location, pat_9_plasma$location)) #10/24
length(intersect(pat_9_ovary$location, pat_9_plasma$location)) #1/9

#pool mutations from 3 mets
pat_9_met_pool <- unique(c(pat_9_ln$location, pat_9_oment$location, pat_9_ovary$location)) #158 unique mutations


pat_9_plasma_found <- pat_9_plasma[pat_9_plasma$location %in% pat_9_met_pool, ] #13 mutations
pat_9_plasma_found_vars <- (pat_9_plasma_found$location) 

# subset to pooled plasma and take a look
pat_9_ln_pooled <- pat_9_ln[pat_9_ln$location %in% pat_9_plasma_found_vars, ]

pat_9_oment_pooled <- pat_9_oment[pat_9_oment$location %in% pat_9_plasma_found_vars, ]

pat_9_ovary_pooled <- pat_9_ovary[pat_9_ovary$location %in% pat_9_plasma_found_vars, ]

pat_9_plasma_pooled <- pat_9_plasma[pat_9_plasma$location %in% pat_9_plasma_found_vars, ]


# subset these to tumor AF for new figure
pat_9_ln_pooled <- pat_9_ln_pooled[, c('location', 'AF')]
colnames(pat_9_ln_pooled) <- c('location', 'AF_ln')

pat_9_oment_pooled <- pat_9_oment_pooled[, c('location', 'AF')]
colnames(pat_9_oment_pooled) <- c('location', 'AF_oment')

pat_9_ovary_pooled <- pat_9_ovary_pooled[, c('location', 'AF')]
colnames(pat_9_ovary_pooled) <- c('location', 'AF_ovary')

pat_9_plasma_pooled <- pat_9_plasma_pooled[, c('location', 'AF')]
colnames(pat_9_plasma_pooled) <- c('location', 'AF_plasma')

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
cell_cols[40:52] <- color.scale(pat_9_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:39] <- color.scale(pat_9_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 13, byrow = FALSE)
pat_9_pooled_t <- data.frame(t(pat_9_muts_pooled))
pat_9_pooled_t <- pat_9_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]

# plot it
# extra space
par(mar=c(6,28.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_9_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_9_muts_pooled[, 4], na.rm = TRUE),max(pat_9_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(0.5,-0.9,5,-0.5,round(c(min(pat_9_muts_pooled[, 4], na.rm = TRUE), max(pat_9_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=1.25, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE),max(pat_9_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(5.5,-0.9,10,-0.5,round(c(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE), max(pat_9_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=6.15, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 5.4, cex = 1.1, font = 2)

# add NA legend
color.legend(12, -0.9, 12.75, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\nNot\nPresent', side=1, line= 3.2, at=11.1, cex = 0.8, font = 2)
legend(x=12.07,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_9_muts_pooled)
mut_col_labels[1:5] <- c('MTOR\np.N999N', 'VAT1L\nc*5485A>G', 'ERBB2\np.P1170A', 'ALK\np.P234P', 'ESR1\np.R245R')
mut_col_labels[6:10] <- c('EGFR\np.T903T', 'TP53\nc.-123C>G', 'NOTCH3\np.P1521P', 'FLT1\nc.*1999G>A', 'NOTCH1\np.D1698D')
mut_col_labels[11:13] <- c('FLT3\np.D324N', 'BRIP1\np.Y1137Y', 'FGFR2\nc.*303G>A')
axis(3, at = (1:ncol(pat_9_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Lymph\nMet   ', 'Omental\nMet   ', 'Ovary\nMet  ')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

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

# another way to build df?
location_list <- list(pat_9_ln$location, pat_9_oment$location, pat_9_ovary$location, pat_9_plasma$location)
upset_df <- fromList(location_list)
colnames(upset_df) <- c('Lymph_Met', 'Omental_Met', 'Ovary_Met', 'Plasma')
#head(upset_df)
head(upset_df)

#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:4])
randomnumber <- round(runif(4, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 8, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(8)
orange_pal <- pal_material('orange', n = 8, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(8)
purple_pal <- pal_material('purple', n= 8, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(8)
bar_colors <- c(blue_pal[c(1,4,5)], orange_pal[c(1,3,3,3)], purple_pal[c(1,2,4)])

upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray60', Lymph_Met = 'white', Omental_Met = 'white', 
                                                                 Ovary_Met = 'white')))),
      
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = TRUE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray60', 'goldenrod4', 'aquamarine3', 'chocolate3'), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      mainbar.y.label = 'Number of Mutations\nin Common', 
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

pat_9_pl_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_pl_om <- pat_9_pl_om$mutation
pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_pl_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$AF
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_om, ]
pat_9_plasma_af <- pat_9_plasma_af$AF
pat_9_pl_om_af <- c(pat_9_oment_met_af, pat_9_plasma_af)
pat_9_pl_om_sd <- sd(pat_9_pl_om_af)
pat_9_pl_om_af <- median(pat_9_pl_om_af)

# pat_9_pl_ov <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
# pat_9_pl_ov <- pat_9_pl_ov$mutation
# pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov, ]
# pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
# pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov, ]
# pat_9_plasma_af <- pat_9_plasma_af$AF
# pat_9_pl_ov_af <- c(pat_9_ovary_met_af, pat_9_plasma_af)
# pat_9_pl_ov_sd <- sd(pat_9_pl_ov_af)
# pat_9_pl_ov_af <- median(pat_9_pl_ov_af)

# pat_9_ov_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
# pat_9_ov_ly <- pat_9_ov_ly$mutation
# pat_9_ln_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_ov_ly, ]
# pat_9_ln_met_af <- pat_9_ln_met_af$AF
# pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_ov_ly, ]
# pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
# pat_9_ov_ly_af <- c(pat_9_ln_met_af, pat_9_ovary_met_af)
# pat_9_ov_ly_sd <- sd(pat_9_ov_ly_af)
# pat_9_ov_ly_af <- median(pat_9_ov_ly_af)

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

# pat_9_pl_ov_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
# pat_9_pl_ov_ly <- pat_9_pl_ov_ly$mutation
# pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_pl_ov_ly, ]
# pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
# pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_pl_ov_ly, ]
# pat_9_plasma_af <- pat_9_plasma_af$AF
# pat_9_ln_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_pl_ov_ly, ]
# pat_9_ln_met_af <- pat_9_ln_met_af$AF
# pat_9_pl_ov_ly_af <- c(pat_9_ovary_met_af, pat_9_plasma_af, pat_9_ln_met_af)
# pat_9_pl_ov_ly_sd <- sd(pat_9_pl_ov_ly_af)
# pat_9_pl_ov_ly_af <- median(pat_9_pl_ov_ly_af)
# 
# pat_9_all <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
# pat_9_all <- pat_9_all$mutation
# pat_9_oment_met_af <- pat_9_oment[pat_9_oment$location %in% pat_9_all, ]
# pat_9_oment_met_af <- pat_9_oment_met_af$AF
# pat_9_ovary_met_af <- pat_9_ovary[pat_9_ovary$location %in% pat_9_all, ]
# pat_9_ovary_met_af <- pat_9_ovary_met_af$AF
# pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$location %in% pat_9_all, ]
# pat_9_plasma_af <- pat_9_plasma_af$AF
# pat_9_lymph_met_af <- pat_9_ln[pat_9_ln$location %in% pat_9_all, ]
# pat_9_lymph_met_af <- pat_9_lymph_met_af$AF
# pat_9_all_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af, pat_9_plasma_af, pat_9_lymph_met_af)
# pat_9_all_sd <- sd(pat_9_all_af)
# pat_9_all_af <- median(pat_9_all_af)


pat_9_afs <- c(pat_9_ln_af, pat_9_om_af, pat_9_ov_af, pat_9_pl_ly_af, pat_9_ov_om_af, pat_9_om_ly_af, 
               pat_9_pl_om_af, pat_9_pl_ly_om_af, pat_9_pl_ov_om_af, pat_9_ly_ov_om_af)
pat_9_sds <- c(pat_9_ln_sd, pat_9_om_sd, pat_9_ov_sd, pat_9_pl_ly_sd, pat_9_ov_om_sd, pat_9_om_ly_sd, 
               pat_9_pl_om_sd, pat_9_pl_ly_om_sd, pat_9_pl_ov_om_sd, pat_9_ly_ov_om_sd)
pat_9_labels <- c('Lymph\nMet', 'Omental\nMet', 'Ovary\nMet', 'Lymph Met\n+ Plasma', 'Omental Met\n+ Ovary Met', 'Omental Met\n+ Lymph Met', 
                  'Omental Met\n+ Plasma', 'Lymph Met\n+ Omental Met\n+ Plasma', 'Omental Met\n+ Ovary Met\n+ Plasma', 
                  'Lymph Met\n+ Ovary Met\n+ Omental Met')
pat_9_afs <- data.frame(pat_9_labels, pat_9_afs, pat_9_sds, bar_colors)
p<-ggplot(data=pat_9_afs, aes(x=pat_9_labels, y=pat_9_afs, fill = bar_colors)) +
  geom_bar(stat="identity", width = 0.5, color = 'black') + scale_x_discrete(limits=pat_9_afs$pat_9_labels) + ylab('Median Mutant Allele Frequency') +
  xlab('Intersections') + theme_bw() + 
  theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) + 
  scale_fill_manual("legend", values = bar_colors[c(1:3, 8, 9, 10, 4, 5, 6, 7)]) +
  geom_errorbar(aes(ymin=pat_9_afs, ymax=pat_9_afs+pat_9_sds), width=.2,
                position=position_dodge(.9)) + theme(axis.text=element_text(size=8, face = 'bold'),
                                                     axis.title=element_text(size=16,face="bold"))

p



## red/black figures to look at how far down plasma can detect ----

pat_9_lymph_met_stats <- pat_9_ln[, c('location', 'AF')]
pat_9_oment_met_stats <- pat_9_oment[, c('location', 'AF')]
pat_9_ovary_met_stats <- pat_9_ovary[, c('location', 'AF')]
pat_9_all <- rbind(pat_9_lymph_met_stats, pat_9_oment_met_stats)
pat_9_all <- rbind(pat_9_all, pat_9_ovary_met_stats)
pat_9_all$AF <- as.numeric(pat_9_all$AF)
pat_9_all <- pat_9_all[order(pat_9_all$AF, decreasing = TRUE), ]
pat_9_all$color <- ifelse(pat_9_all$location %in% pat_9_plasma$location, 'red', 'black')
barplot(pat_9_all$AF, col = pat_9_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 9\n(3 tumors)', ylim = c(0,1.0))


pat_9_plasma_met_stats <- pat_9_plasma[, c('location', 'AF')]
pat_9_plasma_met_stats$AF <- as.numeric(pat_9_plasma_met_stats$AF)
pat_9_plasma_met_stats <- pat_9_plasma_met_stats[order(pat_9_plasma_met_stats$AF, decreasing = TRUE), ]
pat_9_plasma_met_stats$color <- ifelse(pat_9_plasma_met_stats$location %in% pat_9_all$location, 'red', 'black')
barplot(pat_9_plasma_met_stats$AF, col = pat_9_plasma_met_stats$color, ylab = 'Mutant Allele Frequency in Plasma', 
        xlab = 'Variant', main = 'Patient 9\n(3 tumors)', ylim = c(0,1.0))


