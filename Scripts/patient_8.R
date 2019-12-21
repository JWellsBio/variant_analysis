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

## PATIENT 8 ----
pat_8_axillary <- read.delim('Data/Patient_8/pat_8_axillary_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_8_axillary <- mutect_process(pat_8_axillary) #17

pat_8_breast_1 <- read.delim('Data/Patient_8/pat_8_breast_1_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_8_breast_1 <- mutect_process(pat_8_breast_1) #29

pat_8_breast_2 <- read.delim('Data/Patient_8/pat_8_breast_2_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t') 
pat_8_breast_2 <- mutect_process(pat_8_breast_2) #85
pat_8_breast_2 <- pat_8_breast_2[!grepl('clustered_events', pat_8_breast_2$FILTER), ]

pat_8_plasma <- read.delim('Data/Patient_8/pat_8_plasma_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                           stringsAsFactors = FALSE, sep = '\t')
pat_8_plasma <- mutect_process(pat_8_plasma, sample_type = 'plasma') #137 to 139 at 0.01

#all in common?
Reduce(intersect, list(pat_8_axillary$location, pat_8_breast_1$location, pat_8_breast_2$location)) #8 found NONE

## looking at how well plasma detects tumor mutations ----
length(intersect(pat_8_axillary$location, pat_8_plasma$location)) #7/17 41.2%
length(intersect(pat_8_breast_1$location, pat_8_plasma$location)) #8/29 27.6%
length(intersect(pat_8_breast_2$location, pat_8_plasma$location)) #16/85 18.8%


#pool mutations from 3 mets
pat_8_met_pool <- unique(c(pat_8_axillary$location, pat_8_breast_1$location, pat_8_breast_2$location)) #101 unique mutations


pat_8_plasma_found <- pat_8_plasma[pat_8_plasma$location %in% pat_8_met_pool, ] #23 mutations no add'l at 0.01
pat_8_plasma_found_vars <- (pat_8_plasma_found$location)


# subset to pooled plasma and take a look
pat_8_axillary_pooled <- pat_8_axillary[pat_8_axillary$location %in% pat_8_plasma_found_vars, ]

pat_8_breast_1_pooled <- pat_8_breast_1[pat_8_breast_1$location %in% pat_8_plasma_found_vars, ]

pat_8_breast_2_pooled <- pat_8_breast_2[pat_8_breast_2$location %in% pat_8_plasma_found_vars, ]

pat_8_plasma_pooled <- pat_8_plasma[pat_8_plasma$location %in% pat_8_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_8_axillary_pooled <- pat_8_axillary_pooled[, c('location', 'AF')]
colnames(pat_8_axillary_pooled) <- c('location', 'AF_axillary')

pat_8_breast_1_pooled <- pat_8_breast_1_pooled[, c('location', 'AF')]
colnames(pat_8_breast_1_pooled) <- c('location', 'AF_breast_1')

pat_8_breast_2_pooled <- pat_8_breast_2_pooled[, c('location', 'AF')]
colnames(pat_8_breast_2_pooled) <- c('location', 'AF_breast_2')

pat_8_plasma_pooled <- pat_8_plasma_pooled[, c('location', 'AF')]
colnames(pat_8_plasma_pooled) <- c('location', 'AF_plasma')

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
cell_cols[70:92] <- color.scale(pat_8_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:69] <- color.scale(pat_8_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 23, byrow = FALSE)
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
color.legend(1,-0.9,12,-0.5,round(c(min(pat_8_muts_pooled[, 4], na.rm = TRUE), max(pat_8_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=2.5, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_8_muts_pooled[, 1:3], na.rm = TRUE),max(pat_8_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(14,-0.9,26,-0.5,round(c(min(pat_8_muts_pooled[, 1:3], na.rm = TRUE), max(pat_8_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=15.5, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 12.5, cex = 1.1, font = 2)

# add NA legend
color.legend(19, -0.9, 20, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.8, at=17.7, cex = 1.1, font = 2)
legend(x=19.2,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_8_muts_pooled)
mut_col_labels[1:5] <- c('RET\np.L769L', 'LAMP1\np.R186R', 'HSPA12A\nc.-4919C>T', 'ESR1\np.R245R', 'SLX4\np.N457K')
mut_col_labels[6:10] <- c('HNF1A\np.G288G', 'MAP2K2\np.I220I', 'FGFR4\np.R54R', 'FLT1\nc.*265A>G', 'SLX4\np.A952V')
mut_col_labels[11:15] <- c('SLX4\np.A952T', 'ERG\nc.-64C>T', 'FGFR2\nc.*190A>G', 'ROS1\np.L101L', 'FLT1\nc.*3370T>G')
mut_col_labels[16:20] <- c('FANCI\np.C742S', 'CDH1\np.A692A', 'MUTYH\np.Q338H', 'NOTCH3\np.R1857W', 'CHEK2\nc.1117A>G')
mut_col_labels[21:23] <- c('CHEK2\nc.1116C>T', 'ERBB3\np.S1119C', 'ROS1\np.R167Q')
axis(3, at = (1:ncol(pat_8_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Axillary\n(7/17)', 'Breast 1\n(8/29)', 'Breast 2\n(16/85)')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# axillary points
points(x = which(is.na(pat_8_muts_pooled$AF_axillary)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_8_muts_pooled$AF_axillary))), 
       pch = 16)
# breast 1 points
points(x = which(is.na(pat_8_muts_pooled$AF_breast_1)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_8_muts_pooled$AF_breast_1))), 
       pch = 16)
# breast 2 points
points(x = which(is.na(pat_8_muts_pooled$AF_breast_2)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_8_muts_pooled$AF_breast_2))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))


##upset
pat_8_axillary_met_muts <- pat_8_axillary$location
pat_8_axillary_met_muts <- as.data.frame(pat_8_axillary_met_muts)
pat_8_axillary_met_muts$mutation <- rep(1, nrow(pat_8_axillary_met_muts))
colnames(pat_8_axillary_met_muts) <- c('mutation', 'Axillary_Met')

pat_8_breast_1_met_muts <- pat_8_breast_1$location
pat_8_breast_1_met_muts <- as.data.frame(pat_8_breast_1_met_muts)
pat_8_breast_1_met_muts$mutation <- rep(1, nrow(pat_8_breast_1_met_muts))
colnames(pat_8_breast_1_met_muts) <- c('mutation', 'Breast_1_Met')

pat_8_breast_2_met_muts <- pat_8_breast_2$location
pat_8_breast_2_met_muts <- as.data.frame(pat_8_breast_2_met_muts)
pat_8_breast_2_met_muts$mutation <- rep(1, nrow(pat_8_breast_2_met_muts))
colnames(pat_8_breast_2_met_muts) <- c('mutation', 'Breast_2_Met')

pat_8_plasma_muts <- pat_8_plasma$location
pat_8_plasma_muts <- as.data.frame(pat_8_plasma_muts)
pat_8_plasma_muts$mutation <- rep(1, nrow(pat_8_plasma_muts))
colnames(pat_8_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_8_axillary_met_muts, pat_8_breast_1_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_8_breast_2_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_8_plasma_muts, by = 'mutation', all = TRUE)
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

blue_pal <- pal_material('blue', n = 8, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(8)
orange_pal <- pal_material('orange', n = 8, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(8)
purple_pal <- pal_material('purple', n= 8, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(8)
bar_colors <- c(blue_pal[c(1,4,5)], orange_pal[c(1,3,3,3)], purple_pal[c(1,2,4)])

upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray60', Axillary_Met = 'white', Breast_1_Met = 'white', 
                                                                 Breast_2_Met = 'white')))),
      
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray60', 'goldenrod4', 'aquamarine3', 'chocolate3'), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      mainbar.y.label = 'Number of Mutations\nin Common', 
      text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 2.0))





#red/black figures to look at how far down plasma can detect

pat_8_axillary_met_stats <- pat_8_axillary[, c('location', 'AF')]
pat_8_breast_1_met_stats <- pat_8_breast_1[, c('location', 'AF')]
pat_8_breast_2_met_stats <- pat_8_breast_2[, c('location', 'AF')]
pat_8_all <- rbind(pat_8_axillary_met_stats, pat_8_breast_1_met_stats)
pat_8_all <- rbind(pat_8_all, pat_8_breast_2_met_stats)
pat_8_all$AF <- as.numeric(pat_8_all$AF)
pat_8_all <- pat_8_all[order(pat_8_all$AF, decreasing = TRUE), ]
pat_8_all$color <- ifelse(pat_8_all$location %in% pat_8_plasma$location, 'red', 'black')
barplot(pat_8_all$AF, col = pat_8_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 8\n(3 tumors)', ylim = c(0,1.0))

pat_8_plasma_met_stats <- pat_8_plasma[, c('location', 'AF')]
pat_8_plasma_met_stats$AF <- as.numeric(pat_8_plasma_met_stats$AF)
pat_8_plasma_met_stats <- pat_8_plasma_met_stats[order(pat_8_plasma_met_stats$AF, decreasing = TRUE), ]
pat_8_plasma_met_stats$color <- ifelse(pat_8_plasma_met_stats$location %in% pat_8_all$location, 'red', 'black')
barplot(pat_8_plasma_met_stats$AF, col = pat_8_plasma_met_stats$color, ylab = 'Mutant Allele Frequency in Plasma', 
        xlab = 'Variant', main = 'Patient 8\n(3 tumors)', ylim = c(0,1.0))


## correlations ----
pat_8_axillary_plasma <- intersect(pat_8_axillary$location, pat_8_plasma$location)
axillary_corr <- pat_8_axillary[pat_8_axillary$location %in% pat_8_axillary_plasma, ]
plasma_corr <- pat_8_plasma[pat_8_plasma$location %in% pat_8_axillary_plasma, ]
cor.test(axillary_corr$AF, plasma_corr$AF, method = 'spearman') #cor = 0.28 p = 0.55

pat_8_breast_1_plasma <- intersect(pat_8_breast_1$location, pat_8_plasma$location)
breast_1_corr <- pat_8_breast_1[pat_8_breast_1$location %in% pat_8_breast_1_plasma, ]
plasma_corr <- pat_8_plasma[pat_8_plasma$location %in% pat_8_breast_1_plasma, ]
cor.test(breast_1_corr$AF, plasma_corr$AF, method = 'spearman') #cor = 0.099 p = 0.82

pat_8_breast_2_plasma <- intersect(pat_8_breast_2$location, pat_8_plasma$location)
breast_2_corr <- pat_8_breast_2[pat_8_breast_2$location %in% pat_8_breast_2_plasma, ]
plasma_corr <- pat_8_plasma[pat_8_plasma$location %in% pat_8_breast_2_plasma, ]
cor.test(breast_2_corr$AF, plasma_corr$AF, method = 'spearman') #cor = 0.81 p = 1.50e-4
