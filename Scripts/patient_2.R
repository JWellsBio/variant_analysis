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

## Patient 2 ----
pat_2_liver_1 <- read.delim('Data/Patient_2/pat_2_liver_1_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                            stringsAsFactors = FALSE, sep = '\t')
pat_2_liver_1 <- mutect_process(pat_2_liver_1) #15

pat_2_liver_2 <- read.delim('Data/Patient_2/pat_2_liver_2_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                            stringsAsFactors = FALSE, sep = '\t')
pat_2_liver_2 <- mutect_process(pat_2_liver_2) #25

pat_2_breast_1 <- read.delim('Data/Patient_2/pat_2_breast_1_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_2_breast_1 <- mutect_process(pat_2_breast_1) #19

pat_2_breast_2 <- read.delim('Data/Patient_2/pat_2_breast_2_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_2_breast_2 <- mutect_process(pat_2_breast_2) #22

pat_2_plasma <- read.delim('Data/Patient_2/pat_2_plasma_exon_only_mutect_filtered_hg38_lift_ann.tsv', header = TRUE, 
                           stringsAsFactors = FALSE, sep = '\t')
pat_2_plasma <- mutect_process(pat_2_plasma, sample_type = 'plasma') #87


## looking at how well plasma detects tumor mutations ----
length(intersect(pat_2_breast_1$location, pat_2_plasma$location)) #6/19
length(intersect(pat_2_breast_2$location, pat_2_plasma$location)) #6/22
length(intersect(pat_2_liver_1$location, pat_2_plasma$location)) #4/15
length(intersect(pat_2_liver_2$location, pat_2_plasma$location)) #4/25



#pool mutations from 3 mets
#use breast 1 or breast 2
pat_2_met_pool <- unique(c(pat_2_liver_1$location, pat_2_liver_2$location, pat_2_breast_2$location)) #36 unique mutations w br 1 or 34 w br 2 


pat_2_plasma_found <- pat_2_plasma[pat_2_plasma$location %in% pat_2_met_pool, ] #6 mutations w br 1 or 6 w br 2
pat_2_plasma_found_vars <- (pat_2_plasma_found$location)


# subset to pooled plasma and take a look
pat_2_liver_1_pooled <- pat_2_liver_1[pat_2_liver_1$location %in% pat_2_plasma_found_vars, ]

pat_2_liver_2_pooled <- pat_2_liver_2[pat_2_liver_2$location %in% pat_2_plasma_found_vars, ]

#pat_2_breast_1_pooled <- pat_2_breast_1[pat_2_breast_1$location %in% pat_2_plasma_found_vars, ]

pat_2_breast_2_pooled <- pat_2_breast_2[pat_2_breast_2$location %in% pat_2_plasma_found_vars, ]

pat_2_plasma_pooled <- pat_2_plasma[pat_2_plasma$location %in% pat_2_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_2_liver_1_pooled <- pat_2_liver_1_pooled[, c('location', 'AF')]
colnames(pat_2_liver_1_pooled) <- c('location', 'AF_liver_1')

pat_2_liver_2_pooled <- pat_2_liver_2_pooled[, c('location', 'AF')]
colnames(pat_2_liver_2_pooled) <- c('location', 'AF_liver_2')

#pat_2_breast_1_pooled <- pat_2_breast_1_pooled[, c('location', 'AF')]
#colnames(pat_2_breast_1_pooled) <- c('location', 'AF_breast_1')

pat_2_breast_2_pooled <- pat_2_breast_2_pooled[, c('location', 'AF')]
colnames(pat_2_breast_2_pooled) <- c('location', 'AF_breast_2')

pat_2_plasma_pooled <- pat_2_plasma_pooled[, c('location', 'AF')]
colnames(pat_2_plasma_pooled) <- c('location', 'AF_plasma')

# put them together
pat_2_muts_pooled <- merge(pat_2_liver_1_pooled, pat_2_liver_2_pooled, by = 'location', all = TRUE)
#pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_breast_1_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_breast_2_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_plasma_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- pat_2_muts_pooled[order(pat_2_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_2_muts_pooled) <- pat_2_muts_pooled$location
row_muts <- rownames(pat_2_muts_pooled)
pat_2_muts_pooled <- pat_2_muts_pooled[, -1]

pat_2_muts_pooled <- as.data.frame(pat_2_muts_pooled)
rownames(pat_2_muts_pooled) <- row_muts

#plot 

pat_2_3 <- c(pat_2_muts_pooled$AF_liver_1, pat_2_muts_pooled$AF_liver_2, pat_2_muts_pooled$AF_breast_2)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_2_muts_pooled)[1] * dim(pat_2_muts_pooled)[2])
# plasma reds
cell_cols[19:24] <- color.scale(pat_2_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:18] <- color.scale(pat_2_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 6, byrow = FALSE)
pat_2_pooled_t <- data.frame(t(pat_2_muts_pooled))
pat_2_pooled_t <- pat_2_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]

# plot it
# extra space
par(mar=c(6,44.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_2_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_2_muts_pooled[, 4], na.rm = TRUE),max(pat_2_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(0.5,-0.9,4.5,-0.5,round(c(min(pat_2_muts_pooled[, 4], na.rm = TRUE), max(pat_2_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.3, at=1.0, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_2_muts_pooled[, 1:3], na.rm = TRUE),max(pat_2_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(5,-0.9,9,-0.5,round(c(min(pat_2_muts_pooled[, 1:3], na.rm = TRUE), max(pat_2_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.3, at=5.5, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 3.9, at=4.7, cex = 1.1, font = 2)

# add NA legend
color.legend(13.5, -0.9, 14, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.8, at=12.4, cex = 1.1, font = 2)
legend(x=13.55,y=-0.45,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_2_muts_pooled)
mut_col_labels[1:6] <- c('BARD1\nc.1519G>A', 'BARD1\nc.1518T>C', 'FLT1\nc.*1999G>A', 'FGFR2\nc.*313G>A', 'HNF1A\np.S581G', 'FGFR2\nc.*303G>A')
axis(3, at = (1:ncol(pat_2_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver 1', 'Liver 2', 'Breast 2')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_2_muts_pooled$AF_liver_1)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_2_muts_pooled$AF_liver_1))), 
       pch = 16)
#liver 2 points
points(x = which(is.na(pat_2_muts_pooled$AF_liver_2)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_2_muts_pooled$AF_liver_2))), 
       pch = 16)
# breast 2 points
points(x = which(is.na(pat_2_muts_pooled$AF_breast_2)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_2_muts_pooled$AF_breast_1))), 
       pch = 16)




par(mar=c(5.1,4.1,4.1,2.1))


##upset
pat_2_liver_1_met_muts <- pat_2_liver_1$location
pat_2_liver_1_met_muts <- as.data.frame(pat_2_liver_1_met_muts)
pat_2_liver_1_met_muts$mutation <- rep(1, nrow(pat_2_liver_1_met_muts))
colnames(pat_2_liver_1_met_muts) <- c('mutation', 'Liver_1_Met')

pat_2_liver_2_met_muts <- pat_2_liver_2$location
pat_2_liver_2_met_muts <- as.data.frame(pat_2_liver_2_met_muts)
pat_2_liver_2_met_muts$mutation <- rep(1, nrow(pat_2_liver_2_met_muts))
colnames(pat_2_liver_2_met_muts) <- c('mutation', 'Liver_2_Met')

pat_2_breast_2_met_muts <- pat_2_breast_2$location
pat_2_breast_2_met_muts <- as.data.frame(pat_2_breast_2_met_muts)
pat_2_breast_2_met_muts$mutation <- rep(1, nrow(pat_2_breast_2_met_muts))
colnames(pat_2_breast_2_met_muts) <- c('mutation', 'Breast_2_Met')

pat_2_plasma_muts <- pat_2_plasma$location
pat_2_plasma_muts <- as.data.frame(pat_2_plasma_muts)
pat_2_plasma_muts$mutation <- rep(1, nrow(pat_2_plasma_muts))
colnames(pat_2_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_2_liver_1_met_muts, pat_2_liver_2_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_2_breast_2_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_2_plasma_muts, by = 'mutation', all = TRUE)
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
                                                      colors = c(Plasma = 'gray60', Liver_1_Met = 'white', Liver_2_Met = 'white', 
                                                                 Breast_2_Met = 'white')))),
      
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray60', 'goldenrod4', 'aquamarine3', 'chocolate3'), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      mainbar.y.label = 'Number of Mutations\nin Common', 
      text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 2.0))

#red/black figures to look at how far down plasma can detect

pat_2_liver_1_met_stats <- pat_2_liver_1[, c('location', 'AF')]
pat_2_liver_2_met_stats <- pat_2_liver_2[, c('location', 'AF')]
pat_2_breast_2_met_stats <- pat_2_breast_2[, c('location', 'AF')]
pat_2_all <- rbind(pat_2_liver_1_met_stats, pat_2_liver_2_met_stats)
pat_2_all <- rbind(pat_2_all, pat_2_breast_2_met_stats)
pat_2_all$AF <- as.numeric(pat_2_all$AF)
pat_2_all <- pat_2_all[order(pat_2_all$AF, decreasing = TRUE), ]
pat_2_all$color <- ifelse(pat_2_all$location %in% pat_2_plasma$location, 'red', 'black')
barplot(pat_2_all$AF, col = pat_2_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 2\n(3 tumors)', ylim = c(0,1.0))

pat_2_plasma_met_stats <- pat_2_plasma[, c('location', 'AF')]
pat_2_plasma_met_stats$AF <- as.numeric(pat_2_plasma_met_stats$AF)
pat_2_plasma_met_stats <- pat_2_plasma_met_stats[order(pat_2_plasma_met_stats$AF, decreasing = TRUE), ]
pat_2_plasma_met_stats$color <- ifelse(pat_2_plasma_met_stats$location %in% pat_2_all$location, 'red', 'black')
barplot(pat_2_plasma_met_stats$AF, col = pat_2_plasma_met_stats$color, ylab = 'Mutant Allele Frequency in Plasma', 
        xlab = 'Variant', main = 'Patient 2\n(3 tumors)', ylim = c(0,1.0))
