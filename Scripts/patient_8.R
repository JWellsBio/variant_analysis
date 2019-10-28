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
pat_8_axillary <- read.delim('Data/Patient_8/pat_8_axillary_fresh_mutect_filt_hg19_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_8_axillary <- mutect_process(pat_8_axillary) #23

pat_8_breast_1 <- read.delim('Data/Patient_8/pat_8_breast_1_fresh_mutect_filt_hg19_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_8_breast_1 <- mutect_process(pat_8_breast_1) #35

pat_8_breast_2 <- read.delim('Data/Patient_8/pat_8_breast_2_fresh_mutect_filt_hg19_ann.tsv', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t') #512
pat_8_breast_2 <- mutect_process(pat_8_breast_2)

pat_8_plasma <- read.delim('Data/Patient_8/pat_8_plasma_fresh_mutect_filt_hg19_ann.tsv', header = TRUE, 
                           stringsAsFactors = FALSE, sep = '\t')
pat_8_plasma <- mutect_process(pat_8_plasma, sample_type = 'plasma') #141


## looking at how well plasma detects tumor mutations ----
length(intersect(pat_8_axillary$location, pat_8_plasma$location)) #7/23
length(intersect(pat_8_breast_1$location, pat_8_plasma$location)) #8/35
length(intersect(pat_8_breast_2$location, pat_8_plasma$location)) #42/512


#pool mutations from 3 mets
pat_8_met_pool <- unique(c(pat_8_axillary$location, pat_8_breast_1$location, pat_8_breast_2$location)) #525 unique mutations


pat_8_plasma_found <- pat_8_plasma[pat_8_plasma$location %in% pat_8_met_pool, ] #44 mutations
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
cell_cols[133:176] <- color.scale(pat_8_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:132] <- color.scale(pat_8_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 44, byrow = FALSE)
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
