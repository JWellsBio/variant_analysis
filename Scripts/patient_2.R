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
pat_2_liver_1 <- read.delim('Data/Patient_2/pat_2_liver_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                            stringsAsFactors = FALSE, sep = '\t')
pat_2_liver_1 <- mutect_process(pat_2_liver_1) #20

pat_2_liver_2 <- read.delim('Data/Patient_2/pat_2_liver_2_all_int_clean_hg19_ann.txt', header = TRUE, 
                            stringsAsFactors = FALSE, sep = '\t')
pat_2_liver_2 <- mutect_process(pat_2_liver_2) #28

pat_2_breast_1 <- read.delim('Data/Patient_2/pat_2_breast_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_2_breast_1 <- mutect_process(pat_2_breast_1) #34

pat_2_breast_2 <- read.delim('Data/Patient_2/pat_2_breast_2_all_int_clean_hg19_ann.txt', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_2_breast_2 <- mutect_process(pat_2_breast_2) #485

pat_2_plasma <- read.delim('Data/Patient_2/pat_2_plasma_all_int_clean_hg19_ann.txt', header = TRUE, 
                           stringsAsFactors = FALSE, sep = '\t')
pat_2_plasma <- mutect_process(pat_2_plasma, sample_type = 'plasma') #141


## looking at how well plasma detects tumor mutations ----

length(intersect(pat_2_breast_1$location, pat_2_plasma$location)) #9/34
length(intersect(pat_2_breast_2$location, pat_2_plasma$location)) #42/485
length(intersect(pat_2_liver_1$location, pat_2_plasma$location)) #2/20
length(intersect(pat_2_liver_2$location, pat_2_plasma$location)) #7/28



#pool mutations from 3 mets
pat_2_met_pool <- unique(c(pat_2_liver_1$location, pat_2_liver_2$location, pat_2_breast_1$location)) #60 unique mutations


pat_2_plasma_found <- pat_2_plasma[pat_2_plasma$location %in% pat_2_met_pool, ] #15 mutations
pat_2_plasma_found_vars <- (pat_2_plasma_found$location)


# subset to pooled plasma and take a look
pat_2_liver_1_pooled <- pat_2_liver_1[pat_2_liver_1$location %in% pat_2_plasma_found_vars, ]

pat_2_liver_2_pooled <- pat_2_liver_2[pat_2_liver_2$location %in% pat_2_plasma_found_vars, ]

pat_2_breast_1_pooled <- pat_2_breast_1[pat_2_breast_1$location %in% pat_2_plasma_found_vars, ]

#pat_2_breast_2_pooled <- pat_2_breast_2[pat_2_breast_2$location %in% pat_2_plasma_found_vars, ]

pat_2_plasma_pooled <- pat_2_plasma[pat_2_plasma$location %in% pat_2_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_2_liver_1_pooled <- pat_2_liver_1_pooled[, c('location', 'AF')]
colnames(pat_2_liver_1_pooled) <- c('location', 'AF_liver_1')

pat_2_liver_2_pooled <- pat_2_liver_2_pooled[, c('location', 'AF')]
colnames(pat_2_liver_2_pooled) <- c('location', 'AF_liver_2')

pat_2_breast_1_pooled <- pat_2_breast_1_pooled[, c('location', 'AF')]
colnames(pat_2_breast_1_pooled) <- c('location', 'AF_breast_1')

#pat_2_breast_2_pooled <- pat_2_breast_2_pooled[, c('location', 'AF')]
#colnames(pat_2_breast_2_pooled) <- c('location', 'AF_breast_2')
#pat_2_breast_2_pooled$AF_breast_2[4] <- 0.242
#pat_2_breast_2_pooled$AF_breast_2 <- as.numeric(pat_2_breast_2_pooled$AF_breast_2)

pat_2_plasma_pooled <- pat_2_plasma_pooled[, c('location', 'AF')]
colnames(pat_2_plasma_pooled) <- c('location', 'AF_plasma')

# put them together
pat_2_muts_pooled <- merge(pat_2_liver_1_pooled, pat_2_liver_2_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_breast_1_pooled, by = 'location', all = TRUE)
#pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_breast_2_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- merge(pat_2_muts_pooled, pat_2_plasma_pooled, by = 'location', all = TRUE)
pat_2_muts_pooled <- pat_2_muts_pooled[order(pat_2_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_2_muts_pooled) <- pat_2_muts_pooled$location
row_muts <- rownames(pat_2_muts_pooled)
pat_2_muts_pooled <- pat_2_muts_pooled[, -1]

pat_2_muts_pooled <- as.data.frame(pat_2_muts_pooled)
rownames(pat_2_muts_pooled) <- row_muts

#plot 

pat_2_3 <- c(pat_2_muts_pooled$AF_liver_1, pat_2_muts_pooled$AF_liver_2, pat_2_muts_pooled$AF_breast_1)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_2_muts_pooled)[1] * dim(pat_2_muts_pooled)[2])
# plasma reds
cell_cols[46:60] <- color.scale(pat_2_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:45] <- color.scale(pat_2_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 15, byrow = FALSE)
pat_2_pooled_t <- data.frame(t(pat_2_muts_pooled))
pat_2_pooled_t <- pat_2_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]

# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
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
color.legend(6.5, -0.9, 6.75, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.2, at=5.9, cex = 1.1, font = 2)
legend(x=6.525,y=-0.42,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_2_muts_pooled)
mut_col_labels[1:5] <- c('', '', '', '', '')
axis(3, at = (1:ncol(pat_2_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver 1', 'Liver 2', 'Breast 1')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_2_muts_pooled$AF_liver_1)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_2_muts_pooled$AF_liver_1))), 
       pch = 16)

points(x = which(is.na(pat_2_muts_pooled$AF_liver_2)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_2_muts_pooled$AF_liver_2))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_2_muts_pooled$AF_breast_1)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_2_muts_pooled$AF_breast_1))), 
       pch = 16)




par(mar=c(5.1,4.1,4.1,2.1))
