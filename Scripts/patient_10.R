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

## PATIENT 10 ----
pat_10_liver_1 <- read.delim('Data/Patient_10/pat_10_liver_1_all_int_clean_hg19_ann.txt', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_10_liver_1 <- mutect_process(pat_10_liver_1) #40

pat_10_liver_2a <- read.delim('Data/Patient_10/pat_10_liver_2a_all_int_clean_hg19_ann.txt', header = TRUE, 
                              stringsAsFactors = FALSE, sep = '\t')
pat_10_liver_2a <- mutect_process(pat_10_liver_2a) #39

pat_10_liver_5 <- read.delim('Data/Patient_10/pat_10_liver_5_all_int_clean_hg19_ann.txt', header = TRUE, 
                             stringsAsFactors = FALSE, sep = '\t')
pat_10_liver_5 <- mutect_process(pat_10_liver_5) #50

pat_10_plasma <- read.delim('Data/Patient_10/pat_10_plasma_all_int_clean_hg19_ann.txt', header = TRUE, 
                            stringsAsFactors = FALSE, sep = '\t')
pat_10_plasma <- mutect_process(pat_10_plasma) #4299

## looking at how well plasma detects tumor mutations ----

length(intersect(pat_10_liver_1$location, pat_10_plasma$location)) #14/40
length(intersect(pat_10_liver_2a$location, pat_10_plasma$location)) #11/39
length(intersect(pat_10_liver_5$location, pat_10_plasma$location)) #13/50


#pool mutations from 3 mets
pat_10_met_pool <- unique(c(pat_10_liver_1$location, pat_10_liver_2a$location, pat_10_liver_5$location)) #73 unique mutations


pat_10_plasma_found <- pat_10_plasma[pat_10_plasma$location %in% pat_10_met_pool, ] #16 mutations
pat_10_plasma_found_vars <- (pat_10_plasma_found$location) 


pat_10_liver_1_pooled <- pat_10_liver_1[pat_10_liver_1$location %in% pat_10_plasma_found_vars, ]

pat_10_liver_2a_pooled <- pat_10_liver_2a[pat_10_liver_2a$location %in% pat_10_plasma_found_vars, ]

pat_10_liver_5_pooled <- pat_10_liver_5[pat_10_liver_5$location %in% pat_10_plasma_found_vars, ]

pat_10_plasma_pooled <- pat_10_plasma[pat_10_plasma$location %in% pat_10_plasma_found_vars, ]

# subset these to tumor AF for new figure
pat_10_liver_1_pooled <- pat_10_liver_1_pooled[, c('location', 'AF')]
colnames(pat_10_liver_1_pooled) <- c('location', 'AF_liver_1')

pat_10_liver_2a_pooled <- pat_10_liver_2a_pooled[, c('location', 'AF')]
colnames(pat_10_liver_2a_pooled) <- c('location', 'AF_liver_2a')

pat_10_liver_5_pooled <- pat_10_liver_5_pooled[, c('location', 'AF')]
colnames(pat_10_liver_5_pooled) <- c('location', 'AF_liver_5')

pat_10_plasma_pooled <- pat_10_plasma_pooled[, c('location', 'AF')]
colnames(pat_10_plasma_pooled) <- c('location', 'AF_plasma')

# put them together
pat_10_muts_pooled <- merge(pat_10_liver_1_pooled, pat_10_liver_2a_pooled, by = 'location', all = TRUE)
pat_10_muts_pooled <- merge(pat_10_muts_pooled, pat_10_liver_5_pooled, by = 'location', all = TRUE)
pat_10_muts_pooled <- merge(pat_10_muts_pooled, pat_10_plasma_pooled, by = 'location', all = TRUE)
pat_10_muts_pooled <- pat_10_muts_pooled[order(pat_10_muts_pooled$AF_plasma, decreasing = TRUE), ]

rownames(pat_10_muts_pooled) <- pat_10_muts_pooled$location
row_muts <- rownames(pat_10_muts_pooled)
pat_10_muts_pooled <- pat_10_muts_pooled[, -1]

pat_10_muts_pooled <- as.data.frame(pat_10_muts_pooled)
rownames(pat_10_muts_pooled) <- row_muts

#plot 

pat_10_3 <- c(pat_10_muts_pooled$AF_liver_1, pat_10_muts_pooled$AF_liver_2a, pat_10_muts_pooled$AF_liver_5)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_10_muts_pooled)[1] * dim(pat_10_muts_pooled)[2])
# plasma reds
cell_cols[49:64] <- color.scale(pat_10_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:48] <- color.scale(pat_10_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 16, byrow = FALSE)
pat_10_pooled_t <- data.frame(t(pat_10_muts_pooled))
pat_10_pooled_t <- pat_10_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]

# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_10_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

#add legends
legval<-seq(min(pat_10_muts_pooled[, 4], na.rm = TRUE),max(pat_10_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(1,-0.9,11,-0.5,round(c(min(pat_10_muts_pooled[, 4], na.rm = TRUE), max(pat_10_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=2.2, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE),max(pat_10_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(12,-0.9,23,-0.5,round(c(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE), max(pat_10_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=13, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 11.5, cex = 1.1, font = 2)

# add NA legend
color.legend(27, -0.9, 28, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.8, at=24.9, cex = 1.1, font = 2)
legend(x=27.1,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_10_muts_pooled)
mut_col_labels[1:10] <- c('', '', '', '', '', '', 'EP300 p.Q2223P', 'SLX4 p.N457K', 'FGFR4 p.G388R', '')
mut_col_labels[11:16] <- c('SLX4 p.A952V', 'SLX4 p.A952T', '', '', '', '')
axis(3, at = (1:ncol(pat_10_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.8, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver 1', 'Liver 2a', 'Liver 5')
axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_10_muts_pooled$AF_liver_1)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_10_muts_pooled$AF_liver_1))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_10_muts_pooled$AF_liver_2a)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_10_muts_pooled$AF_liver_2a))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_10_muts_pooled$AF_liver_5)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_10_muts_pooled$AF_liver_5))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))