## setting data for 4-dimension figure showing what is found
## data calls on other scripts already being called

variants_found <- c(nrow(pat_9_ln), nrow(pat_9_oment), nrow(pat_9_ovary), nrow(pat_8_axillary), nrow(pat_8_breast_1), nrow(pat_8_breast_2), 
                    nrow(pat_2_liver_1), nrow(pat_2_liver_2), nrow(pat_2_breast_1), nrow(pat_10_liver_1), nrow(pat_10_liver_2a), 
                    nrow(pat_10_liver_5), nrow(pat_ema_heart), nrow(pat_ema_l_kidney), nrow(pat_ema_r_kidney), nrow(pat_ema_liver_1), 
                    nrow(pat_ema_liver_2), nrow(pat_ema_oment_1), nrow(pat_ema_oment_2))
found_in_plasma <- c(nrow(pat_9_ln_pooled), nrow(pat_9_oment_pooled), nrow(pat_9_ovary_pooled), nrow(pat_8_axillary_pooled), nrow(pat_8_breast_1_pooled), 
                     nrow(pat_8_breast_2_pooled), nrow(pat_2_liver_1_pooled), nrow(pat_2_liver_2_pooled), nrow(pat_2_breast_1_pooled), 
                     nrow(pat_10_liver_1_pooled), nrow(pat_10_liver_2a_pooled), nrow(pat_10_liver_5_pooled), nrow(pat_ema_heart_pooled), 
                     nrow(pat_ema_l_kidney_pooled), nrow(pat_ema_r_kidney_pooled), nrow(pat_ema_liver_1_pooled), 
                     nrow(pat_ema_liver_2_pooled), nrow(pat_ema_oment_1_pooled), nrow(pat_ema_oment_2_pooled))
sample_name <- c(rep('Patient 9', 3), rep('Patient 8', 3), rep('Patient 2', 3), rep('Patient 10', 3), rep('Patient EMA', 7))
mean_MAF <- c(mean(pat_9_ln_pooled$AF), mean(pat_9_oment_pooled$AF), mean(pat_9_ovary_pooled$AF), mean(pat_8_axillary_pooled$AF), mean(pat_8_breast_1_pooled$AF), 
              mean(pat_8_breast_2_pooled$AF), mean(pat_2_liver_1_pooled$AF), mean(pat_2_liver_2_pooled$AF), mean(pat_2_breast_1_pooled$AF), 
              mean(pat_10_liver_1_pooled$AF), mean(pat_10_liver_2a_pooled$AF), mean(pat_10_liver_5_pooled$AF), mean(pat_ema_heart_pooled$AF), 
              mean(pat_ema_l_kidney_pooled$AF), mean(pat_ema_r_kidney_pooled$AF), mean(pat_ema_liver_1_pooled$AF), 
              mean(pat_ema_liver_2_pooled$AF), mean(pat_ema_oment_1_pooled$AF), mean(pat_ema_oment_2_pooled$AF))
graph_df <- data.frame(variants_found, found_in_plasma, sample_name, mean_MAF)
graph_df$percent_found <- graph_df$found_in_plasma/graph_df$variants_found


graph_df <- graph_df[graph_df$variants_found < 300, ]

library(ggplot2)
ggplot(graph_df, aes(x=variants_found, y=percent_found)) + geom_point(aes(size=graph_df$mean_MAF, shape=as.factor(graph_df$sample_name))) + theme_bw()
