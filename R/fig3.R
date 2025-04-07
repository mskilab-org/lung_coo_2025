library(skitools)
library(data.table)
library(skitools)
library(skidb)
library(plyr)
library(dplyr)
library(MASS)
library(wesanderson)
library(readxl)
library(ggalluvial)
library(wesanderson)
library(readxl)
library(ggsankey)
library(effects)
library(stats)
library(forcats) 
library(ggpubr)
library(ggforce)

mat.in.2 <- readRDS('../data/mat.in.2.rds')
col_fun2 <- readRDS('../data/col_fun2.rds')
# ------------------------------------------------------------------------------------------------
# figure 3A
# ------------------------------------------------------------------------------------------------
column_ha3_1 = HeatmapAnnotation(
  Tobacco = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,29], ylim = c(0, 200000), gp = gpar(col = "red", fill = "#FF0000") , axis = FALSE),
  ID3 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,32], ylim = c(0, 10000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  
  SBS1 = anno_barplot( mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,25], ylim = c(0, 3000), gp = gpar(col = "darkgreen", fill = "darkgreen") , axis = FALSE) ,
  SBS5 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,26], ylim = c(0, 60000), gp = gpar(col = "darkgreen", fill = "darkgreen") , axis = FALSE) ,
  
  Apobec = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,27], ylim = c(0, 40000), gp = gpar(col = "orange", fill = "#FF8000") , axis = FALSE) ,
  
  ID1 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,30], ylim = c(0, 4000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  ID2 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,31], ylim = c(0, 7000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  ID12 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,33], ylim = c(0, 3200), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  
  purity = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,38], ylim = c(0, 1), gp = gpar(col = "green", fill = "green") , axis = FALSE),
  TMB =  anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,37], ylim = c(0, 100), gp = gpar(col = "brown", fill = "brown") , axis = FALSE),
  Smoking_Status = mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,34],
  
  col = list(
    Smoking_Status = c("Smoker" = "orange", "Never Smoker" = "blue") #,
  ),
  show_annotation_name = FALSE
)


column_ha3_2 = HeatmapAnnotation(
  Tobacco = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,29], ylim = c(0, 200000), gp = gpar(col = "red", fill = "#FF0000") ),
  ID3 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,32], ylim = c(0, 10000), gp = gpar(col = "#9900CC", fill = "#9900CC") ) ,
  
  SBS1 = anno_barplot( mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,25], ylim = c(0, 3000), gp = gpar(col = "darkgreen", fill = "darkgreen") ) ,
  SBS5 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,26], ylim = c(0, 60000), gp = gpar(col = "darkgreen", fill = "darkgreen") ) ,
  
  Apobec = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,27], ylim = c(0, 40000), gp = gpar(col = "orange", fill = "#FF8000") ) ,
  
  ID1 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,30], ylim = c(0, 4000), gp = gpar(col = "#9900CC", fill = "#9900CC") ) ,
  ID2 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,31], ylim = c(0, 7000), gp = gpar(col = "#9900CC", fill = "#9900CC") ) ,
  ID12 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,33], ylim = c(0, 3200), gp = gpar(col = "#9900CC", fill = "#9900CC") ) ,
  
  purity = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,38], ylim = c(0, 1), gp = gpar(col = "green", fill = "green") ),
  TMB =  anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,37], ylim = c(0, 100), gp = gpar(col = "brown", fill = "brown") ),
  Smoking_Status = mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,34],
  
  col = list(
    Smoking_Status = c("Smoker" = "orange", "Never Smoker" = "blue") #,
  ),
  show_annotation_name = FALSE
)


column_ha3_3 = HeatmapAnnotation(
  Tobacco = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,29], ylim = c(0, 200000), gp = gpar(col = "red", fill = "#FF0000") , axis = FALSE),
  ID3 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,32], ylim = c(0, 10000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  
  SBS1 = anno_barplot( mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,25], ylim = c(0, 3000), gp = gpar(col = "darkgreen", fill = "darkgreen") , axis = FALSE) ,
  SBS5 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,26], ylim = c(0, 60000), gp = gpar(col = "darkgreen", fill = "darkgreen") , axis = FALSE) ,
  
  Apobec = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,27], ylim = c(0, 40000), gp = gpar(col = "orange", fill = "#FF8000") , axis = FALSE) ,
  
  ID1 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,30], ylim = c(0, 4000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  ID2 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,31], ylim = c(0, 7000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  ID12 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,33], ylim = c(0, 3200), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  
  purity = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,38], ylim = c(0, 1), gp = gpar(col = "green", fill = "green") , axis = FALSE),
  TMB =  anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,37], ylim = c(0, 100), gp = gpar(col = "brown", fill = "brown") , axis = FALSE),
  Smoking_Status = mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,34],
  
  col = list(
    Smoking_Status = c("Smoker" = "orange", "Never Smoker" = "blue") #,
  ),
  show_annotation_name = FALSE
)



column_ha3_4 = HeatmapAnnotation(
  #snv = anno_barplot(mat.in.2[,26], gp = gpar(col = "black", fill = "#4575B4")),
  Tobacco = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,29], ylim = c(0, 200000), gp = gpar(col = "red", fill = "#FF0000") , axis = FALSE),
  ID3 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,32], ylim = c(0, 10000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  
  SBS1 = anno_barplot( mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,25], ylim = c(0, 3000), gp = gpar(col = "darkgreen", fill = "darkgreen") , axis = FALSE) ,
  SBS5 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,26], ylim = c(0, 60000), gp = gpar(col = "darkgreen", fill = "darkgreen") , axis = FALSE) ,
  
  Apobec = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,27], ylim = c(0, 40000), gp = gpar(col = "orange", fill = "#FF8000") , axis = FALSE) ,
  
  ID1 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,30], ylim = c(0, 4000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  ID2 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,31], ylim = c(0, 7000), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  ID12 = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,33], ylim = c(0, 3200), gp = gpar(col = "#9900CC", fill = "#9900CC") , axis = FALSE) ,
  
  purity = anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,38], ylim = c(0, 1), gp = gpar(col = "green", fill = "green") , axis = FALSE),
  TMB =  anno_barplot(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,37], ylim = c(0, 100), gp = gpar(col = "brown", fill = "brown") , axis = FALSE),
  Smoking_Status = mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,34],
  
  col = list(
    Smoking_Status = c("Smoker" = "orange", "Never Smoker" = "blue") #,
  ),
  show_annotation_name = TRUE
)


set.seed(90210)
ht1 <- Heatmap(t(mat.in.2[mat.in.2$cluster_sk_non_apobec == 1,2:24]), name = "HT1 Relative Risk", col = col_fun2, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 10), 
               column_names_side = c("bottom"), show_column_names = FALSE, column_km = 1, column_km_repeats = 100, top_annotation = column_ha3_1,
               show_parent_dend_line = FALSE, column_gap = unit(c(4), "mm"), column_title = NULL, border = TRUE)

set.seed(90210)
ht2 <- Heatmap(t(mat.in.2[mat.in.2$cluster_sk_non_apobec == 2,2:24]), name = "HT2 Relative Risk", col = col_fun2, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 10), 
               column_names_side = c("bottom"), show_column_names = FALSE, column_km = 1, column_km_repeats = 100, top_annotation = column_ha3_2,
               show_parent_dend_line = FALSE, column_gap = unit(c(4), "mm"), column_title = NULL, border = TRUE)

set.seed(90210)
ht3 <- Heatmap(t(mat.in.2[mat.in.2$cluster_sk_non_apobec == 3,2:24]), name = "HT3 Relative Risk", col = col_fun2, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 10), 
               column_names_side = c("bottom"), show_column_names = FALSE, column_km = 1, column_km_repeats = 100, top_annotation = column_ha3_3,
               show_parent_dend_line = FALSE, column_gap = unit(c(4), "mm"), column_title = NULL, border = TRUE)

set.seed(90210)
ht4 <- Heatmap(t(mat.in.2[mat.in.2$cluster_sk_non_apobec == 4,2:24]), name = "HT4 Relative Risk", col = col_fun2, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 10), 
               column_names_side = c("bottom"), show_column_names = FALSE, column_km = 1, column_km_repeats = 100, top_annotation = column_ha3_4,
               show_parent_dend_line = FALSE, column_gap = unit(c(4), "mm"), column_title = NULL, border = TRUE)


ht_list =  ht2+ ht1 + ht3 + ht4 
draw(ht_list, row_title = "Heatmap list", column_title = "Heatmap list")


# ------------------------------------------------------------------------------------------------
# figure 3B
# ------------------------------------------------------------------------------------------------
# Distal Lung

fig3B_data <- readRDS('../data/fig3B_data.rds')
ggplot(fig3B_data[cluster == 'Distal'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0(in.title, 'LUAD - Cluster 1')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# Ambiguous
ggplot(fig3B_data[cluster == 'Ambiguous'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0(in.title, 'LUAD - Cluster 2')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# Proximal_1

ggplot(fig3B_data[cluster == 'Proximal_1'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0(in.title, 'LUAD - Cluster 3')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# Proximal_2

ggplot(fig3B_data[cluster == 'Proximal_2'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0(in.title, 'LUAD - Cluster 4')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# ------------------------------------------------------------------------------------------------
# figure 3C
# ------------------------------------------------------------------------------------------------

# tmb

ik = 'tmb_75_pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")


# Tobacco  
extd_data_3 <- fread('../data/extd_data_3.csv')
ik = 'Tobacco_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))
prox.mut = res.plot[clust.int == 'Proximal',]$nmut
prox.wt = res.plot[clust.int == 'Proximal',]$tot - res.plot[clust.int == 'Proximal',]$nmut
dist.mut = res.plot[clust.int == 'Distal',]$nmut
dist.wt = res.plot[clust.int == 'Distal',]$tot - res.plot[clust.int == 'Distal',]$nmut
prox.dist.fisher = matrix(c(prox.mut, dist.mut, prox.wt, dist.wt), nrow = 2, byrow = TRUE) %>% fisher.test %>% dflm

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  # labs(title = paste0(ik, ' - p = ', round(prox.dist.fisher$p,6), ', e = ', round(prox.dist.fisher$estimate,3)), x = '', y = '') + theme_classic() +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# ID3

ik = 'ID3_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# Apobec

ik = 'Apobec_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# SBS1

ik = 'SBS1_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# SBS5

ik = 'SBS5_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# ID1

ik = 'ID1_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# ID12

ik = 'ID12_75pct'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# ------------------------------------------------------------------------------------------------
# figure 3D
# ------------------------------------------------------------------------------------------------

ik = 'Never_Smoker'
freq.plot.tmp.f = extd_data_3[!is.na(CellOfOrigin)]
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# ------------------------------------------------------------------------------------------------
# figure 3E
# ------------------------------------------------------------------------------------------------

# KRAS

ik = 'KRAS'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# TP53

ik = 'TP53'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# STK11

ik = 'STK11'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# EGFR


ik = 'EGFR'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")

# SMARCA4

ik = 'SMARCA4'
freq.plot.tmp.f = extd_data_3 
freq.plot.tmp.f$clust.int = freq.plot.tmp.f$CellOfOrigin
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Proximal','Distal', 'Ambiguous'))

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'olivedrab3', 'darkgoldenrod3', 'darkgrey')) +
  labs(title = paste0(ik), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = '')) + theme(legend.position = "top")


