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






# ------------------------------------------------------------------------------------------------
# EDF 6A
# ------------------------------------------------------------------------------------------------
sims_acc_results <- fread('../data/edf6_luad_sims_results.csv')

p1_4 <- ggplot(sims_acc_results[, .(celltype, accuracy_level_finest)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level_finest)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,1) + 
  ggtitle("Finest level match") +
  xlab("True cell type") + ylab("Accuracy")

p1_3 <- ggplot(sims_acc_results[, .(celltype, accuracy_level4)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level4)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Level 4 match") +
  xlab("True cell type") + ylab("Accuracy")

p1_2 <- ggplot(sims_acc_results[, .(celltype, accuracy_level3)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level3)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Level 3 match") +
  xlab("True cell type") + ylab("Accuracy")

p1_1 <- ggplot(sims_acc_results[, .(celltype, accuracy_level2)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level2)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Level 2 match") +
  xlab("True cell type") + ylab("Accuracy")

p1_0 <- ggplot(sims_acc_results[, .(celltype, accuracy_leveldnd)] %>% as.data.frame(), aes(x = celltype, y = accuracy_leveldnd)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Distal vs non-distal match") +
  xlab("True cell type") + ylab("Accuracy")

p_acc <- ggpubr::ggarrange(p1_0, p1_1, p1_2, p1_3, p1_4, ncol = 3, nrow = 2)
print(p_acc)








# ------------------------------------------------------------------------------------------------
# EDF 8A
# ------------------------------------------------------------------------------------------------

col_fun2 <- readRDS('../data/col_fun2.rds')
mat.in.2 <- fread('../data/mat.in.2.lusc.csv')
mat.in.2 <- mat.in.2 %>% as.data.frame()
column_ha2 = HeatmapAnnotation(
  sbs4 = anno_barplot(mat.in.2[,25], gp = gpar(col = "red", fill = "#FF0000")),
  sbs1 = anno_barplot( mat.in.2[,26], gp = gpar(col = "green", fill = "#99FF33")) ,
  sbs5 = anno_barplot(mat.in.2[,27], gp = gpar(col = "green", fill = "#99FF33")) ,
  sbs2 = anno_barplot(mat.in.2[,28], gp = gpar(col = "orange", fill = "#FF8000")) ,
  sbs13 = anno_barplot(mat.in.2[,29], gp = gpar(col = "orange", fill = "#FF8000")) ,
  id1 = anno_barplot(mat.in.2[,31], gp = gpar(col = "green", fill = "#99FF33")) ,
  id2 = anno_barplot(mat.in.2[,32], gp = gpar(col = "violet", fill = "#9900CC")) ,
  id3 = anno_barplot(mat.in.2[,33], gp = gpar(col = "red", fill = "#FF0000")) ,
  id12 = anno_barplot(mat.in.2[,34], gp = gpar(col = "violet", fill = "#9900CC")) ,
  smoker = mat.in.2[,35],
  col = list(
    smoker = c("Smoker" = "orange", "Never Smoker" = "blue") 
  )
)

set.seed(90210)
Heatmap(t(mat.in.2[,1:23]), name = "Relative Risk", col = col_fun2, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 10), 
            column_names_side = c("bottom"), show_column_names = FALSE, column_km = 4, column_km_repeats = 100, top_annotation = column_ha2,
            show_parent_dend_line = FALSE, column_gap = unit(c(4), "mm"), column_title = NULL, border = TRUE)




# ------------------------------------------------------------------------------------------------
# EDF 9C
# ------------------------------------------------------------------------------------------------

edf9 <- readRDS('../data/edf9.rds')

ik <-  'Papillary'


freq.plot.tmp.f = edf9

freq.plot.tmp.f$clust.int = freq.plot.tmp.f$T_Identity
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Non-distal identity','Distal identity'))
non.dist.mut = res.plot[clust.int == 'Non-distal identity',]$nmut
non.dist.wt = res.plot[clust.int == 'Non-distal identity',]$tot - res.plot[clust.int == 'Non-distal identity',]$nmut
dist.mut = res.plot[clust.int == 'Distal identity',]$nmut
dist.wt = res.plot[clust.int == 'Distal identity',]$tot - res.plot[clust.int == 'Distal identity',]$nmut

prox.dist.fisher = matrix(c(non.dist.mut, dist.mut, non.dist.wt, dist.wt), nrow = 2, byrow = TRUE) %>% fisher.test %>% dflm


ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'forestgreen', 'darkgoldenrod3')) +
  labs(title = paste0('Papillary fraction of LUAD identity'), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = 'Identity')) + theme(legend.position = "top")






# ------------------------------------------------------------------------------------------------
# EDF 9D
# ------------------------------------------------------------------------------------------------
# by Identity 

edf9 <- readRDS('../data/edf9.rds')

ik <-  'NSCLC_NOS'

freq.plot.tmp.f = edf9

freq.plot.tmp.f$clust.int = freq.plot.tmp.f$T_Identity
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('Non-distal identity','Distal identity'))
non.dist.mut = res.plot[clust.int == 'Non-distal identity',]$nmut
non.dist.wt = res.plot[clust.int == 'Non-distal identity',]$tot - res.plot[clust.int == 'Non-distal identity',]$nmut
dist.mut = res.plot[clust.int == 'Distal identity',]$nmut
dist.wt = res.plot[clust.int == 'Distal identity',]$tot - res.plot[clust.int == 'Distal identity',]$nmut

prox.dist.fisher = matrix(c(non.dist.mut, dist.mut, non.dist.wt, dist.wt), nrow = 2, byrow = TRUE) %>% fisher.test %>% dflm


ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'forestgreen', 'darkgoldenrod3', '#bdbdbd')) +
  labs(title = paste0('NSCLC-NOS fraction of LUAD identity'), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = 'Identity')) + theme(legend.position = "top")





# ------------------------------------------------------------------------------------------------
# EDF 9E
# ------------------------------------------------------------------------------------------------
# NSCLC-NOS vs TP53

edf9 <- readRDS('../data/edf9.rds')

ik <-  'NSCLC_NOS'

freq.plot.tmp.f = edf9

freq.plot.tmp.f$clust.int = freq.plot.tmp.f$TP53_mut
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('TP53 MUT','WT'))
non.dist.mut = res.plot[clust.int == 'TP53 MUT',]$nmut
non.dist.wt = res.plot[clust.int == 'TP53 MUT',]$tot - res.plot[clust.int == 'TP53 MUT',]$nmut
dist.mut = res.plot[clust.int == 'WT',]$nmut
dist.wt = res.plot[clust.int == 'WT',]$tot - res.plot[clust.int == 'WT',]$nmut

prox.dist.fisher = matrix(c(non.dist.mut, dist.mut, non.dist.wt, dist.wt), nrow = 2, byrow = TRUE) %>% fisher.test %>% dflm


ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'salmon3', 'snow4')) +
  labs(title = paste0('NSCLC-NOS histology fraction'), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = 'TP53 status')) + theme(legend.position = "top")

# ------------------------------------------------------------------------------------------------
# EDF 9F
# ------------------------------------------------------------------------------------------------

edf9 <- readRDS('../data/edf9.rds')

ik = 'mut_tp53_hist'
freq.plot.tmp.f = edf9

freq.plot.tmp.f$clust.int = freq.plot.tmp.f$NSCLC_NOS_VS_Others
freq.plot.tmp.f = freq.plot.tmp.f[!duplicated(freq.plot.tmp.f),]
col.num = which(colnames(freq.plot.tmp.f) == ik)
freq.plot.tmp.f$mut.int = freq.plot.tmp.f[,..col.num]
freq.plot.tmp.f = freq.plot.tmp.f[!is.na(mut.int),]
freq.plot.tmp.f$mut.int = as.numeric(freq.plot.tmp.f$mut.int)
res.plot = freq.plot.tmp.f[, prop.test(sum(mut.int), .N) %>% dflm %>% cbind(nmut = sum(mut.int), tot = .N), by = .(clust.int)][, fracmut := estimate]
res.plot$clust.int = factor(res.plot$clust.int, levels = c('NSCLC_NOS','Non-NSCLC_NOS'))
nsclc.nos.mut = res.plot[clust.int == 'NSCLC_NOS',]$nmut
nsclc.nos.wt = res.plot[clust.int == 'NSCLC_NOS',]$tot - res.plot[clust.int == 'NSCLC_NOS',]$nmut
others.mut = res.plot[clust.int == 'Non-NSCLC_NOS',]$nmut
others.wt = res.plot[clust.int == 'Non-NSCLC_NOS',]$tot - res.plot[clust.int == 'Non-NSCLC_NOS',]$nmut

nsclc.nos.others.fisher = matrix(c(nsclc.nos.mut, others.mut, nsclc.nos.wt, others.wt), nrow = 2, byrow = TRUE) %>% fisher.test %>% dflm

ggplot(res.plot, aes(x = clust.int, y = fracmut, fill = clust.int)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  scale_fill_manual(values = c(  'orchid3',   'snow4')) +
  labs(title = paste0('Fraction of samples with TP53 mut'), x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 20, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = clust.int, y = ci.upper + 0.05, label = paste0(nmut, '/', tot)), size = 7) + ylim(0,1.05) +
  guides(fill = guide_legend(title = 'Histology')) + theme(legend.position = "top")


