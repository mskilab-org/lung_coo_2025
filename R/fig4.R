## Figure 4A ##
library(skitools)

extd_data_3 <- fread('../data/Fig4/4ABC_extd_data_3.csv')
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

## Figure 4B ##
#Load libraries.
library(skitools)

# tmb
extd_data_3 <- fread('../data/Fig4/4ABC_extd_data_3.csv')
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

## Figure 4C ##
#Load libraries.
library(skitools)

# KRAS
extd_data_3 <- fread('../data/Fig4/4ABC_extd_data_3.csv')
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
