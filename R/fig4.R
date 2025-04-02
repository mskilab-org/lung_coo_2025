library(skitools)
library(ggalluvial)
library(wesanderson)

# fig 4A

pts_coo_id <- readRDS('../data/pts_coo_id.rds')
pts_coo_id$Origin = factor(pts_coo_id$Origin, levels = c(  "Distal Lung_TP53 MUT", "Ambiguous_TP53 MUT", "Proximal Lung_TP53 MUT",  "Distal Lung_WT", "Ambiguous_WT", "Proximal Lung_WT"))
pts_coo_id$Identity = factor(pts_coo_id$Identity, levels = c("Distal Lung", "Non-Distal Lung"))

ggplot(pts_coo_id,
       aes(y = Patient_Frequency,
           axis1 = Origin,
           axis2 = Identity)) +
  theme_bw() +
  geom_alluvium(aes(fill = Origin), width = 1/12) +
  geom_stratum(width = 1/12, discern = TRUE, fill = "grey90", color = "black") +
  geom_label(stat = "stratum", size = 3.5, aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Origin", "Identity"), expand = c(0.05, 0.05)) +
  scale_fill_manual(values = c(wes_palettes$GrandBudapest2[4], wes_palettes$Chevalier1[1], wes_palettes$GrandBudapest2[3], wes_palettes$Chevalier1[2], wes_palettes$GrandBudapest2[2], wes_palettes$Chevalier1[3], wes_palettes$AsteroidCity3[1], wes_palettes$AsteroidCity1[1])) +
  ggtitle("Fig 4A - Origin --> Identity") + coord_flip()



# Fig 4B

pts_coo_id_4B <- readRDS( '../data/pts_coo_id_4B.rds')

res.plot =  pts_coo_id_4B[, prop.test(sum(Identity == 'Distal Lung'), .N) %>% dflm %>% cbind(nprox = sum(Identity == 'Distal Lung'), tot = .N), by = .(TP53_mut = ifelse(TP53_mut, 'TP53 MUT', 'WT'), Origin)][, fracprox := estimate]
res.plot$Origin = factor(res.plot$Origin, levels = c('Non-Distal Lung','Distal Lung'))
res.plot$TP53_mut = factor(res.plot$TP53_mut, levels = c('TP53 MUT','WT'))
ggplot(res.plot, aes(x = TP53_mut, y = fracprox, fill = TP53_mut)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.15) + theme_bw() + 
  facet_grid(~Origin) +
  scale_fill_manual(values = c(wes_palettes$Royal1[2], wes_palettes$Royal1[1])) +
  labs(title = '', x = '', y = '') + theme_classic() +
  theme(plot.title = element_text(size = 0, face = 'bold'),
        axis.text.x = element_text(size = 0, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 22, angle = 0, hjust = 1),
        axis.title.x = element_text(size = 0, face = 'plain'),
        axis.title.y = element_text(size = 5, face = 'bold'),
        axis.ticks.x = element_blank()) + 
  geom_text(mapping = aes(x = TP53_mut, y = ci.upper + 0.05, label = paste0(nprox, '/', tot)), size = 7) +
  guides(fill = guide_legend(title = 'Fig 4B - Distal fraction')) + theme(legend.position = "bottom")

