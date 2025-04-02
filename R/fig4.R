
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
  ggtitle("") + coord_flip()



