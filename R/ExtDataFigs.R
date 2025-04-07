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
# EDF 8A
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






