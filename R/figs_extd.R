# ------------------------------------------------------------------------------------------------
# EDF 1A
# ------------------------------------------------------------------------------------------------
#Load the required libraries.
library(skitools)
library(data.table)
library(skidb)
library(plyr)
library(dplyr)
library(MASS)
library(wesanderson)
library(readxl)
library(ggalluvial)
library(readxl)
library(effects)
library(stats)
library(forcats) 
library(ggpubr)
library(ggforce)
library(ComplexHeatmap)
library(Seurat)

#Load the expression matrix and the cluster data.
expDat = readRDS("../data/EDF/EDF1/1A_supExp.rds")
clustDat = readRDS("../data/EDF/EDF1/1A_supClusters.rds")

#Construct the seurat object.
seu = CreateSeuratObject(counts = expDat)
seu@meta.data$seurat_clusters = clustDat

#Generate the Dotplot.
p = DotPlot(seu, features = rownames(expDat), group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

# ------------------------------------------------------------------------------------------------
# EDF 1B & 1C
# ------------------------------------------------------------------------------------------------

#Load required libraries.
library(skitools)
library(Seurat)
library(ggpubr)

#Load seurat metadata (with clustering results and cell names), load per cell aneploidy scores, and load shannon entropy scores per cluster.
tmp.plot = data.table(readRDS("../data/EDF/EDF1/1BC_seuMetaData.rds"))
tmp.l.norm = readRDS("../data/EDF/EDF1/1BC_anneuploidy_zscore.rds")
exp_shannon = readRDS("../data/EDF/EDF1/1BC_clusterEntropy.rds")

#Create table of total aneuploidy score per cell. Add to the metadata table.
aneu.score = data.table(cell = rownames(tmp.l.norm), Aneuplody.Score = apply(tmp.l.norm,1,function(x) sum(abs(x))))
tmp.plot = merge(tmp.plot,aneu.score,by = 'cell')

#Compute aneuploidy score per cluster. 
tmp.df = aggregate(Aneuplody.Score ~ seurat_clusters, data = tmp.plot, FUN = mean)
colnames(tmp.df)[2] = 'Aneuploidy.Score.mean'
#Add diversity (entropy) value per cluster. 
tmp.df$S.Diversity.num = exp_shannon

#Label carcinoma clusters based on high aneuploidy and low entropy.
tmp.df$cat = ifelse(((tmp.df$S.Diversity.num < 10) &(tmp.df$Aneuploidy.Score.mean > 40)),"Carcinoma","Non-Carcinoma")
tmp.nw = tmp.plot[,.(Tissue,seurat_clusters)]
tmp.nw = as.data.frame(tmp.nw)
tmp.dt = tmp.df
tmp.df = as.data.frame(tmp.df)
tmp.df = merge(tmp.df,tmp.nw,by='seurat_clusters')
tmp.df.nw = tmp.df %>%  distinct(seurat_clusters, .keep_all = TRUE)
tmp.df.nw = as.data.table(tmp.df.nw)

#Generate scatterplot for Sup 1B.
p = ggscatter(tmp.df.nw, x = "S.Diversity.num", y = "Aneuploidy.Score.mean",
              color = 'cat',
              repel = TRUE, xlab = "Shannon's Diversity", ylab = 'Aneuploidy Score (mean)')
print(p)
#Link carcinoma and non-carcinoma annotations to each cell.
tmp.plot$cat = ""
for(i in 1:nrow(tmp.df.nw)){
    tmp.plot[tmp.plot$seurat_clusters==tmp.df.nw$seurat_clusters[i],]$cat = tmp.df.nw$cat[i]  
    }
#Organize cells by annotation.
tmp.plot$cat = factor(tmp.plot$cat, levels = c("Non-Carcinoma","Carcinoma"))
tmp.plot = tmp.plot[order(tmp.plot$cat),]                   

#Organize rows in aneuploidy output for heatmap.
tmp.l.nw = tmp.l.norm[match(tmp.plot$cell,rownames(tmp.l.norm)),]

#Set heatmap annotation and color.
row.ha = HeatmapAnnotation(Tissue = tmp.plot$cat, name = 'Tissue', width = unit(2, "cm"), show_legend = TRUE, show_annotation_name = FALSE, col = list(Tissue = c("Carcinoma" = "#F8766D", "Non-Carcinoma" = "#00BA38")), annotation_legend_param = list("Normal", "Tumor Epithelial", "Carcinoma"), which = "row")
library(circlize)
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
#Plot heatmap.
p = Heatmap(tmp.l.nw, name = 'z-score', col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 8), column_names_side = c("bottom"), show_column_names = TRUE, show_row_names = FALSE, column_title_gp = gpar(fontsize = 20, fontface = "bold")) + row.ha
print(p)

# ------------------------------------------------------------------------------------------------
# EDF 2
# ------------------------------------------------------------------------------------------------
#Load requires libraries.
library(skitools)
library(data.tree)

#Load sikkema metadata table. Format data so that each row represents a full hierarchy for a given celltype.
sikkema_meta = read.table("../data/EDF/EDF2/2_sikk_meta_epi.csv",header=T,sep=",")
sikkema_meta = sikkema_meta[,c("ann_level_1","ann_level_2","ann_level_3","ann_level_4","ann_finest_level")]

#Adjust AT0 annotation to conform with Alveolar annotation.
sikkema_meta[sikkema_meta$ann_finest_level=="AT0",]$ann_level_2 = "Alveolar epithelium"
sikkema_meta[sikkema_meta$ann_finest_level=="AT0",]$ann_level_3 = "AT2"
sikkema_meta[sikkema_meta$ann_finest_level=="AT0",]$ann_level_4 = "AT0"

#Populate AT1 and AT2 annotations for level 3 (i.e.: replace curent None value with the cell name).
sikkema_meta[sikkema_meta$ann_finest_level=="AT2",]$ann_level_4 = "AT2"
sikkema_meta[sikkema_meta$ann_finest_level=="AT1",]$ann_level_4 = "AT1"

formatted_data <- apply(sikkema_meta, 1, function(row) {   
  paste0("level1: ", row[1], ", level2: ", row[2], ", level3: ", row[3], 
         ", level4: ", row[4], ", level5: ", row[5], ", finest_level: ", row[6]) 
})

# Convert the formatted data to a data frame (optional but makes it easier to work with).
formatted_data <- data.frame(formatted_data)

# Add a new column 'pathString' by concatenating non-NA values for each row with '/'
sikkema_meta$pathString <- apply(sikkema_meta, 1, function(row) {   
  paste(na.omit(row), collapse = "/")  # Concatenate non-NA values with "/"
})

#Convert cell level hierarchy into graph.
hierarchy_tree <- as.Node(sikkema_meta)
tree_graph <- ToDiagrammeRGraph(hierarchy_tree)
tree_graph <- DiagrammeR::add_global_graph_attrs(
  graph = tree_graph,
  attr = "rankdir",
  value = "LR",   # 'LR' stands for Left-to-Right layout
  attr_type = "graph"
)
tree_graph <- DiagrammeR::add_global_graph_attrs(
  graph = tree_graph,
  attr = "fontsize",
  value = "30",   # Adjust the font size (default is usually smaller)
  attr_type = "node"
)

# Render the left-to-right tree plot.
DiagrammeR::render_graph(tree_graph)

# ------------------------------------------------------------------------------------------------
# EDF 3
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

#Load file with cell metadata and identity call results.
normal_labs = fread("../data/EDF/EDF3/3A_combined_emb_for_normal_for_homogenous_sikk_ep_nw5.csv")
#Filter results based on uncertainty of identity calls. If uncertainty is lower than 0.5, label as unknown.
normal_labs$ann_finest_lev_transferred_label_filtered = ifelse(normal_labs$ann_finest_lev_transfer_uncert < 0.5,normal_labs$ann_finest_lev_transferred_label_unfiltered,"Unknown")
#Filter for identity call results and true identity labels.
normal.mat = normal_labs[,c("V1","ann_finest_lev_transferred_label_filtered","major.ident.f","ref_or_query")]

#Filter for cells used as query (i.e.: not used for training). Remove all cells with unknown labels, then adjust labels based on celltype identity calls.
query_mat = normal.mat[which(normal.mat$ref_or_query=="query"),]
query_mat = query_mat[which(query_mat$ann_finest_lev_transferred_label_filtered!="Unknown"),]
AT2 = c("AT2","AT0","AT2 proliferating")
AT1 = c("AT1")
smg = c("SMG serous (bronchial)","SMG mucous","SMG duct","SMG serous (nasal)")
ciliated = c("Multiciliated (non-nasal)","Multiciliated (nasal)","Deuterosomal")
basal = c("Basal resting","Suprabasal","Hillock-like")
secretory = c("Goblet (nasal)","Club (nasal)","Club (non-nasal)","pre-TB secretory","Goblet (bronchial)","Goblet (subsegmental)")
other = c("Tuft","Ionocyte","Neuroendocrine")
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% AT2,"transfered_labs"]="AT2"
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% AT1,"transfered_labs"]="AT1"
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% smg,"transfered_labs"]="SMG"
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% ciliated,"transfered_labs"]="ciliated"
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% basal,"transfered_labs"]="basal"
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% rare,"transfered_labs"]="other"
query_mat[query_mat$ann_finest_lev_transferred_label_filtered %in% secretory,"transfered_labs"]="secretory"
query_mat[query_mat$major.ident.f %in% c("Alveolar.Type.II.cells"),"orig_labs"]="AT2"
query_mat[query_mat$major.ident.f %in% c("Club.cells"),"orig_labs"]="secretory"
query_mat[query_mat$major.ident.f %in% c("Alveolar.Type.I.cells"),"orig_labs"]="AT1"
query_mat[query_mat$major.ident.f %in% c("Basal.cells"),"orig_labs"]="basal"
query_mat[query_mat$major.ident.f %in% c("Ciliated.cells"),"orig_labs"]="ciliated"
query_mat[query_mat$major.ident.f %in% c("Neuroendocrine.cells"),"orig_labs"]="other"

#Compute accuracy by comparing identity calls vs true identities.
AT2_mat = query_mat[which(query_mat$orig_labs=="AT2"),]
AT2_acc = sum(AT2_mat$transfered_labs==AT2_mat$orig_labs)/nrow(AT2_mat)

AT1_mat = query_mat[which(query_mat$orig_labs=="AT1"),]
AT1_acc = sum(AT1_mat$transfered_labs==AT1_mat$orig_labs)/nrow(AT1_mat)

secretory_mat = query_mat[which(query_mat$orig_labs=="secretory"),]
secretory_acc = sum(secretory_mat$transfered_labs==secretory_mat$orig_labs)/nrow(secretory_mat)

basal_mat = query_mat[which(query_mat$orig_labs=="basal"),]
basal_acc = sum(basal_mat$transfered_labs==basal_mat$orig_labs)/nrow(basal_mat)

rare_mat = query_mat[which(query_mat$orig_labs=="other"),]
rare_acc = sum(rare_mat$transfered_labs==rare_mat$orig_labs)/nrow(rare_mat)

ciliated_mat = query_mat[which(query_mat$orig_labs=="ciliated"),]
ciliated_acc = sum(ciliated_mat$transfered_labs==ciliated_mat$orig_labs)/nrow(ciliated_mat)

#Combine results into a matrix and plot.
acc_mat = c(AT2_acc,AT1_acc,secretory_acc,basal_acc,ciliated_acc,rare_acc)
acc_mat = as.matrix(acc_mat)
rownames(acc_mat) = c("AT2","AT1","secretory","basal","ciliated","other")
acc_mat = as.data.frame(acc_mat)
acc_mat$cat = rownames(acc_mat)

p = ggplot(acc_mat, aes(x = cat, y = V1)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() + coord_flip()
print(p)
                                                                             
# ------------------------------------------------------------------------------------------------
# EDF 3B - Left Panel
# ------------------------------------------------------------------------------------------------
#Load required libraries.
library(skitools)

#Load GRanges with centroid expression data.
juan.genes.gr.nw = readRDS("../data/EDF/EDF3/3BLP_genes_gr_LUAD.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)

#Define quantiles for partitioning the genes by basal resting expression.
sup.quant <- quantile(juan.genes.gr.nw.dt$Basal.resting)
sup.quant.2 <- quantile(juan.genes.gr.nw.dt$Basal.resting, probs = c(0.33, 0.67, 1))

#Define a data.frame starting from the expression information.
sikk_sup <- juan.genes.gr.nw.dt$Basal.resting %>% as.data.frame()
#Add quantile information and gene name.
colnames(sikk_sup)[1] = 'Goblet'
sikk_sup$gene <- juan.genes.gr.nw.dt$gene_name
lowsup.quant <- sikk_sup[which(sikk_sup$Goblet <= -6.888385), ]
midsup.quant <- sikk_sup[which(sikk_sup$Goblet > -6.888385  & sikk_sup$Goblet <= -6.101867) , ]
highsup.quant <- sikk_sup[which(sikk_sup$Goblet > -6.101867) , ]
lowsup.quant$Quart = 'Q1'
midsup.quant$Quart = 'Q2'
highsup.quant$Quart = 'Q3'
sup.quart.mat <- rbind(lowsup.quant,midsup.quant,highsup.quant)

#Load GRanges with snv counts information per gene for LUAD. Convert to mutational density.
mgenes= readRDS(file.path(params$fileLoc,"EDF/EDF3/3B_mgenes.rds"))
mgenes2 <- copy(mgenes)
juan.genes.gr.nw.dt$densityLUSC <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*53)

#Make sure quantile data.frame gene order matches snv counts GRanges gene order. Remove NAs. Then add gene mutational density information for LUAD to the data.frame.
supLUSC_centwrtTMB = sup.quart.mat[match(juan.genes.gr.nw.dt$gene_name,sup.quart.mat$gene),]
supLUSC_centwrtTMB = supLUSC_centwrtTMB[complete.cases(supLUSC_centwrtTMB),]
sup.luscmat = cbind(supLUSC_centwrtTMB,juan.genes.gr.nw.dt$densityLUSC)
colnames(sup.luscmat)[[4]] = "densityLUSC"

#Compute mean and standard deviation values for mutational density per expression quantile.
sup.luscmat.Q1 = sup.luscmat[which(sup.luscmat$Quart=="Q1"),]
sup.luscmat.Q1mean = mean(sup.luscmat.Q1$densityLUSC)
sup.luscmat.Q1se = sd(sup.luscmat.Q1$densityLUSC)/sqrt(length(sup.luscmat.Q1$densityLUSC)) 
sup.luscmat.Q2 = sup.luscmat[which(sup.luscmat$Quart=="Q2"),]
sup.luscmat.Q2mean = mean(sup.luscmat.Q2$densityLUSC)
sup.luscmat.Q2se = sd(sup.luscmat.Q2$densityLUSC)/sqrt(length(sup.luscmat.Q2$densityLUSC))
sup.luscmat.Q3 = sup.luscmat[which(sup.luscmat$Quart=="Q3"),]
sup.luscmat.Q3mean = mean(sup.luscmat.Q3$densityLUSC)
sup.luscmat.Q3se = sd(sup.luscmat.Q3$densityLUSC)/sqrt(length(sup.luscmat.Q3$densityLUSC))
#Save mean and sd values as a new data.frame.
mean_sup_mat_lusc = as.data.frame(rbind(sup.luscmat.Q1mean,sup.luscmat.Q2mean,sup.luscmat.Q3mean))
se_sup_mat_lusc = as.data.frame(rbind(sup.luscmat.Q1se,sup.luscmat.Q2se,sup.luscmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
sup.mean.mat = cbind(quart.mat,mean_sup_mat_lusc,se_sup_mat_lusc)
rownames(sup.mean.mat) = c(1,2,3)
colnames(sup.mean.mat) = c("Quart","Mean_TMB","SE")

#Do the same process for AT2 expression. Define quantiles.
AT2.quant <- quantile(juan.genes.gr.nw.dt$AT2)
AT2.quant.2 <- quantile(juan.genes.gr.nw.dt$AT2, probs = c(0.33, 0.67, 1))

#Define a data.frame starting from the expression information.
sikk_AT2 <- juan.genes.gr.nw.dt$AT2 %>% as.data.frame()
#Add quantile information and gene name.
colnames(sikk_AT2)[1] = 'AT2'
sikk_AT2$gene <- juan.genes.gr.nw.dt$gene_name
lowAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 <= -6.886128), ]
midAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > -6.886128  & sikk_AT2$AT2 <= -5.963367) , ]
highAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > -5.963367) , ]
lowAT2.quant$Quart = 'Q1'
midAT2.quant$Quart = 'Q2'
highAT2.quant$Quart = 'Q3'
AT2.quart.mat <- rbind(lowAT2.quant,midAT2.quant,highAT2.quant)

#Make sure quantile data.frame gene order matches snv counts GRanges gene order. Remove NAs. Then add gene mutational density information for LUAD to the data.frame.
AT2lusc_centwrtTMB = AT2.quart.mat[match(juan.genes.gr.nw.dt$gene_name,AT2.quart.mat$gene),]
AT2lusc_centwrtTMB = AT2lusc_centwrtTMB[complete.cases(AT2lusc_centwrtTMB),]  # remove NAs
AT2.luscmat = cbind(AT2lusc_centwrtTMB,juan.genes.gr.nw.dt$densityLUSC)
colnames(AT2.luscmat)[[4]] = "densityLUSC"

#Compute mean and standard deviation values for mutational density per expression quantile.
AT2.luscmat.Q1 = AT2.luscmat[which(AT2.luscmat$Quart=="Q1"),]
AT2.luscmat.Q1mean = mean(AT2.luscmat.Q1$densityLUSC)
AT2.luscmat.Q1se = sd(AT2.luscmat.Q1$densityLUSC)/sqrt(length(AT2.luscmat.Q1$densityLUSC)) 
AT2.luscmat.Q2 = AT2.luscmat[which(AT2.luscmat$Quart=="Q2"),]
AT2.luscmat.Q2mean = mean(AT2.luscmat.Q2$densityLUSC)
AT2.luscmat.Q2se = sd(AT2.luscmat.Q2$densityLUSC)/sqrt(length(AT2.luscmat.Q2$densityLUSC))
AT2.luscmat.Q3 = AT2.luscmat[which(AT2.luscmat$Quart=="Q3"),]
AT2.luscmat.Q3mean = mean(AT2.luscmat.Q3$densityLUSC)
AT2.luscmat.Q3se = sd(AT2.luscmat.Q3$densityLUSC)/sqrt(length(AT2.luscmat.Q3$densityLUSC))
mean_AT2_mat = as.data.frame(rbind(AT2.luscmat.Q1mean,AT2.luscmat.Q2mean,AT2.luscmat.Q3mean))
#Save mean and sd values as a new data.frame.
se_AT2_mat = as.data.frame(rbind(AT2.luscmat.Q1se,AT2.luscmat.Q2se,AT2.luscmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

#Combine data.frames with mean and standard deviation per expression quartile for AT2 and basal resting cells.
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
AT2.mean.mat = cbind(quart.mat,mean_AT2_mat,se_AT2_mat)
rownames(AT2.mean.mat) = c(1,2,3)
colnames(AT2.mean.mat) = c("Quart","Mean_TMB","SE")
AT2.mean.mat$celltype = "AT2"
sup.mean.mat$celltype="Basalresting"

#Plot the graph.
mean.mat.lusc = rbind(AT2.mean.mat,sup.mean.mat)
mean.mat.lusc$celltype = as.factor(mean.mat.lusc$celltype)
p = ggplot(mean.mat.lusc, aes(x = Quart, y = Mean_TMB, group = celltype, color = celltype)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.2) + 
  geom_errorbar(aes(ymin = Mean_TMB - SE, ymax = Mean_TMB + SE), width = 0.2, size = 0.9, color = 'black') +
  theme(axis.text = element_text(size = 20)) + 
  theme(axis.title = element_text(size = 20)) + 
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

print(p)
                                                                             
# ------------------------------------------------------------------------------------------------
# EDF 3B - Right Panel
# ------------------------------------------------------------------------------------------------
#Load required libraries.
library(skitools)

#Load GRanges with centroid expression data.
juan.genes.gr.nw = readRDS("../data/EDF/EDF3/3BRP_genes_gr_LUSC.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)

#Define quantiles for partitioning the genes by basal resting expression.
sup.quant <- quantile(juan.genes.gr.nw.dt$Basal.resting)
sup.quant.2 <- quantile(juan.genes.gr.nw.dt$Basal.resting, probs = c(0.33, 0.67, 1))

#Define a data.frame starting from the expression information.
sikk_sup <- juan.genes.gr.nw.dt$Basal.resting %>% as.data.frame()
#Add quantile information and gene name.
colnames(sikk_sup)[1] = 'Goblet'
sikk_sup$gene <- juan.genes.gr.nw.dt$gene_name
lowsup.quant <- sikk_sup[which(sikk_sup$Goblet <= -6.888385), ]
midsup.quant <- sikk_sup[which(sikk_sup$Goblet > -6.888385  & sikk_sup$Goblet <= -6.101867) , ]
highsup.quant <- sikk_sup[which(sikk_sup$Goblet > -6.101867) , ]
lowsup.quant$Quart = 'Q1'
midsup.quant$Quart = 'Q2'
highsup.quant$Quart = 'Q3'
sup.quart.mat <- rbind(lowsup.quant,midsup.quant,highsup.quant)

#Load GRanges with snv counts information per gene for LUAD. Convert to mutational density.
mgenes= readRDS(file.path(params$fileLoc,"EDF/EDF4/4AB_mgenes.rds"))
mgenes2 <- copy(mgenes)
juan.genes.gr.nw.dt$densityLUSC <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*53)

#Make sure quantile data.frame gene order matches snv counts GRanges gene order. Remove NAs. Then add gene mutational density information for LUAD to the data.frame.
supLUSC_centwrtTMB = sup.quart.mat[match(juan.genes.gr.nw.dt$gene_name,sup.quart.mat$gene),]
supLUSC_centwrtTMB = supLUSC_centwrtTMB[complete.cases(supLUSC_centwrtTMB),]
sup.luscmat = cbind(supLUSC_centwrtTMB,juan.genes.gr.nw.dt$densityLUSC)
colnames(sup.luscmat)[[4]] = "densityLUSC"

#Compute mean and standard deviation values for mutational density per expression quantile.
sup.luscmat.Q1 = sup.luscmat[which(sup.luscmat$Quart=="Q1"),]
sup.luscmat.Q1mean = mean(sup.luscmat.Q1$densityLUSC)
sup.luscmat.Q1se = sd(sup.luscmat.Q1$densityLUSC)/sqrt(length(sup.luscmat.Q1$densityLUSC)) 
sup.luscmat.Q2 = sup.luscmat[which(sup.luscmat$Quart=="Q2"),]
sup.luscmat.Q2mean = mean(sup.luscmat.Q2$densityLUSC)
sup.luscmat.Q2se = sd(sup.luscmat.Q2$densityLUSC)/sqrt(length(sup.luscmat.Q2$densityLUSC))
sup.luscmat.Q3 = sup.luscmat[which(sup.luscmat$Quart=="Q3"),]
sup.luscmat.Q3mean = mean(sup.luscmat.Q3$densityLUSC)
sup.luscmat.Q3se = sd(sup.luscmat.Q3$densityLUSC)/sqrt(length(sup.luscmat.Q3$densityLUSC))
#Save mean and sd values as a new data.frame.
mean_sup_mat_lusc = as.data.frame(rbind(sup.luscmat.Q1mean,sup.luscmat.Q2mean,sup.luscmat.Q3mean))
se_sup_mat_lusc = as.data.frame(rbind(sup.luscmat.Q1se,sup.luscmat.Q2se,sup.luscmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
sup.mean.mat = cbind(quart.mat,mean_sup_mat_lusc,se_sup_mat_lusc)
rownames(sup.mean.mat) = c(1,2,3)
colnames(sup.mean.mat) = c("Quart","Mean_TMB","SE")

#Do the same process for AT2 expression. Define quantiles.
AT2.quant <- quantile(juan.genes.gr.nw.dt$AT2)
AT2.quant.2 <- quantile(juan.genes.gr.nw.dt$AT2, probs = c(0.33, 0.67, 1))

#Define a data.frame starting from the expression information.
sikk_AT2 <- juan.genes.gr.nw.dt$AT2 %>% as.data.frame()
#Add quantile information and gene name.
colnames(sikk_AT2)[1] = 'AT2'
sikk_AT2$gene <- juan.genes.gr.nw.dt$gene_name
lowAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 <= -6.886128), ]
midAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > -6.886128  & sikk_AT2$AT2 <= -5.963367) , ]
highAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > -5.963367) , ]
lowAT2.quant$Quart = 'Q1'
midAT2.quant$Quart = 'Q2'
highAT2.quant$Quart = 'Q3'
AT2.quart.mat <- rbind(lowAT2.quant,midAT2.quant,highAT2.quant)

#Make sure quantile data.frame gene order matches snv counts GRanges gene order. Remove NAs. Then add gene mutational density information for LUAD to the data.frame.
AT2lusc_centwrtTMB = AT2.quart.mat[match(juan.genes.gr.nw.dt$gene_name,AT2.quart.mat$gene),]
AT2lusc_centwrtTMB = AT2lusc_centwrtTMB[complete.cases(AT2lusc_centwrtTMB),]
AT2.luscmat = cbind(AT2lusc_centwrtTMB,juan.genes.gr.nw.dt$densityLUSC)
colnames(AT2.luscmat)[[4]] = "densityLUSC"

#Compute mean and standard deviation values for mutational density per expression quantile.
AT2.luscmat.Q1 = AT2.luscmat[which(AT2.luscmat$Quart=="Q1"),]
AT2.luscmat.Q1mean = mean(AT2.luscmat.Q1$densityLUSC)
AT2.luscmat.Q1se = sd(AT2.luscmat.Q1$densityLUSC)/sqrt(length(AT2.luscmat.Q1$densityLUSC)) 
AT2.luscmat.Q2 = AT2.luscmat[which(AT2.luscmat$Quart=="Q2"),]
AT2.luscmat.Q2mean = mean(AT2.luscmat.Q2$densityLUSC)
AT2.luscmat.Q2se = sd(AT2.luscmat.Q2$densityLUSC)/sqrt(length(AT2.luscmat.Q2$densityLUSC))
AT2.luscmat.Q3 = AT2.luscmat[which(AT2.luscmat$Quart=="Q3"),]
AT2.luscmat.Q3mean = mean(AT2.luscmat.Q3$densityLUSC)
AT2.luscmat.Q3se = sd(AT2.luscmat.Q3$densityLUSC)/sqrt(length(AT2.luscmat.Q3$densityLUSC))
mean_AT2_mat = as.data.frame(rbind(AT2.luscmat.Q1mean,AT2.luscmat.Q2mean,AT2.luscmat.Q3mean))
#Save mean and sd values as a new data.frame.
se_AT2_mat = as.data.frame(rbind(AT2.luscmat.Q1se,AT2.luscmat.Q2se,AT2.luscmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

#Combine data.frames with mean and standard deviation per expression quartile for AT2 and basal resting cells.
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
AT2.mean.mat = cbind(quart.mat,mean_AT2_mat,se_AT2_mat)
rownames(AT2.mean.mat) = c(1,2,3)
colnames(AT2.mean.mat) = c("Quart","Mean_TMB","SE")
AT2.mean.mat$celltype = "AT2"
sup.mean.mat$celltype="Basalresting"

#Plot the graph.
mean.mat.lusc = rbind(AT2.mean.mat,sup.mean.mat)
mean.mat.lusc$celltype = as.factor(mean.mat.lusc$celltype)
p = ggplot(mean.mat.lusc, aes(x = Quart, y = Mean_TMB, group = celltype, color = celltype)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.2) + 
  geom_errorbar(aes(ymin = Mean_TMB - SE, ymax = Mean_TMB + SE), width = 0.2, size = 0.9, color = 'black') +
  theme(axis.text = element_text(size = 20)) + 
  theme(axis.title = element_text(size = 20)) + 
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

print(p)

# ------------------------------------------------------------------------------------------------
# EDF 3C - Left Panel
# ------------------------------------------------------------------------------------------------
#Load libraries
library(skitools)

#Load GRanges with mutational counts per patient.
#The complete GRanges is too heavy for Github, so load the different parts as data.table, join and convert top GRanges.
dt1 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP1.rds")
dt2 = readRDS("../data/EDF/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP2.rds")
dt3 = readRDS("../data/EDF/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP3.rds")
dt = dt2gr(rbind(dt1,dt2,dt3))

#Gather columns with tobacco counts per patient. Save with pair ID as data.table.
count.int = grep("snv.tob.count.", colnames(mcols(dt)), value = TRUE)
count.int = as.data.frame(count.int)
count.int$pair = gsub("snv.tob.count.","",count.int$count.int)
count.int$pair = toupper(count.int$pair)

#Load list of luad_pairs. 
luad_pairs = readRDS("../data/EDF/EDF5/5ABC_new.pairs.luad.rds")
juan.genes.dt = gr2dt(dt)
juan.genes.dt = as.data.frame(juan.genes.dt)

#Filter GRanges for luad pairs. Gather sum of tobacco mutations per gene for LUADs.
LUAD_patnt_mat = intersect(luad_pairs,count.int$pair)
LUAD_patnt_mat = as.matrix(LUAD_patnt_mat)
LUAD_patnt_mat = LUAD_patnt_mat[complete.cases(LUAD_patnt_mat),]
count.int.mat = count.int[match(LUAD_patnt_mat,count.int$pair),]
juan.genes.dt_nw = juan.genes.dt[,match(count.int.mat$count.int,colnames(juan.genes.dt))]
snv_count = as.matrix(rowSums(juan.genes.dt_nw))
juan_dt1 = juan.genes.dt[,1:6]
juan_dt1$snv.count = snv_count
colnames(juan_dt1)[7] = "snv.count"

#Load centroid expression per gene. Append average of centroid expression per gene to data.frame with LUAD tobacco counts per gene.
sikk.cent = readRDS("../data/EDF/EDF5/5ABC_epcells.rds")
sikk.cent_nw = sikk.cent[match(juan_dt1$gene_name,rownames(sikk.cent)),]
sikk.cent_nw = as.data.frame(rowMeans(sikk.cent_nw))
juan_dt2 = cbind(juan_dt1,sikk.cent_nw)
colnames(juan_dt2)[7] = "snv.count"
colnames(juan_dt2)[8] = "exp"

#Compute expression quantiles for the gene expression distribution, and assign quantile information for each gene.
ep.quant <- quantile(juan_dt2$exp)
sikk_ep <- juan_dt2$exp %>% as.data.frame()
colnames(sikk_ep)[1] = 'Epcells'
sikk_ep$gene <- juan_dt2$gene_name
lowep.quant <- sikk_ep[which(sikk_ep$Epcells <= -6.810384), ]
midep.quant <- sikk_ep[which(sikk_ep$Epcells > -6.810384  & sikk_ep$Epcells <= -5.882191) , ]
highep.quant <- sikk_ep[which(sikk_ep$Epcells > -5.882191) , ]
lowep.quant$Quart = 'Q1'
midep.quant$Quart = 'Q2'
highep.quant$Quart = 'Q3'
ep.quart.mat <- rbind(lowep.quant,midep.quant,highep.quant)

#Per expression quantile (Q1, Q2, and Q3), compute mean and standard deviation of gene mutational density. 
juan_dt2$densityLUAD <- (10**6)*juan_dt2$snv.count/(juan_dt2$width*246)
epLUAD_centwrtTMB = ep.quart.mat[match(juan_dt2$gene_name,ep.quart.mat$gene),]
epLUAD_centwrtTMB = epLUAD_centwrtTMB[complete.cases(epLUAD_centwrtTMB),]  # remove NAs
ep.luadmat = cbind(epLUAD_centwrtTMB,juan_dt2$densityLUAD)
colnames(ep.luadmat)[[4]] = "densityLUAD"
ep.luadmat.Q1 = ep.luadmat[which(ep.luadmat$Quart=="Q1"),]
ep.luadmat.Q1mean = mean(ep.luadmat.Q1$densityLUAD)
ep.luadmat.Q1se = sd(ep.luadmat.Q1$densityLUAD)/sqrt(length(ep.luadmat.Q1$densityLUAD)) 
ep.luadmat.Q2 = ep.luadmat[which(ep.luadmat$Quart=="Q2"),]
ep.luadmat.Q2mean = mean(ep.luadmat.Q2$densityLUAD)
ep.luadmat.Q2se = sd(ep.luadmat.Q2$densityLUAD)/sqrt(length(ep.luadmat.Q2$densityLUAD))
ep.luadmat.Q3 = ep.luadmat[which(ep.luadmat$Quart=="Q3"),]
ep.luadmat.Q3mean = mean(ep.luadmat.Q3$densityLUAD)
ep.luadmat.Q3se = sd(ep.luadmat.Q3$densityLUAD)/sqrt(length(ep.luadmat.Q3$densityLUAD))
mean_ep_mat = as.data.frame(rbind(ep.luadmat.Q1mean,ep.luadmat.Q2mean,ep.luadmat.Q3mean))
se_ep_mat = as.data.frame(rbind(ep.luadmat.Q1se,ep.luadmat.Q2se,ep.luadmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

#Generate new data.frame with mean and sd of mutational density per quantile.
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
ep.mean.mat = cbind(quart.mat,mean_ep_mat,se_ep_mat)
rownames(ep.mean.mat) = c(1,2,3)
colnames(ep.mean.mat) = c("Quart","Mean_TMB","SE")
ep.mean.mat$celltype="SBS4"

#Generate plot of mutational density per expression quantile.
mean.mat.luad = ep.mean.mat
mean.mat.luad$celltype = as.factor(mean.mat.luad$celltype)
p = ggplot(mean.mat.luad, aes(x = Quart, y = Mean_TMB, group = celltype, color = celltype)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = Mean_TMB - SE, ymax = Mean_TMB + SE), width = 0.2, size = 0.9, color = 'black') +
  theme(axis.text = element_text(size = 20)) + 
  theme(axis.title = element_text(size = 20)) + 
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )+ labs(x = "Expression", y = "Mean TMB for tobacco")+ theme(legend.position = "none")

print(p)
# ------------------------------------------------------------------------------------------------
# EDF 3C - Center Panel
# ------------------------------------------------------------------------------------------------
#Load libraries
library(skitools)

#Load GRanges with mutational counts per patient.
#The complete GRanges is too heavy for Github, so load the different parts as data.table, join and convert top GRanges.
dt1 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP1.rds")
dt2 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP2.rds")
dt3 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP3.rds")
dt = dt2gr(rbind(dt1,dt2,dt3))

#Gather columns with aging counts per patient. Save with pair ID as data.table.
count.int = grep("snv.aging.count.", colnames(mcols(dt)), value = TRUE)
count.int = as.data.frame(count.int)
count.int$pair = gsub("snv.aging.count.","",count.int$count.int)
count.int$pair = toupper(count.int$pair)

#Load list of luad_pairs. 
luad_pairs = readRDS("../data/EDF/EDF3/3C_new.pairs.luad.rds")
juan.genes.dt = gr2dt(dt)
juan.genes.dt = as.data.frame(juan.genes.dt)

#Filter GRanges for luad pairs. Gather sum of tobacco mutations per gene for LUADs.
LUAD_patnt_mat = intersect(luad_pairs,count.int$pair)
LUAD_patnt_mat = as.matrix(LUAD_patnt_mat)
LUAD_patnt_mat = LUAD_patnt_mat[complete.cases(LUAD_patnt_mat),]
count.int.mat = count.int[match(LUAD_patnt_mat,count.int$pair),]
juan.genes.dt_nw = juan.genes.dt[,match(count.int.mat$count.int,colnames(juan.genes.dt))]
snv_count = as.matrix(rowSums(juan.genes.dt_nw))
juan_dt1 = juan.genes.dt[,1:6]
juan_dt1$snv.count = snv_count
colnames(juan_dt1)[7] = "snv.count"

#Load centroid expression per gene. Append average of centroid expression per gene to data.frame with LUAD tobacco counts per gene.
sikk.cent = readRDS("../data/EDF/EDF3/3C_epcells.rds")
sikk.cent_nw = sikk.cent[match(juan_dt1$gene_name,rownames(sikk.cent)),]
sikk.cent_nw = as.data.frame(rowMeans(sikk.cent_nw))
juan_dt2 = cbind(juan_dt1,sikk.cent_nw)
colnames(juan_dt2)[7] = "snv.count"
colnames(juan_dt2)[8] = "exp"

#Compute expression quantiles for the gene expression distribution, and assign quantile information for each gene.
ep.quant <- quantile(juan_dt2$exp)
sikk_ep <- juan_dt2$exp %>% as.data.frame()
colnames(sikk_ep)[1] = 'Epcells'
sikk_ep$gene <- juan_dt2$gene_name
lowep.quant <- sikk_ep[which(sikk_ep$Epcells <= -6.810384), ]
midep.quant <- sikk_ep[which(sikk_ep$Epcells > -6.810384  & sikk_ep$Epcells <= -5.882191) , ]
highep.quant <- sikk_ep[which(sikk_ep$Epcells > -5.882191) , ]
lowep.quant$Quart = 'Q1'
midep.quant$Quart = 'Q2'
highep.quant$Quart = 'Q3'
ep.quart.mat <- rbind(lowep.quant,midep.quant,highep.quant)

#Per expression quantile (Q1, Q2, and Q3), compute mean and standard deviation of gene mutational density. 
juan_dt2$densityLUAD <- (10**6)*juan_dt2$snv.count/(juan_dt2$width*246)
epLUAD_centwrtTMB = ep.quart.mat[match(juan_dt2$gene_name,ep.quart.mat$gene),]
epLUAD_centwrtTMB = epLUAD_centwrtTMB[complete.cases(epLUAD_centwrtTMB),]  # remove NAs
ep.luadmat = cbind(epLUAD_centwrtTMB,juan_dt2$densityLUAD)
colnames(ep.luadmat)[[4]] = "densityLUAD"
ep.luadmat.Q1 = ep.luadmat[which(ep.luadmat$Quart=="Q1"),]
ep.luadmat.Q1mean = mean(ep.luadmat.Q1$densityLUAD)
ep.luadmat.Q1se = sd(ep.luadmat.Q1$densityLUAD)/sqrt(length(ep.luadmat.Q1$densityLUAD)) 
ep.luadmat.Q2 = ep.luadmat[which(ep.luadmat$Quart=="Q2"),]
ep.luadmat.Q2mean = mean(ep.luadmat.Q2$densityLUAD)
ep.luadmat.Q2se = sd(ep.luadmat.Q2$densityLUAD)/sqrt(length(ep.luadmat.Q2$densityLUAD))
ep.luadmat.Q3 = ep.luadmat[which(ep.luadmat$Quart=="Q3"),]
ep.luadmat.Q3mean = mean(ep.luadmat.Q3$densityLUAD)
ep.luadmat.Q3se = sd(ep.luadmat.Q3$densityLUAD)/sqrt(length(ep.luadmat.Q3$densityLUAD))
mean_ep_mat = as.data.frame(rbind(ep.luadmat.Q1mean,ep.luadmat.Q2mean,ep.luadmat.Q3mean))
se_ep_mat = as.data.frame(rbind(ep.luadmat.Q1se,ep.luadmat.Q2se,ep.luadmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

#Generate new data.frame with mean and sd of mutational density per quantile.
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
ep.mean.mat = cbind(quart.mat,mean_ep_mat,se_ep_mat)
rownames(ep.mean.mat) = c(1,2,3)
colnames(ep.mean.mat) = c("Quart","Mean_TMB","SE")
ep.mean.mat$celltype="Aging"

#Generate plot of mutational density per expression quantile.
mean.mat.luad = ep.mean.mat
mean.mat.luad$celltype = as.factor(mean.mat.luad$celltype)
p = ggplot(mean.mat.luad, aes(x = Quart, y = Mean_TMB, group = celltype, color = celltype)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = Mean_TMB - SE, ymax = Mean_TMB + SE), width = 0.2, size = 0.9, color = 'black') +
  theme(axis.text = element_text(size = 20)) + 
  theme(axis.title = element_text(size = 20)) + 
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )+ labs(x = "Expression", y = "Mean TMB for aging")+ theme(legend.position = "none")

print(p)
# ------------------------------------------------------------------------------------------------
# EDF 3C - Right Panel
# ------------------------------------------------------------------------------------------------
#Load libraries
library(skitools)

#Load GRanges with mutational counts per patient.
#The complete GRanges is too heavy for Github, so load the different parts as data.table, join and convert top GRanges.
dt1 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP1.rds")
dt2 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP2.rds")
dt3 = readRDS("../data/EDF/EDF3/3C_DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExprP3.rds")
dt = dt2gr(rbind(dt1,dt2,dt3))

#Gather columns with apobec counts per patient. Save with pair ID as data.table.
count.int = grep("snv.apobec.count.", colnames(mcols(dt)), value = TRUE)
count.int = as.data.frame(count.int)
count.int$pair = gsub("snv.apobec.count.","",count.int$count.int)
count.int$pair = toupper(count.int$pair)

#Load list of luad_pairs. 
luad_pairs = readRDS("../data/EDF/EDF5/5ABC_new.pairs.luad.rds"))
juan.genes.dt = gr2dt(dt)
juan.genes.dt = as.data.frame(juan.genes.dt)

#Filter GRanges for luad pairs. Gather sum of tobacco mutations per gene for LUADs.
LUAD_patnt_mat = intersect(luad_pairs,count.int$pair)
LUAD_patnt_mat = as.matrix(LUAD_patnt_mat)
LUAD_patnt_mat = LUAD_patnt_mat[complete.cases(LUAD_patnt_mat),]
count.int.mat = count.int[match(LUAD_patnt_mat,count.int$pair),]
juan.genes.dt_nw = juan.genes.dt[,match(count.int.mat$count.int,colnames(juan.genes.dt))]
snv_count = as.matrix(rowSums(juan.genes.dt_nw))
juan_dt1 = juan.genes.dt[,1:6]
juan_dt1$snv.count = snv_count
colnames(juan_dt1)[7] = "snv.count"

#Load centroid expression per gene. Append average of centroid expression per gene to data.frame with LUAD tobacco counts per gene.
sikk.cent = readRDS("../data/EDF/EDF5/5ABC_epcells.rds"))
sikk.cent_nw = sikk.cent[match(juan_dt1$gene_name,rownames(sikk.cent)),]
sikk.cent_nw = as.data.frame(rowMeans(sikk.cent_nw))
juan_dt2 = cbind(juan_dt1,sikk.cent_nw)
colnames(juan_dt2)[7] = "snv.count"
colnames(juan_dt2)[8] = "exp"

#Compute expression quantiles for the gene expression distribution, and assign quantile information for each gene.
ep.quant <- quantile(juan_dt2$exp)
sikk_ep <- juan_dt2$exp %>% as.data.frame()
colnames(sikk_ep)[1] = 'Epcells'
sikk_ep$gene <- juan_dt2$gene_name
lowep.quant <- sikk_ep[which(sikk_ep$Epcells <= -6.810384), ]
midep.quant <- sikk_ep[which(sikk_ep$Epcells > -6.810384  & sikk_ep$Epcells <= -5.882191) , ]
highep.quant <- sikk_ep[which(sikk_ep$Epcells > -5.882191) , ]
lowep.quant$Quart = 'Q1'
midep.quant$Quart = 'Q2'
highep.quant$Quart = 'Q3'
ep.quart.mat <- rbind(lowep.quant,midep.quant,highep.quant)

#Per expression quantile (Q1, Q2, and Q3), compute mean and standard deviation of gene mutational density. 
juan_dt2$densityLUAD <- (10**6)*juan_dt2$snv.count/(juan_dt2$width*246)
epLUAD_centwrtTMB = ep.quart.mat[match(juan_dt2$gene_name,ep.quart.mat$gene),]
epLUAD_centwrtTMB = epLUAD_centwrtTMB[complete.cases(epLUAD_centwrtTMB),]  # remove NAs
ep.luadmat = cbind(epLUAD_centwrtTMB,juan_dt2$densityLUAD)
colnames(ep.luadmat)[[4]] = "densityLUAD"
ep.luadmat.Q1 = ep.luadmat[which(ep.luadmat$Quart=="Q1"),]
ep.luadmat.Q1mean = mean(ep.luadmat.Q1$densityLUAD)
ep.luadmat.Q1se = sd(ep.luadmat.Q1$densityLUAD)/sqrt(length(ep.luadmat.Q1$densityLUAD)) 
ep.luadmat.Q2 = ep.luadmat[which(ep.luadmat$Quart=="Q2"),]
ep.luadmat.Q2mean = mean(ep.luadmat.Q2$densityLUAD)
ep.luadmat.Q2se = sd(ep.luadmat.Q2$densityLUAD)/sqrt(length(ep.luadmat.Q2$densityLUAD))
ep.luadmat.Q3 = ep.luadmat[which(ep.luadmat$Quart=="Q3"),]
ep.luadmat.Q3mean = mean(ep.luadmat.Q3$densityLUAD)
ep.luadmat.Q3se = sd(ep.luadmat.Q3$densityLUAD)/sqrt(length(ep.luadmat.Q3$densityLUAD))
mean_ep_mat = as.data.frame(rbind(ep.luadmat.Q1mean,ep.luadmat.Q2mean,ep.luadmat.Q3mean))
se_ep_mat = as.data.frame(rbind(ep.luadmat.Q1se,ep.luadmat.Q2se,ep.luadmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

#Generate new data.frame with mean and sd of mutational density per quantile.
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
ep.mean.mat = cbind(quart.mat,mean_ep_mat,se_ep_mat)
rownames(ep.mean.mat) = c(1,2,3)
colnames(ep.mean.mat) = c("Quart","Mean_TMB","SE")
ep.mean.mat$celltype="Apobec"

#Generate plot of mutational density per expression quantile.
mean.mat.luad = ep.mean.mat
mean.mat.luad$celltype = as.factor(mean.mat.luad$celltype)
p = ggplot(mean.mat.luad, aes(x = Quart, y = Mean_TMB, group = celltype, color = celltype)) + 
  geom_point(size = 4) + 
  geom_line(size = 1.2) +
  geom_errorbar(aes(ymin = Mean_TMB - SE, ymax = Mean_TMB + SE), width = 0.2, size = 0.9, color = 'black') +
  theme(axis.text = element_text(size = 20)) + 
  theme(axis.title = element_text(size = 20)) + 
  theme(
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )+ labs(x = "Expression", y = "Mean TMB for apobec")+ theme(legend.position = "none")

print(p)

# ------------------------------------------------------------------------------------------------
# EDF 3D
# ------------------------------------------------------------------------------------------------
#Load libraries
library(skitools)

# Accuracy for predicting true cell type at different levels 
sims_acc_results <- fread('../data/EDF/EDF3/3D_luad_sims_results.csv')

# accuracy at level finest
p1_5 <- ggplot(sims_acc_results[, .(celltype, accuracy_level_finest)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level_finest)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,1) + 
  ggtitle("Finest level match") +
  xlab("True cell type") + ylab("Accuracy")

# accuracy at level 4 resolution
p1_4 <- ggplot(sims_acc_results[, .(celltype, accuracy_level4)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level4)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Level 4 match") +
  xlab("True cell type") + ylab("Accuracy")

# accuracy at level 3 resolution
p1_3 <- ggplot(sims_acc_results[, .(celltype, accuracy_level3)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level3)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Level 3 match") +
  xlab("True cell type") + ylab("Accuracy")

# accuracy at level 2 resolution
p1_2 <- ggplot(sims_acc_results[, .(celltype, accuracy_level2)] %>% as.data.frame(), aes(x = celltype, y = accuracy_level2)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Level 2 match") +
  xlab("True cell type") + ylab("Accuracy")

# accuracy at level distal - non - distal resolution
p1_1 <- ggplot(sims_acc_results[, .(celltype, accuracy_leveldnd)] %>% as.data.frame(), aes(x = celltype, y = accuracy_leveldnd)) +
  geom_point(size=2, shape=23, color = 'darkgreen') + 
  theme_classic() +
  theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + ylim(0,1) + 
  ggtitle("Distal vs non-distal match") +
  xlab("True cell type") + ylab("Accuracy")

p_acc <- ggpubr::ggarrange(p1_1, p1_2, p1_3, p1_4, p1_5, ncol = 3, nrow = 2)
print(p_acc)

# ------------------------------------------------------------------------------------------------
# EDF 3E
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

#Plots can be made for univariate and multivariate simultaneously or individually. Comments are made when relevant to do one or the other.

#Load simulation result files. Add Method column to identify if results are uni or multivariate.
resDT = readRDS("../data/EDF/EDF3/3E_testRun_withOverlap_Cov.rds")
resDT$Method = "Univariate"

#If results are plotted in the same graph, uncomment the following lines, adjusting the Method column value accordingly. 
#MulDT = readRDS(file.path(params$fileLoc,"EDF/EDF3/3E_testRun_withOverlapMultivariatewithCovs.rds"))
# MulDT$Method = "Multivariate"
#resDT = rbind(resDT,MulDT,fill=TRUE)

#Load csv with Sikkema annotations.
annotSikk = read.csv("../data/EDF/EDF3/3E_sikk_meta_epi3.csv")
#Add proximal and distal annotations.
annotSikk$ann_level_PD = ifelse(grepl("^AT",annotSikk$ann_finest_level),"Proximal","Distal")

#Pick a level (i.e.: column).
#Check the overlapping cell calls. Are all the overlapping cells from the same level? If not just put ambiguous. If yes use the level label.
resDT$Overlap_Level = ""
resDT$True_Level = ""
#Define a column to determine if sim had the right call (1) or not (0).
resDT$Lev_Hit = 0

#For the true cell just add the corresponding level value.
#Check if the 2 correspond. If yes add 1 to a new level accuracy column. If not just add a 0.  
for(i in 1:nrow(resDT)){
#    if(i%%100==0){
#        print(i)
#     }
    resDT$True_Level[i] = annotSikk[annotSikk$ann_finest_level==resDT$thisCell[i],]$ann_level_PD
    calledLevels = unique(annotSikk[annotSikk$ann_finest_level %in% strsplit(resDT$allOverlap[i],",")[[1]],]$ann_level_PD)
    if(length(calledLevels)>1){
        resDT$Overlap_Level[i] = "Ambiguous"
     }else if(length(calledLevels)==0){
         resDT$Overlap_Level[i] = "No call"
    }else{
         resDT$Overlap_Level[i] = calledLevels
     }
    resDT$Lev_Hit[i] = ifelse(resDT$Overlap_Level[i] == resDT$True_Level[i],1,0)
}

#Get average accuracy per TMB using the Lev_Hit column (accuracy = number of 1s in said column vs number of results for that TMB).
resDT$Lev_Hit = as.numeric(resDT$Lev_Hit)

#Results are grouped eithet by Method and TMB or only TMB depending on the desired graph.
#average_df = resDT %>% group_by(Method,TMB) %>% summarise(average_value = mean(Lev_Hit, na.rm = TRUE))
average_df = resDT %>% group_by(TMB) %>% dplyr::summarise(average_value = mean(Lev_Hit, na.rm = TRUE))  
average_df = as.data.frame(average_df)  

#Line graph is plotted. Depending on whether both uni and multivariate are plotted use the first line or the second (comment/uncomment adequately).
p =  ggplot(data=average_df, aes(x=log10(TMB), y=average_value)) + geom_line(color="blue")+ geom_point() + ggtitle("Accuracy per TMB Multivariate - Proximal vs Distal" ) + xlab("Log10(TMB) (mut/mb)") + ylab("Accuracy") + scale_y_continuous(limits=c(0,1)) +   theme(text = element_text(size = 20)) 
#p =  ggplot(data=average_df, aes(x=log10(TMB), y=average_value, color = Method, group=Method)) + geom_line(color="red")+ geom_point(aes(shape=Method),size=2) + ggtitle("Accuracy per TMB - Celltype" ) + xlab("Log10(TMB) (mut/mb)") + ylab("Accuracy") + scale_y_continuous(limits=c(0,1)) +   theme(text = element_text(size = 20)) +  scale_color_manual(values = c("Univariate" = "red", "Multivariate" = "blue")) 
print(p)

# ------------------------------------------------------------------------------------------------
# EDF 4A
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

# Heat map for LUSC 
nsclc_col <- readRDS('../data/EDF/EDF4/4A_nsclc_hm_col.rds')
hm_lusc <- fread('../data/EDF/EDF4/4A_hm_lusc.csv')
hm_lusc <- hm_lusc %>% as.data.frame()
column_ha2 = HeatmapAnnotation(
  sbs4 = anno_barplot(hm_lusc[,25], gp = gpar(col = "red", fill = "#FF0000")),
  sbs1 = anno_barplot( hm_lusc[,26], gp = gpar(col = "green", fill = "#99FF33")) ,
  sbs5 = anno_barplot(hm_lusc[,27], gp = gpar(col = "green", fill = "#99FF33")) ,
  sbs2 = anno_barplot(hm_lusc[,28], gp = gpar(col = "orange", fill = "#FF8000")) ,
  sbs13 = anno_barplot(hm_lusc[,29], gp = gpar(col = "orange", fill = "#FF8000")) ,
  id1 = anno_barplot(hm_lusc[,31], gp = gpar(col = "green", fill = "#99FF33")) ,
  id2 = anno_barplot(hm_lusc[,32], gp = gpar(col = "violet", fill = "#9900CC")) ,
  id3 = anno_barplot(hm_lusc[,33], gp = gpar(col = "red", fill = "#FF0000")) ,
  id12 = anno_barplot(hm_lusc[,34], gp = gpar(col = "violet", fill = "#9900CC")) ,
  smoker = hm_lusc[,35],
  col = list(
    smoker = c("Smoker" = "orange", "Never Smoker" = "blue") 
  )
)

set.seed(90210)
Heatmap(t(hm_lusc[,1:23]), name = "Relative Risk", col = nsclc_col, cluster_rows = TRUE, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 15), column_names_gp = gpar(fontsize = 10), 
            column_names_side = c("bottom"), show_column_names = FALSE, column_km = 4, column_km_repeats = 100, top_annotation = column_ha2,
            show_parent_dend_line = FALSE, column_gap = unit(c(4), "mm"), column_title = NULL, border = TRUE)
# ------------------------------------------------------------------------------------------------
# EDF 4B
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)
# LUSC - Cluster 1 Distal Lung
edf_8_data <- readRDS('../data/EDF/EDF4/4B_data.rds')
ggplot(edf_8_data[cluster == 'Distal'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0('LUSC Cluster 1')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# LUSC - Cluster 2 Proximal 1
ggplot(edf_8_data[cluster == 'Proximal_1'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0('LUSC Cluster 2')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# LUSC - Cluster 3  Proximal_2
ggplot(edf_8_data[cluster == 'Proximal_2'], aes(x = reorder(celltype,estimate), y = estimate, fill = Cell_Class)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("grey")) +
  geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.65, width = 0.25, color = 'black') +
  scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(-0.25, 0.35)) + # LUAD
  scale_fill_manual("legend", values = c("Proximal" = "olivedrab3", "Distal" = "darkgoldenrod3")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Relative Risk") +
  ggtitle(paste0('LUSC Cluster 3')) +
  theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) + guides(fill=guide_legend(title="Cell types")) +
  coord_flip()

# ------------------------------------------------------------------------------------------------
# EDF 5A
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

pts_coo_id_edf9a <- readRDS('../data/EDF/EDF5/5A_pts_coo_id.rds')
pts_coo_id_edf9a$Lineage_plasticity <- ''
pts_coo_id_edf9a[Identity == 'Distal Lung' & Origin == 'Distal Lung' ]$Lineage_plasticity <- 'Lineage conserved'
pts_coo_id_edf9a[Identity == 'Proximal Lung' & Origin == 'Distal Lung' ]$Lineage_plasticity <- 'Lineage plasticity'
pts_coo_id_edf9a[Identity == 'Distal Lung' & Origin == 'Non-Distal Lung' ]$Lineage_plasticity <- 'Lineage plasticity'
pts_coo_id_edf9a[Identity == 'Proximal Lung' & Origin == 'Non-Distal Lung' ]$Lineage_plasticity <- 'Lineage conserved'

res.plot =  pts_coo_id_edf9a[, prop.test(sum(Lineage_plasticity == 'Lineage plasticity'), .N) %>% dflm %>% cbind(nprox = sum(Lineage_plasticity == 'Lineage plasticity'), tot = .N), by = .(TP53_mut = ifelse(TP53_mut, 'TP53 MUT', 'WT'), Origin)][, fracprox := estimate]
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
  guides(fill = guide_legend(title = 'Fig 4C - Lineage plasticity fraction')) + theme(legend.position = "bottom")

# ------------------------------------------------------------------------------------------------
# EDF 5B
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)
library(ggsankey)

tmp.plot.luad.df3 <- readRDS('../data/EDF/EDF5/5B_data.rds')

tmp.plot.luad.pl <- ggplot(tmp.plot.luad.df3, aes(x = x,                        
                                                  next_x = next_x,                                     
                                                  node = node,
                                                  next_node = next_node,        
                                                  fill = factor(node),
                                                  label = paste0(node, " = ", n)))             # This Creates a label for each node

tmp.plot.luad.pl <- tmp.plot.luad.pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                                                  node.color = "black",     # This is your node color        
                                                  show.legend = TRUE)        # This determines if you want your legend to show

tmp.plot.luad.pl <- tmp.plot.luad.pl + geom_sankey_label(Size = 3,
                                                         color = "black", 
                                                         fill = "white") # This specifies the Label format for each node 
tmp.plot.luad.pl <- tmp.plot.luad.pl + theme_bw()
tmp.plot.luad.pl <- tmp.plot.luad.pl + theme(legend.position = 'none')
tmp.plot.luad.pl <- tmp.plot.luad.pl + theme(axis.title = element_blank(),
                                             axis.text.y = element_blank(),
                                             axis.ticks = element_blank(),
                                             panel.grid = element_blank())
tmp.plot.luad.pl <- tmp.plot.luad.pl + labs(title = "Origin - Identity - Histology")
tmp.plot.luad.pl <- tmp.plot.luad.pl + labs(subtitle = "LUAD")
tmp.plot.luad.pl <- tmp.plot.luad.pl + labs(fill = 'Nodes')
tmp.plot.luad.pl
# ------------------------------------------------------------------------------------------------
# EDF 5C
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

edf9 <- readRDS('../data/EDF/EDF5/5CDEF_data.rds')

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
# EDF 5D
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

# by Identity 
edf9 <- readRDS('../data/EDF/EDF5/5CDEF_data.rds')

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
# EDF 5E
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

# NSCLC-NOS vs TP53
edf9 <- readRDS('../data/EDF/EDF5/5CDEF_data.rds')
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
# EDF 5F
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)

edf9 <- readRDS('../data/EDF/EDF5/5CDEF_data.rds')

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

# ------------------------------------------------------------------------------------------------
# EDF 5G
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)
#Load oncoprint data and plot.
edf5g_a <- readRDS('../data/EDF/EDF5/5G_dataA.rds')
edf5g_b <- readRDS('../data/EDF/EDF5/5G_dataB.rds')

oncoprint(edf5g_b, oncotab = edf5g_a,  genes = c('NKX2-1','SMARCA4','STK11','APC','KEAP1','ALK','MCL1','MAP2K1','EGFR','FGFR1','FOXP1','CDK6','BCL2L1', 'TP53','RB1','CDKN2A','TERT','MYC','ARID1A','PTEN','CCND1','PIK3CA','ERBB2','CCNE1','NF1'), sort.genes = FALSE, colnames.fontsize = 20, 
          rownames.fontsize = 15, signature.thresh = 20, track.height = 1.5, split.gap = 0.5, wes = TRUE, track.gap = track.height/2, drop = TRUE, drop.genes = TRUE, svevents = FALSE,  sort.tumors = FALSE, return.oncotab = FALSE, height = 12, width = 25,  filename = 'EDF9g_onco.pdf')

# ------------------------------------------------------------------------------------------------
# EDF 6A
# ------------------------------------------------------------------------------------------------
#Load libraries
library(skitools)
library(parallel)
library(MASS)

#Load csv with relative risk information per centroid for every sample.
rel_risk = read.table(file.path(params$fileLoc,"EDF/EDF6/6A_luad_rel_risk.csv"),sep=",",header=T,row.names=1)

#Group samples based on Distal, Proximal, or Ambiguous COO calls. Compute mean relative risk for each.
#Distal
distal = rel_risk[which(rel_risk$infer_coo_non_apobec_sigp=="Distal Lung"),]
distal_ep = distal[,1:23]
distal_ep = t(distal_ep)
distal_ep_mean = as.matrix(rowMeans(distal_ep))
distal_ep_mean = as.data.frame(distal_ep_mean)
distal_ep_mean$celltype = rownames(distal_ep_mean)
#Proximal
proximal = rel_risk[which(rel_risk$infer_coo_non_apobec_sigp=="Proximal Lung"),]
proximal_ep = proximal[,1:23]
proximal_ep = t(proximal_ep)
proximal_ep_mean = as.matrix(rowMeans(proximal_ep))
proximal_ep_mean = as.data.frame(proximal_ep_mean)
proximal_ep_mean$celltype = rownames(proximal_ep_mean)
#Ambiguous
ambi = rel_risk[which(rel_risk$infer_coo_non_apobec_sigp=="Ambiguous"),]
ambi_ep = ambi[,1:23]
ambi_ep = t(ambi_ep)
ambi_ep_mean = as.matrix(rowMeans(ambi_ep))
ambi_ep_mean = as.data.frame(ambi_ep_mean)
ambi_ep_mean$celltype = rownames(ambi_ep_mean)

#Load WCM-1 GLM output. Get vector with relative risks per centroid.
wcm = readRDS(file.path(params$fileLoc,"EDF/EDF6/6A_cov_res.rds"))
wcm_mat = wcm[,c("celltype","estimate")]
wcm_mat = as.data.frame(wcm_mat)
rownames(wcm_mat) = wcm_mat$celltype
wcm_mat1 = as.data.frame(wcm_mat[,2])
rownames(wcm_mat1) = rownames(wcm_mat)
colnames(wcm_mat1) = "estimate"

#Compute euclidean distances between average relative risk vectors of distal, proximal, or ambiguous samples.
CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
distal_dist = CalculateEuclideanDistance(distal_ep_mean$V1, wcm_mat1$estimate)
prox_dist = CalculateEuclideanDistance(proximal_ep_mean$V1, wcm_mat1$estimate)
ambi_dist = CalculateEuclideanDistance(ambi_ep_mean$V1, wcm_mat1$estimate)
distances = c(distal_dist = CalculateEuclideanDistance(distal_ep_mean$V1, wcm_mat1$estimate),
              prox_dist = CalculateEuclideanDistance(proximal_ep_mean$V1, wcm_mat1$estimate),
              ambi_dist = CalculateEuclideanDistance(ambi_ep_mean$V1, wcm_mat1$estimate))

#Results are saved in the distances data.frame object (figure was manually made based on this object).
distances = as.data.frame(distances)
distances$Vector = rownames(distances)
DT::datatable(distances, options = list(scrollX = TRUE), caption = "Distances Table")
                     
# ------------------------------------------------------------------------------------------------
# EDF 6B 
# ------------------------------------------------------------------------------------------------
```{r, extfig10A, results = "asis", echo = TRUE, out.width= "100%", fig.align = "center"}
#Load libraries.
library(skitools)

tmp.plot = readRDS("../data/EDF/EDF6/6BC_data.rds")

highlight_color <- "darkslategray3"  # Color for WCM-1
grey_color <- "grey"        # Color for other patients

tmp.plot$color <- ifelse(tmp.plot$patient == "WCM-1", highlight_color, grey_color)

p <- ggplot(tmp.plot, aes(x, y, colour = color), size = 3) + 
  geom_point() + 
  theme_bw() + 
  theme(legend.position = "top") +
  scale_colour_identity()
print(p)
                     
# ------------------------------------------------------------------------------------------------
# EDF 6C
# ------------------------------------------------------------------------------------------------
#Load libraries.
library(skitools)
library(RColorBrewer)

tmp.plot = readRDS("../data/EDF/EDF6/6BC_data.rds")

my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(tmp.plot$cluster)))
my_colors[1] <- 'turquoise4'
my_colors[2] <- 'mediumseagreen'
my_colors[3] <- 'darkred'
names(my_colors) <- unique(tmp.plot$cluster)
tmp.plot$color <- ifelse(tmp.plot$patient == "WCM-1", my_colors[tmp.plot$cluster], "grey")
p <- ggplot(tmp.plot, aes(x = x, y = y, colour = factor(cluster)), size = 3) + 
  geom_point(aes(color = ifelse(patient == "WCM-1", as.character(cluster), "Other"))) + 
  scale_color_manual(values = c(my_colors, Other = "grey"), 
                     breaks = c(names(my_colors), "Other"),
                     labels = c(names(my_colors), "Other Patients")) +
  theme_bw() + 
  theme(legend.position = "top") 
print(p)
