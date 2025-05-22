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
library(ggsankey)
library(effects)
library(stats)
library(forcats) 
library(ggpubr)
library(ggforce)
library(ComplexHeatmap)

# ------------------------------------------------------------------------------------------------
# EDF 1A
# ------------------------------------------------------------------------------------------------

#Load the required libraries.
library(skitools)
library(Seurat)

#Load the expression matrix and the cluster data.
expDat = readRDS("../data/sup1AExp.rds")
clustDat = readRDS("../data/sup1AClusters.rds")

#Construct the seurat object.
seu = CreateSeuratObject(counts = expDat)
seu@meta.data$seurat_clusters = clustDat

#Generate the Dotplot.
p = DotPlot(seu, features = rownames(expDat), group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ppdf(print(p),filename="../data/Supp1A.pdf")


# ------------------------------------------------------------------------------------------------
# EDF 1B & 1C
# ------------------------------------------------------------------------------------------------

#Load required libraries.
library(skitools)
library(Seurat)
library(ggpubr)

#Load seurat metadata (with clustering results and cell names), load per cell aneploidy scores, and load shannon entropy scores per cluster.
tmp.plot = data.table(readRDS("~/projects/scLung/db/seuMetaData.rds"))
tmp.l.norm = readRDS("~/../spanja/Projects/COO/anneuploidy zscore.rds")
exp_shannon = readRDS("/gpfs/commons/home/jandrademartínez/projects/scLung/db/clusterEntropy.rds")

#Create table of total aneuploidy score per cell. Add to the metadata table.
aneu.score = data.table(cell = rownames(tmp.l.norm), Aneuplody.Score = apply(tmp.l.norm,1,function(x) sum(abs(x))))
tmp.plot = merge(tmp.plot,aneu.score,by = 'cell')

#Compute aneuploidy score per cluster. 
tmp.df = aggregate(Aneuplody.Score ~ seurat_clusters, data = tmp.plot, FUN = mean)
colnames(tmp.df)[2] = 'Aneuploidy.Score.mean'
#Add diversity (entropy) value per cluster. 
tmp.df$S.Diversity.num = exp_shannon


#Label carcinoma cluster based on high aneuploidy and low entropy.
tmp.df$cat = ifelse(((tmp.df$S.Diversity.num < 10) &(tmp.df$Aneuploidy.Score.mean > 40)),"Carcinoma","Others")
tmp.nw = tmp.plot[,.(Tissue,seurat_clusters)]
tmp.nw = as.data.frame(tmp.nw)
tmp.dt = tmp.df
tmp.df = as.data.frame(tmp.df)
tmp.df = merge(tmp.df,tmp.nw,by='seurat_clusters')
tmp.df.nw = tmp.df %>%  distinct(seurat_clusters, .keep_all = TRUE)
tmp.df.nw = as.data.table(tmp.df.nw)
#Add annotations for non-carcinoma clusters based on their associated tissue (tumor associated epithlial or normal).
tmp.df.nw[which((tmp.df.nw$Tissue == "Tumor" | tmp.df.nw$Tissue == "Metastasis") & (tmp.df.nw$cat == "Others")),"cat"] = "Tumor Associated Epithelial"
tmp.df.nw[which(tmp.df.nw$Tissue== "Adjacent Lung") & (tmp.df.nw$cat == "Others"),"cat"] = "Normal"
tmp.df.nw[13,'cat']="Carcinoma"

#Generate scatterplot for Sup 1B.
p = ggscatter(tmp.df.nw, x = "S.Diversity.num", y = "Aneuploidy.Score.mean",
              color = 'cat',
              repel = TRUE, xlab = "Shannon's Diversity", ylab = 'Aneuploidy Score (mean)')
ppdf(print(p), cex = 0.7,filename="../data/Supp1C.pdf")

#Link carcinoma and non-carcinoma annotations to each cell.
tmp.plot$cat = ""
for(i in 1:nrow(tmp.df.nw)){
    tmp.plot[tmp.plot$seurat_clusters==tmp.df.nw$seurat_clusters[i],]$cat = tmp.df.nw$cat[i]  
    }
#Organize cells by annotation.
tmp.plot$cat = factor(tmp.plot$cat, levels = c("Normal","Tumor Associated Epithelial","Carcinoma"))
tmp.plot = tmp.plot[order(tmp.plot$cat),]                   

#Organize rows in aneuploidy output for heatmap.
tmp.l.nw = tmp.l.norm[match(tmp.plot$cell,rownames(tmp.l.norm)),]

#Set heatmap annotation and color.
row.ha = HeatmapAnnotation(Tissue = tmp.plot$cat, name = 'Tissue', width = unit(2, "cm"), show_legend = TRUE, show_annotation_name = FALSE, col = list(Tissue = c("Carcinoma" = "#F8766D", "Normal" = "#00BA38", "Tumor Associated Epithelial" = "#619CFF")), annotation_legend_param = list("Normal", "Tumor Epithelial", "Carcinoma"), which = "row")
library(circlize)
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
#Plot heatmap.
p = Heatmap(tmp.l.nw, name = 'z-score', col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 8), column_names_side = c("bottom"), show_column_names = TRUE, show_row_names = FALSE, column_title_gp = gpar(fontsize = 20, fontface = "bold")) + row.ha
ppdf(print(p), cex = c(1,2),filename="../data/Supp1C.pdf")

# ------------------------------------------------------------------------------------------------
# EDF 2
# ------------------------------------------------------------------------------------------------

#Load requires libraries.
library(skitools)
library(data.tree)

#Load sikkema metadata table. Format data so that each row represents a full hierarchy for a given celltype.
sikkema_meta = read.table("../data/sikk_meta_epi.csv",header=T,sep=",")
sikkema_meta = sikkema_meta[,c("ann_level_1","ann_level_2","ann_level_3","ann_level_4","ann_finest_level")]

#Adjust AT0 annotation to conform with Alveolar annotation.
sikkema_meta[sikkema_meta$ann_finest_level=="AT0",]$ann_level_2 = "Alveolar epithelium"
sikkema_meta[sikkema_meta$ann_finest_level=="AT0",]$ann_level_3 = "AT2"
sikkema_meta[sikkema_meta$ann_finest_level=="AT0",]$ann_level_4 = "AT2"

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

# Render the left-to-right tree plot and export to pdf.
DiagrammeR::export_graph(
  tree_graph,
  file_name = "../data//hierarchical_dendrogram_v3.pdf",  # specify file name
  file_type = "pdf"
)
 

# ------------------------------------------------------------------------------------------------
# EDF 3
# ------------------------------------------------------------------------------------------------

#Load libraries.
library(skitools)

#Load file with cell metadata and identity call results.
normal_labs = fread("../data/combined_emb_for_normal_for_homogenous_sikk_ep_nw5.csv")
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
# EDF 4A
# ------------------------------------------------------------------------------------------------

#Load required libraries.
library(skitools)

#Load GRanges with centroid expression data.
juan.genes.gr.nw = readRDS("../data/genes_gr_LUAD.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)

#Define quantiles for partitioning the genes by AT2 expression.
AT0.quant <- quantile(juan.genes.gr.nw.dt$AT2)
AT0.quant.2 <- quantile(juan.genes.gr.nw.dt$AT2, probs = c(0.33, 0.67, 1))

#Define a data.frame starting from the expression information. 
sikk_AT0 <- juan.genes.gr.nw.dt$AT2 %>% as.data.frame()
#Add quantile information and gene name.
colnames(sikk_AT0)[1] = 'AT2'
sikk_AT0$gene <- juan.genes.gr.nw.dt$gene_name
lowAT0.quant <- sikk_AT0[which(sikk_AT0$AT2 <= -6.886128), ]
midAT0.quant <- sikk_AT0[which(sikk_AT0$AT2 > -6.886128  & sikk_AT0$AT2 <= -5.963367) , ]
highAT0.quant <- sikk_AT0[which(sikk_AT0$AT2 > -5.963367) , ]
lowAT0.quant$Quart = 'Q1'
midAT0.quant$Quart = 'Q2'
highAT0.quant$Quart = 'Q3'
AT0.quart.mat <- rbind(lowAT0.quant,midAT0.quant,highAT0.quant)

#Load GRanges with snv counts information per gene for LUAD. Convert to mutational density.
mgenes= readRDS("../data/mgenes.rds")
mgenes2 <- copy(mgenes)
juan.genes.gr.nw.dt$densityLUAD <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*246)
saveRDS(mgenes,"~/mgenes.rds")

#Make sure quantile data.frame gene order matches snv counts GRanges gene order. Remove NAs. Then add gene mutational density information for LUAD to the data.frame.
AT0LUAD_centwrtTMB = AT0.quart.mat[match(juan.genes.gr.nw.dt$gene_name,AT0.quart.mat$gene),]
AT0LUAD_centwrtTMB = AT0LUAD_centwrtTMB[complete.cases(AT0LUAD_centwrtTMB),]  
AT0.luadmat = cbind(AT0LUAD_centwrtTMB,juan.genes.gr.nw.dt$densityLUAD)
colnames(AT0.luadmat)[[4]] = "densityLUAD"

#Compute mean and standard deviation values for mutational density per expression quantile.
AT0.luadmat.Q1 = AT0.luadmat[which(AT0.luadmat$Quart=="Q1"),]
AT0.luadmat.Q1mean = mean(AT0.luadmat.Q1$densityLUAD)
AT0.luadmat.Q1se = sd(AT0.luadmat.Q1$densityLUAD)/sqrt(length(AT0.luadmat.Q1$densityLUAD)) 
AT0.luadmat.Q2 = AT0.luadmat[which(AT0.luadmat$Quart=="Q2"),]
AT0.luadmat.Q2mean = mean(AT0.luadmat.Q2$densityLUAD)
AT0.luadmat.Q2se = sd(AT0.luadmat.Q2$densityLUAD)/sqrt(length(AT0.luadmat.Q2$densityLUAD))
AT0.luadmat.Q3 = AT0.luadmat[which(AT0.luadmat$Quart=="Q3"),]
AT0.luadmat.Q3mean = mean(AT0.luadmat.Q3$densityLUAD)
AT0.luadmat.Q3se = sd(AT0.luadmat.Q3$densityLUAD)/sqrt(length(AT0.luadmat.Q3$densityLUAD))
#Save mean and sd values as a new data.frame.
mean_AT0_mat = as.data.frame(rbind(AT0.luadmat.Q1mean,AT0.luadmat.Q2mean,AT0.luadmat.Q3mean))
se_AT0_mat = as.data.frame(rbind(AT0.luadmat.Q1se,AT0.luadmat.Q2se,AT0.luadmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
AT0.mean.mat = cbind(quart.mat,mean_AT0_mat,se_AT0_mat)
rownames(AT0.mean.mat) = c(1,2,3)
colnames(AT0.mean.mat) = c("Quart","Mean_TMB","SE")

#Now repeat the complete process but using the  basal-resting centroid expression data.
#Compute expression quantiles.
suprabasal.quant <- quantile(juan.genes.gr.nw.dt$Basal.resting)
suprabasal.quant.2 <- quantile(juan.genes.gr.nw.dt$Basal.resting, probs = c(0.33, 0.67, 1))
sikk_sup <- juan.genes.gr.nw.dt$Basal.resting %>% as.data.frame()
colnames(sikk_sup)[1] = 'goblet'
sikk_sup$gene <- juan.genes.gr.nw.dt$gene_name

#Assign quantile labels to each gene.
lowsup.quant <- sikk_sup[which(sikk_sup$goblet <= -6.888385), ]
midsup.quant <- sikk_sup[which(sikk_sup$goblet > -6.888385  & sikk_sup$goblet <= -6.101867) , ]
highsup.quant <- sikk_sup[which(sikk_sup$goblet > -6.101867) , ]
lowsup.quant$Quart = 'Q1'
midsup.quant$Quart = 'Q2'
highsup.quant$Quart = 'Q3'
sup.quart.mat <- rbind(lowsup.quant,midsup.quant,highsup.quant)

#Make a data.frame with expression quantile and mutational density. Compute mean and standard deviation of mutational density per expression quantile. 
supLUAD_centwrtTMB = sup.quart.mat[match(juan.genes.gr.nw.dt$gene_name,sup.quart.mat$gene),]
supLUAD_centwrtTMB = supLUAD_centwrtTMB[complete.cases(supLUAD_centwrtTMB),]
sup.luadmat = cbind(supLUAD_centwrtTMB,juan.genes.gr.nw.dt$densityLUAD)
colnames(sup.luadmat)[[4]] = "densityLUAD"
sup.luadmat.Q1 = sup.luadmat[which(sup.luadmat$Quart=="Q1"),]
sup.luadmat.Q1mean = mean(sup.luadmat.Q1$densityLUAD)
sup.luadmat.Q1se = sd(sup.luadmat.Q1$densityLUAD)/sqrt(length(sup.luadmat.Q1$densityLUAD)) 
sup.luadmat.Q2 = sup.luadmat[which(sup.luadmat$Quart=="Q2"),]
sup.luadmat.Q2mean = mean(sup.luadmat.Q2$densityLUAD)
sup.luadmat.Q2se = sd(sup.luadmat.Q2$densityLUAD)/sqrt(length(sup.luadmat.Q2$densityLUAD))
sup.luadmat.Q3 = sup.luadmat[which(sup.luadmat$Quart=="Q3"),]
sup.luadmat.Q3mean = mean(sup.luadmat.Q3$densityLUAD)
sup.luadmat.Q3se = sd(sup.luadmat.Q3$densityLUAD)/sqrt(length(sup.luadmat.Q3$densityLUAD))
mean_sup_mat = as.data.frame(rbind(sup.luadmat.Q1mean,sup.luadmat.Q2mean,sup.luadmat.Q3mean))
se_sup_mat = as.data.frame(rbind(sup.luadmat.Q1se,sup.luadmat.Q2se,sup.luadmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

#Combine data.frames with mean and standard deviation per expression quartile for AT2 and basal resting cells. 
quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
sup.mean.mat = cbind(quart.mat,mean_sup_mat,se_sup_mat)
rownames(sup.mean.mat) = c(1,2,3)
colnames(sup.mean.mat) = c("Quart","Mean_TMB","SE")
sup.mean.mat$celltype = "basalresting"
AT0.mean.mat$celltype="AT2"
mean.mat.luad = rbind(AT0.mean.mat,sup.mean.mat)
mean.mat.luad$celltype = as.factor(mean.mat.luad$celltype)

#Plot the graph.
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
  )

print(p)
                                                                             
# ------------------------------------------------------------------------------------------------
# EDF 4B
# ------------------------------------------------------------------------------------------------

#Load required libraries.
library(skitools)

#Load GRanges with centroid expression data.
juan.genes.gr.nw = readRDS("../data/genes_gr_LUSC.rds")
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
mgenes= readRDS("../data/mgenes.rds")
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
# EDF 6A
# ------------------------------------------------------------------------------------------------

# Accuracy for predicting true cell type at different levels 

sims_acc_results <- fread('../data/edf6_luad_sims_results.csv')

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
# EDF 7
# ------------------------------------------------------------------------------------------------

## EDF 7 --Accuracy Plots ##
#Load libraries.
library(skitools)

#Plots can be made for univariate and multivariate simultaneously or individually. Comments are made when relevant to do one or the other.

#Load simulation result files. Add Method column to identify if results are uni or multivariate.
resDT = readRDS("../data/edf7_testRun_withOverlap_Cov.rds")
resDT$Method = "Univariate"

#If results are plotted in the same graph, uncomment the following lines, adjusting the Method column value accordingly. 
#MulDT = readRDS("~/projects/scLung/db/edf7_testRun_withOverlapMultivariatewithCovs.rds")
# MulDT$Method = "Multivariate"
#resDT = rbind(resDT,MulDT,fill=TRUE)

#Load csv with Sikkema annotations.
annotSikk = read.csv("../data/edf7_sikk_meta_epi3.csv")
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
    if(i%%100==0){
        print(i)
     }
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
average_df = resDT %>% group_by(TMB) %>% summarise(average_value = mean(Lev_Hit, na.rm = TRUE))  
average_df = as.data.frame(average_df)  

#Line graph is plotted. Depending on whether both uni and multivariate are plotted use the first line or the second (comment/uncomment adequately).
p =  ggplot(data=average_df, aes(x=log10(TMB), y=average_value)) + geom_line(color="blue")+ geom_point() + ggtitle("Accuracy per TMB Multivariate - Proximal vs Distal" ) + xlab("Log10(TMB) (mut/mb)") + ylab("Accuracy") + scale_y_continuous(limits=c(0,1)) +   theme(text = element_text(size = 20)) 
#p =  ggplot(data=average_df, aes(x=log10(TMB), y=average_value, color = Method, group=Method)) + geom_line(color="red")+ geom_point(aes(shape=Method),size=2) + ggtitle("Accuracy per TMB - Celltype" ) + xlab("Log10(TMB) (mut/mb)") + ylab("Accuracy") + scale_y_continuous(limits=c(0,1)) +   theme(text = element_text(size = 20)) +  scale_color_manual(values = c("Univariate" = "red", "Multivariate" = "blue")) 
ppdf(print(p),file="edf7.pdf")





# ------------------------------------------------------------------------------------------------
# EDF 8A
# ------------------------------------------------------------------------------------------------
# Heat map for LUSC 

nsclc_col <- readRDS('../data/nsclc_hm_col.rds')
hm_lusc <- fread('../data/hm_lusc.csv')
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
# EDF 8B
# ------------------------------------------------------------------------------------------------

# LUSC - Cluster 1 Distal Lung
edf_8_data <- readRDS('../data/edf8.rds')
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
# EDF 9A
# ------------------------------------------------------------------------------------------------



pts_coo_id_edf9a <- readRDS( '../data/pts_coo_id_4B.rds')
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
# EDF 9B
# ------------------------------------------------------------------------------------------------
tmp.plot.luad.df3 <- readRDS('../data/edf9b.rds')

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



# ------------------------------------------------------------------------------------------------
# EDF 9G
# ------------------------------------------------------------------------------------------------
library(skitools)

edf9g_a <- readRDS('../data/edf9g_a.rds')
edf9g_b <- readRDS('../data/edf9g_b.rds')

oncoprint(edf9g_b, oncotab = edf9g_a,  genes = c('NKX2-1','SMARCA4','STK11','APC','KEAP1','ALK','MCL1','MAP2K1','EGFR','FGFR1','FOXP1','CDK6','BCL2L1', 'TP53','RB1','CDKN2A','TERT','MYC','ARID1A','PTEN','CCND1','PIK3CA','ERBB2','CCNE1','NF1'), sort.genes = FALSE, colnames.fontsize = 20, 
          rownames.fontsize = 15, signature.thresh = 20, track.height = 1.5, split.gap = 0.5, wes = TRUE, track.gap = track.height/2, drop = TRUE, drop.genes = TRUE, svevents = FALSE,  sort.tumors = FALSE,
          return.oncotab = FALSE, height = 12, width = 25,  filename = 'EDF9g_onco.pdf')

# ------------------------------------------------------------------------------------------------
# EDF 10 A 
# ------------------------------------------------------------------------------------------------
tmp.plot = readRDS("../data/edf10.rds")

highlight_color <- "darkslategray3"  # Color for WCM-1
grey_color <- "grey"        # Color for other patients

tmp.plot$color <- ifelse(tmp.plot$patient == "WCM-1", highlight_color, grey_color)

p <- ggplot(tmp.plot, aes(x, y, colour = color), size = 3) + 
  geom_point() + 
  theme_bw() + 
  theme(legend.position = "top") +
  scale_colour_identity()
ppdf(print(p),filename="UMAP with WCM-1 patient knn6_nw1.pdf")

# ------------------------------------------------------------------------------------------------
# EDF 10 B
# ------------------------------------------------------------------------------------------------
library(RColorBrewer)

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
ppdf(print(p),filename="UMAP with WCM-1 with clusters in knn6.pdf")
