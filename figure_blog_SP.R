## Figure 1C ###
#Load libraries.
library(skitools)
library(dplyr)

#Load identity calls per cell and patient.
combined_emb = read.table("~/../spanja/Projects/COO/combined_emb_for_carc_for_homogenous_sikk_ep_500_nw.csv",sep=",",header=TRUE)
query_samp = combined_emb[which(combined_emb$ref_or_query=="query"),]
cellcount_mat = query_samp[,c("X","Patient")]
aa = as.data.frame(table(cellcount_mat$Patient))
grtr_than_50 = aa[which(aa$Freq>50),]
query_samp_nw = query_samp[which(query_samp$Patient %in% grtr_than_50$Var1),]
query_samp_nw$ann_finest_lev_transferred_label_filtered = ifelse(query_samp_nw$ann_finest_lev_transfer_uncert<0.5,query_samp_nw$ann_finest_lev_transferred_label_unfiltered,"Unknown")
query_samp_nw = query_samp_nw[which(query_samp_nw$ann_finest_lev_transferred_label_filtered!="Unknown"),]

                                        #Select patients to plot (LUAD is the variable name but it includes LUSC and LCNEC patients as well).
LUAD = c('LX681','UHL-7','WCM-3','UHL-3','LX666','WCM-2','LX701','LX679','WCM-4','UHL-5','UHL-8','LX661','LX699','LX676','UHL-1','UHL-4','UHL-6','WCM-1','PS05','PS06','PS09','PS03','PS04','PS10','UHL-2','LX680')
LUAD_query = query_samp_nw[which(query_samp_nw$Patient %in% LUAD),]
distal = c("AT0","AT1","AT2","AT2 proliferating")
rare = c("Ionocyte","Tuft","Neuroendocrine")
distal_rare = c(distal,rare)
epicells = c("AT1","AT2","AT2 proliferating","AT0","Suprabasal","Basal resting","Hillock-like","Multiciliated (non-nasal)","Multiciliated (nasal)","Deuterosomal","Neuroendocrine","Ionocyte","Tuft","Goblet (nasal)","Club (nasal)","Club (non-nasal)","pre-TB secretory","Goblet (bronchial)","Goblet (subsegmental)","SMG serous (bronchial)","SMG mucous","SMG duct","SMG serous (nasal)")
proximal = setdiff(epicells,distal_rare)
LUAD_query$cat = ifelse(LUAD_query$ann_finest_lev_transferred_label_filtered %in% distal,"distal",LUAD_query$ann_finest_lev_transferred_label_filtered)

#Group cells in categories for plotting.
LUAD_query$cat = ifelse(LUAD_query$cat=="rare",LUAD_query$ann_finest_lev_transferred_label_filtered,LUAD_query$cat)
AT1 = c("AT1")
AT2 = c("AT2","AT2 proliferating","AT0")
basal = c("Suprabasal","Basal resting","Hillock-like")
multicilated = c("Multiciliated (non-nasal)","Multiciliated (nasal)","Deuterosomal")
rare = c("Neuroendocrine","Ionocyte","Tuft")
secretory = c("Goblet (nasal)","Club (non-nasal)","Club (nasal)","pre-TB secretory","Goblet (bronchial)","Goblet (subsegmental)")
smg = c("SMG serous (bronchial)","SMG mucous","SMG duct","SMG serous (nasal)")
LUAD_query[LUAD_query$ann_finest_lev_transferred_label_filtered %in% AT1,"cat1"]="AT1"
LUAD_query[LUAD_query$ann_finest_lev_transferred_label_filtered %in% AT2,"cat1"]="AT2"
LUAD_query[LUAD_query$ann_finest_lev_transferred_label_filtered %in% basal,"cat1"]="Basal"
LUAD_query[LUAD_query$ann_finest_lev_transferred_label_filtered %in% multicilated,"cat1"]="multiciliated"
LUAD_query[LUAD_query$ann_finest_lev_transferred_label_filtered %in% secretory,"cat1"]="secretory"
LUAD_query[LUAD_query$ann_finest_lev_transferred_label_filtered %in% smg,"cat1"]="smg"
LUAD_query$cat1 = ifelse(LUAD_query$ann_finest_lev_transferred_label_filtered %in% rare,LUAD_query$cat,LUAD_query$cat1)

#Consolidate results into percentages per cell group per patient. 
dft2 <- data.frame(LUAD_query) %>% group_by(Patient, cat1) %>% summarise(count = n()) %>%             # Count the rows per group 
  mutate(total = sum(count),             # Calculate total counts per patient
         percentage = (count / total) * 100) %>%  # Calculate percentage
  dplyr::select(-total)

#Set patient order and plot.
dft2$Patient <- factor(dft2$Patient, levels = c('LX681','UHL-7','WCM-3','UHL-3','LX666','WCM-2','LX701','LX679','WCM-4','UHL-5','UHL-8','LX661','LX699','LX676','UHL-1','UHL-4','UHL-6','WCM-1','PS05','PS06','PS09','PS03','PS04','PS10','UHL-2','LX680'))

p = ggplot(data=dft2, aes(x=Patient, y=percentage, fill=cat1)) + 
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("AT1"= '#1F77B4',"AT2"="#FF7F0E","Basal"="#279E68",
                             "multiciliated"="#D62728","Neuroendocrine"="#AA0","Tuft" = "salmon","Ionocyte"="grey",
                             "secretory"="#8C564B","smg"="#E377C2"))
ppdf(print(p),filename="Fig1C.pdf")


## SUpplementary Fig 1 A

seu = readRDS("~/../spanja/Projects/COO/seu_data_zhang_and_kofi.rds")

seu = FindNeighbors(seu, k.param = 15, reduction = "pca", dims = 1:50, verbose = TRUE)
seu = FindClusters(seu, resolution = 0.5, group.singletons = TRUE, verbose = TRUE)

epGenes=c("EPCAM","AGER","CLDN18","SFTPA1","SFTPC","KRT5","KRT17","FOXJ1","FHAD1","LYPD2","CLDN10","TFF1","TFF2","ASCL1","CALCA")

expression_data = GetAssayData(seu, assay = "RNA")[epGenes,]
saveRDS(expression_data,"/gpfs/commons/home/jandrademartínez/projects/scLung/db/sup1AExp.rds")
saveRDS(seu@meta.data$seurat_clusters,"/gpfs/commons/home/jandrademartínez/projects/scLung/db/sup1AClusters.rds")

#Load the required libraries.
library(skitools)
library(Seurat)

#Load the expression matrix and the cluster data.
expDat = readRDS("/gpfs/commons/home/jandrademartínez/projects/scLung/db/sup1AExp.rds")
clustDat = readRDS("/gpfs/commons/home/jandrademartínez/projects/scLung/db/sup1AClusters.rds")

#Construct the seurat object.
seu = CreateSeuratObject(counts = expDat)
seu@meta.data$seurat_clusters = clustDat

#Generate the Dotplot.
p = DotPlot(seu, features = rownames(expDat), group.by="seurat_clusters") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ppdf(print(p),filename="Supp1A.pdf")

## Supplementary Fig 1B & 1C
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
ppdf(print(p), cex = 0.7,filename="shannon diversity vs aneuploidy score seurat clusters_nw1b.pdf")

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
ppdf(print(p), cex = c(1,2),filename="heatmap_for_anneuploidy_v3.pdf")

## Supplementary Fig 2

#Load requires libraries.
library(skitools)
library(data.tree)

#Load sikkema metadata table. Format data so that each row represents a full hierarchy for a given celltype.
sikkema_meta = read.table("/gpfs/commons/groups/imielinski_lab/projects/scLung/db/sikk_meta_epi.csv",header=T,sep=",")
sikkema_meta = sikkema_meta[,c("ann_level_1","ann_level_2","ann_level_3","ann_level_4","ann_finest_level")]
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

DiagrammeR::export_graph(
  tree_graph,
  file_name = "/gpfs/commons/groups/imielinski_lab/projects/scLung/db/hierarchical_dendrogram_v3.pdf",  # specify file name
  file_type = "pdf"
)
 
## Supplementary Fig 3

#Load libraries.
library(skitools)

#Load file with cell metadata and identity call results.
normal_labs = fread("/gpfs/commons/groups/imielinski_lab/projects/scLung/db/scLung_RMDFiles/combined_emb_for_normal_for_homogenous_sikk_ep_nw5.csv")
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


## Supplementary Fig 4 A
juan.genes.gr.nw = readRDS("~/Projects/COO/genes_gr_LUAD.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)

AT0.quant <- quantile(juan.genes.gr.nw.dt$AT2)
AT0.quant.2 <- quantile(juan.genes.gr.nw.dt$AT2, probs = c(0.33, 0.67, 1))
sikk_AT0 <- juan.genes.gr.nw.dt$AT2 %>% as.data.frame()
colnames(sikk_AT0)[1] = 'AT2'
sikk_AT0$gene <- juan.genes.gr.nw.dt$gene_name
lowAT0.quant <- sikk_AT0[which(sikk_AT0$AT2 <= -6.886128), ]
midAT0.quant <- sikk_AT0[which(sikk_AT0$AT2 > -6.886128  & sikk_AT0$AT2 <= -5.963367) , ]
highAT0.quant <- sikk_AT0[which(sikk_AT0$AT2 > -5.963367) , ]
dim(lowAT0.quant )
lowAT0.quant$Quart = 'Q1'
dim(midAT0.quant )
midAT0.quant$Quart = 'Q2'
dim(highAT0.quant )
highAT0.quant$Quart = 'Q3'
AT0.quart.mat <- rbind(lowAT0.quant,midAT0.quant,highAT0.quant)
head(AT0.quart.mat)
mgenes= readRDS("~/data/scLung/mgenes.rds")
mgenes2 <- copy(mgenes)
juan.genes.gr.nw.dt$densityLUAD <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*246)
#juan_genes_upd$densityLUSC <- (10**6)*juan_genes_upd$snv.count.lusc_num_54/(juan_genes_upd$width*54)
#juan_genes_upd$densitySCLC <- (10**6)*juan_genes_upd$snv.count.hmf.sclc_num_54/(juan_genes_upd$width*54)

AT0LUAD_centwrtTMB = AT0.quart.mat[match(juan.genes.gr.nw.dt$gene_name,AT0.quart.mat$gene),]
AT0LUAD_centwrtTMB = AT0LUAD_centwrtTMB[complete.cases(AT0LUAD_centwrtTMB),]  # remove NAs
#AT2TMB = mgenes2[match(AT2LUAD_centwrtTMB$gene,mgenes2$gene_name),]
AT0.luadmat = cbind(AT0LUAD_centwrtTMB,juan.genes.gr.nw.dt$densityLUAD)
colnames(AT0.luadmat)[[4]] = "densityLUAD"
AT0.luadmat.Q1 = AT0.luadmat[which(AT0.luadmat$Quart=="Q1"),]
AT0.luadmat.Q1mean = mean(AT0.luadmat.Q1$densityLUAD)
AT0.luadmat.Q1se = sd(AT0.luadmat.Q1$densityLUAD)/sqrt(length(AT0.luadmat.Q1$densityLUAD)) 
AT0.luadmat.Q2 = AT0.luadmat[which(AT0.luadmat$Quart=="Q2"),]
AT0.luadmat.Q2mean = mean(AT0.luadmat.Q2$densityLUAD)
AT0.luadmat.Q2se = sd(AT0.luadmat.Q2$densityLUAD)/sqrt(length(AT0.luadmat.Q2$densityLUAD))
AT0.luadmat.Q3 = AT0.luadmat[which(AT0.luadmat$Quart=="Q3"),]
AT0.luadmat.Q3mean = mean(AT0.luadmat.Q3$densityLUAD)
AT0.luadmat.Q3se = sd(AT0.luadmat.Q3$densityLUAD)/sqrt(length(AT0.luadmat.Q3$densityLUAD))
mean_AT0_mat = as.data.frame(rbind(AT0.luadmat.Q1mean,AT0.luadmat.Q2mean,AT0.luadmat.Q3mean))
se_AT0_mat = as.data.frame(rbind(AT0.luadmat.Q1se,AT0.luadmat.Q2se,AT0.luadmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
AT0.mean.mat = cbind(quart.mat,mean_AT0_mat,se_AT0_mat)
rownames(AT0.mean.mat) = c(1,2,3)
colnames(AT0.mean.mat) = c("Quart","Mean_TMB","SE")

suprabasal.quant <- quantile(juan.genes.gr.nw.dt$Basal.resting)

suprabasal.quant.2 <- quantile(juan.genes.gr.nw.dt$Basal.resting, probs = c(0.33, 0.67, 1))
sikk_sup <- juan.genes.gr.nw.dt$Basal.resting %>% as.data.frame()
colnames(sikk_sup)[1] = 'goblet'
sikk_sup$gene <- juan.genes.gr.nw.dt$gene_name
suprabasal.quant
# 33%          67%         100%
# 1.113533e-05 5.937357e-04 6.213304e+00
# lowAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 <= 3.647752e-06), ]
# midAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > 3.647752e-06 & sikk_AT2$AT2 <= 1.082263e-03) , ]
# highAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > 1.082263e-03) , ]
lowsup.quant <- sikk_sup[which(sikk_sup$goblet <= -6.888385), ]
midsup.quant <- sikk_sup[which(sikk_sup$goblet > -6.888385  & sikk_sup$goblet <= -6.101867) , ]
highsup.quant <- sikk_sup[which(sikk_sup$goblet > -6.101867) , ]
dim(lowsup.quant )
lowsup.quant$Quart = 'Q1'
dim(midsup.quant )
midsup.quant$Quart = 'Q2'
dim(highsup.quant )
highsup.quant$Quart = 'Q3'
sup.quart.mat <- rbind(lowsup.quant,midsup.quant,highsup.quant)
#head(pretb.quart.mat)
mgenes= readRDS("~/data/scLung/mgenes.rds")
mgenes2 <- copy(mgenes)
juan.genes.gr.nw.dt$densityLUAD <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*246)
#juan_genes_upd$densityLUSC <- (10**6)*juan_genes_upd$snv.count.lusc_num_54/(juan_genes_upd$width*54)
#juan_genes_upd$densitySCLC <- (10**6)*juan_genes_upd$snv.count.hmf.sclc_num_54/(juan_genes_upd$width*54)



#saveRDS(mgenes2, '~/data/scLung/mgenes2.rds')
#AT2.quart.mat %>% class
supLUAD_centwrtTMB = sup.quart.mat[match(juan.genes.gr.nw.dt$gene_name,sup.quart.mat$gene),]
supLUAD_centwrtTMB = supLUAD_centwrtTMB[complete.cases(supLUAD_centwrtTMB),]  # remove NAs
#AT2TMB = mgenes2[match(AT2LUAD_centwrtTMB$gene,mgenes2$gene_name),]
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

quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
sup.mean.mat = cbind(quart.mat,mean_sup_mat,se_sup_mat)
rownames(sup.mean.mat) = c(1,2,3)
colnames(sup.mean.mat) = c("Quart","Mean_TMB","SE")
sup.mean.mat$celltype = "basalresting"
AT0.mean.mat$celltype="AT2"

mean.mat.luad = rbind(AT0.mean.mat,sup.mean.mat)
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
  )

ppdf(print(p),filename="basalresting_AT2_TCNER_LUAD_fig2_line_dot_v2.pdf",width=10,height=10)

## Supplementary Fig 4B

juan.genes.gr.nw = readRDS("~/Projects/COO/genes_gr_LUSC.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)
sup.quant <- quantile(juan.genes.gr.nw.dt$Basal.resting)
sup.quant.2 <- quantile(juan.genes.gr.nw.dt$Basal.resting, probs = c(0.33, 0.67, 1))
sikk_sup <- juan.genes.gr.nw.dt$Basal.resting %>% as.data.frame()
colnames(sikk_sup)[1] = 'Goblet'
sikk_sup$gene <- juan.genes.gr.nw.dt$gene_name
sup.quant
lowsup.quant <- sikk_sup[which(sikk_sup$Goblet <= -6.888385), ]
midsup.quant <- sikk_sup[which(sikk_sup$Goblet > -6.888385  & sikk_sup$Goblet <= -6.101867) , ]
highsup.quant <- sikk_sup[which(sikk_sup$Goblet > -6.101867) , ]
dim(lowsup.quant )
lowsup.quant$Quart = 'Q1'
dim(midsup.quant )
midsup.quant$Quart = 'Q2'
dim(highsup.quant )
highsup.quant$Quart = 'Q3'
sup.quart.mat <- rbind(lowsup.quant,midsup.quant,highsup.quant)
head(sup.quart.mat)
mgenes= readRDS("~/data/scLung/mgenes.rds")
mgenes2 <- copy(mgenes)
#juan.genes.gr.nw.dt$densityLUAD <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*246)
juan.genes.gr.nw.dt$densityLUSC <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*53)
#juan_genes_upd$densitySCLC <- (10**6)*juan_genes_upd$snv.count.hmf.sclc_num_54/(juan_genes_upd$width*54)

#saveRDS(mgenes2, '~/data/scLung/mgenes2.rds')
#AT2.quart.mat %>% class
supLUSC_centwrtTMB = sup.quart.mat[match(juan.genes.gr.nw.dt$gene_name,sup.quart.mat$gene),]
supLUSC_centwrtTMB = supLUSC_centwrtTMB[complete.cases(supLUSC_centwrtTMB),]  # remove NAs
#AT2TMB = mgenes2[match(AT2LUAD_centwrtTMB$gene,mgenes2$gene_name),]
sup.luscmat = cbind(supLUSC_centwrtTMB,juan.genes.gr.nw.dt$densityLUSC)
colnames(sup.luscmat)[[4]] = "densityLUSC"
sup.luscmat.Q1 = sup.luscmat[which(sup.luscmat$Quart=="Q1"),]
sup.luscmat.Q1mean = mean(sup.luscmat.Q1$densityLUSC)
sup.luscmat.Q1se = sd(sup.luscmat.Q1$densityLUSC)/sqrt(length(sup.luscmat.Q1$densityLUSC)) 
sup.luscmat.Q2 = sup.luscmat[which(sup.luscmat$Quart=="Q2"),]
sup.luscmat.Q2mean = mean(sup.luscmat.Q2$densityLUSC)
sup.luscmat.Q2se = sd(sup.luscmat.Q2$densityLUSC)/sqrt(length(sup.luscmat.Q2$densityLUSC))
sup.luscmat.Q3 = sup.luscmat[which(sup.luscmat$Quart=="Q3"),]
sup.luscmat.Q3mean = mean(sup.luscmat.Q3$densityLUSC)
sup.luscmat.Q3se = sd(sup.luscmat.Q3$densityLUSC)/sqrt(length(sup.luscmat.Q3$densityLUSC))
mean_sup_mat_lusc = as.data.frame(rbind(sup.luscmat.Q1mean,sup.luscmat.Q2mean,sup.luscmat.Q3mean))
se_sup_mat_lusc = as.data.frame(rbind(sup.luscmat.Q1se,sup.luscmat.Q2se,sup.luscmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
sup.mean.mat = cbind(quart.mat,mean_sup_mat_lusc,se_sup_mat_lusc)
rownames(sup.mean.mat) = c(1,2,3)
colnames(sup.mean.mat) = c("Quart","Mean_TMB","SE")

AT2.quant <- quantile(juan.genes.gr.nw.dt$AT2)
AT2.quant.2 <- quantile(juan.genes.gr.nw.dt$AT2, probs = c(0.33, 0.67, 1))
sikk_AT2 <- juan.genes.gr.nw.dt$AT2 %>% as.data.frame()
colnames(sikk_AT2)[1] = 'AT2'
sikk_AT2$gene <- juan.genes.gr.nw.dt$gene_name
AT2.quant

lowAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 <= -6.886128), ]
midAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > -6.886128  & sikk_AT2$AT2 <= -5.963367) , ]
highAT2.quant <- sikk_AT2[which(sikk_AT2$AT2 > -5.963367) , ]
dim(lowAT2.quant )
lowAT2.quant$Quart = 'Q1'
dim(midAT2.quant )
midAT2.quant$Quart = 'Q2'
dim(highAT2.quant )
highAT2.quant$Quart = 'Q3'
AT2.quart.mat <- rbind(lowAT2.quant,midAT2.quant,highAT2.quant)
#head(pretb.quart.mat)
mgenes= readRDS("~/data/scLung/mgenes.rds")
mgenes2 <- copy(mgenes)
#juan.genes.gr.nw.dt$densityLUAD <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*246)
juan.genes.gr.nw.dt$densityLUSC <- (10**6)*juan.genes.gr.nw.dt$snv.count/(juan.genes.gr.nw.dt$width*53)
#juan_genes_upd$densitySCLC <- (10**6)*juan_genes_upd$snv.count.hmf.sclc_num_54/(juan_genes_upd$width*54)



#saveRDS(mgenes2, '~/data/scLung/mgenes2.rds')
#AT2.quart.mat %>% class
AT2lusc_centwrtTMB = AT2.quart.mat[match(juan.genes.gr.nw.dt$gene_name,AT2.quart.mat$gene),]
AT2lusc_centwrtTMB = AT2lusc_centwrtTMB[complete.cases(AT2lusc_centwrtTMB),]  # remove NAs
#AT2TMB = mgenes2[match(AT2LUAD_centwrtTMB$gene,mgenes2$gene_name),]
AT2.luscmat = cbind(AT2lusc_centwrtTMB,juan.genes.gr.nw.dt$densityLUSC)
colnames(AT2.luscmat)[[4]] = "densityLUSC"
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
se_AT2_mat = as.data.frame(rbind(AT2.luscmat.Q1se,AT2.luscmat.Q2se,AT2.luscmat.Q3se))
Q1.quart = "Q1"
Q2.quart = "Q2"
Q3.quart = "Q3"

quart.mat = as.data.frame(rbind(Q1.quart,Q2.quart,Q3.quart))
AT2.mean.mat = cbind(quart.mat,mean_AT2_mat,se_AT2_mat)
rownames(AT2.mean.mat) = c(1,2,3)
colnames(AT2.mean.mat) = c("Quart","Mean_TMB","SE")
AT2.mean.mat$celltype = "AT2"
sup.mean.mat$celltype="Basalresting"

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

ppdf(print(p), filename = "Basalresting_AT2_TCNER_LUSC_point_line_plot_v4.pdf",width=10,height=10)

## Fig 2 A

library(skitools)
library(MASS)
library(parallel)
juan.genes.gr = readRDS("~/data/scLung/Sikkema/infraobjs/DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExpr.rds")
count.int = grep("snv.non.apobec.count", colnames(mcols(juan.genes.gr)), value = TRUE)
count.int = as.data.frame(count.int)
count.int$pair = gsub("snv.non.apobec.count.","",count.int$count.int)
count.int$pair = toupper(count.int$pair)
luad_pairs = readRDS("~/data/scLung/Sikkema/pairs/new.pairs.luad.rds")
juan.genes.dt = gr2dt(juan.genes.gr)
juan.genes.dt = as.data.frame(juan.genes.dt)
epicells = c(

LUAD_patnt_mat = intersect(luad_pairs,count.int$pair)
LUAD_patnt_mat = as.matrix(LUAD_patnt_mat)
LUAD_patnt_mat = LUAD_patnt_mat[complete.cases(LUAD_patnt_mat),]
count.int.mat = count.int[match(LUAD_patnt_mat,count.int$pair),]
juan.genes.dt_nw = juan.genes.dt[,match(count.int.mat$count.int,colnames(juan.genes.dt))]
snv_count = as.matrix(rowSums(juan.genes.dt_nw))
juan_dt1 = juan.genes.dt[,1:6]
juan_dt1$snv.count = snv_count
sikk.cent = readRDS("~/Projects/COO/sikk_cent.rds")
sikk.cent_nw = sikk.cent[match(juan_dt1$gene_name,rownames(sikk.cent)),]
juan_dt2 = cbind(juan_dt1,sikk.cent_nw,juan.genes.dt[,6558:6562])
juan_dt2 = as.data.table(juan_dt2)
juan.genes.gr.nw = dt2gr(juan_dt2)


centroids= epicells
juan.genes.gr.nw$avgexpr = apply(mcols(juan.genes.gr.nw)[, colnames(mcols(juan.genes.gr.nw)) %in% centroids],1,mean)
juan.genes.dt1 = gr2dt(juan.genes.gr.nw)
juan.genes.dt1 = juan.genes.dt1[!is.na(gc),]
juan.genes.dt1 = juan.genes.dt1[!is.na(rept),]
#juan.genes.dt = dt2gr(juan.genes.dt)

## include only lung cell types for individual samples
cmut = grep('snv.count', names(mcols(juan.genes.gr.nw)), value = TRUE)

cts = epicells
combos = expand.grid(cmut, cts) %>% as.data.table


res = mclapply(1:nrow(combos), function(i)
{
  tt = combos[i,]$Var1 %>% as.character;
  celltype = combos[i,]$Var2 %>% as.character
  message(paste(tt, celltype));
  dat = data.table(snv = mcols(juan.genes.gr.nw)[[tt]],
                   gc = juan.genes.gr.nw$gc,
                   rept = juan.genes.gr.nw$rept,
                   avg = juan.genes.gr.nw$avgexpr,
                   cov =  mcols(juan.genes.gr.nw)[[celltype]],
                   exon.frac = log(juan.genes.gr.nw$exon.frac),
                   wid = log(width(juan.genes.gr.nw)))
  if(sum(dat$snv) < 1)
  {
    message(" -- No mutation data -- ")
    res = data.frame(name = c(NA,NA), method = c(NA,NA), p = c(1,1), estimate = c(1,1), ci.lower = c(1,1), ci.upper = c(1,1), effect = c(1, 1))
  } else {
    #res = glm(snv ~ cov + exon.frac + avg + rept + gc + offset(wid), family ="poisson",data = dat) %>% dflm
    res = glm.nb(snv ~ cov + avg + exon.frac + rept + gc + offset(wid), data = dat) %>% dflm
  }
  res$tt = tt
  res$celltype = celltype
  return(res)
}, mc.cores = 35, mc.preschedule = TRUE)

res = rbindlist(res)
res[, fdr := p.adjust(p), by = name] ## include BH correction
#res = readRDS("~/data/scLung/Sikkema/resforfirsquartile.rds")
#res$tt = paste0(res$tt,".LUSC")


## create grouped naming annotation
res$plot.label = "LUAD"
res$plot.label = gsub("snv.count.","",res$plot.label)
#res$plot.label = paste0(toupper(gsub("_num_"," (",res$plot.label)),')')

res$combination = paste0(toupper(gsub("_num_"," (",res$tt)), ')')
res$combination = gsub(")","",res$combination)
res$combination = paste0(res$combination," - ", res$celltype)
#res$tt = gsub("_num.*","",res$tt)
unique(res$plot.label)


## set estimates value for each celltype below fdr threshold to null value of 1
dt = dcast.data.table(res[grep('cov', name), ][, estimate2 := ifelse(fdr<0.1, estimate, 1)], tt ~ celltype, value.var = 'estimate2')

dt.plot = melt(dt, id.vars = c("tt"))
colnames(dt.plot)[2] = "celltype"
dt.plot = as.data.table(dt.plot)
dt.plot[, combination := paste0(tt,"_",celltype)]

tmp.res = res[grep('cov', name), ][,.(combination,tt,celltype,estimate,ci.lower,ci.upper,p,fdr,plot.label)]
colnames(tmp.res)[1] = "combination.formal"
tmp.res[, combination := paste0(tt,"_",celltype)]

tmp.res$tt = NULL
tmp.res$celltype = NULL

dt.plot = merge(dt.plot, tmp.res, by = "combination", all.x = TRUE)
unique(dt.plot$plot.label)

##dt.plot$celltype = factor(dt.plot$celltype, levels = c(lung.cluster, stromal.cluster))
dt.plot$tumor.type = "LUAD"
dt.plot$tumor.type = toupper(dt.plot$tumor.type)
dt.plot$combo.label = paste0(dt.plot$tumor.type,"_",dt.plot$celltype)

## center values at 1
dt.plot$estimate = dt.plot$estimate - 1
dt.plot$ci.lower = dt.plot$ci.lower - 1
dt.plot$ci.upper = dt.plot$ci.upper - 1


## add annotations for distal and proximal
dt.plot$Cell_Type = dt.plot$celltype
#dt.plot$Cell_Type = ifelse(dt.plot$celltype %in% c('Alveolar.Type.I.cells','Alveolar.Type.II.cells'),'Distal Lung Cell Types', 'Proximal Lung Cell Types')
#dt.plot$Cell_Type = factor(dt.plot$Cell_Type, levels = c('Distal Lung Cell Types', 'Proximal Lung Cell Types'))


## Relative risk barplots
in.pattern.tt = c('LUAD')
in.pattern.snv = "snv.count"
##in.pattern.snv = "snv.tob.count"

ppdf(
  for (ik in in.pattern.tt)
  {
    #ik = "HMF.SCLC"
    message('\nPlotting mutational density for: ', ik, ' samples')
    in.ct.sub = epicells
    in.pattern = in.pattern.snv
    in.count.sub = grep(paste0(in.pattern),dt.plot$tt, value = TRUE)
    dt.plot.sub = dt.plot[tumor.type %in% ik & celltype %in% in.ct.sub & tt %in% in.count.sub,]
    setorder(dt.plot.sub,celltype)
    print(dt.plot.sub)
    in.title = unique(gsub(" - .*","",dt.plot.sub$combination.formal))
    ## set upper and lower limits
    in.axis.breaks = 0.05
    in.axis.min = -0.25
    in.axis.max = 10.05
    
    p = ggplot(dt.plot.sub, aes(x = reorder(celltype,estimate), y = estimate)) +
      geom_bar(stat = "identity", alpha = 0.9) +
      scale_fill_manual(values = c("grey")) +
      geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.85, width = 0.25, color = 'black') +
      #scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(in.axis.min, in.axis.max)) +
      scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1) +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      xlab("") +
      ylab("Relative Risk") +
      ggtitle(paste0(in.title)) +
      theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) +
      coord_flip()
    print(p)
  }, file= "new_non_apobec_LUAD.pdf",cex = c(0.95,1.25))


## Fig 2B
juan.genes.gr = readRDS("~/data/scLung/Sikkema/infraobjs/DeconstructSigsV3_MutDensity_Patients_SikkCentroids_Vs_Genes_All_Celltypes_PreAvgExpr.rds")
count.int = grep("snv.non.apobec.count", colnames(mcols(juan.genes.gr)), value = TRUE)
count.int = as.data.frame(count.int)
count.int$pair = gsub("snv.non.apobec.count.","",count.int$count.int)
count.int$pair = toupper(count.int$pair)
#luad_pairs = readRDS("~/data/scLung/Sikkema/pairs/new.pairs.luad.rds")
lusc_pairs = readRDS("~/data/scLung/Sikkema/pairs/new.pairs.lusc.rds")
juan.genes.dt = gr2dt(juan.genes.gr)
juan.genes.dt = as.data.frame(juan.genes.dt)

#LUAD_patnt = pairs.in[which(unlist(pairs.in$tumor_type)=="LUSC"),]
LUAD_patnt_mat = intersect(lusc_pairs,count.int$pair)
LUAD_patnt_mat = as.matrix(LUAD_patnt_mat)
LUAD_patnt_mat = LUAD_patnt_mat[complete.cases(LUAD_patnt_mat),]
count.int.mat = count.int[match(LUAD_patnt_mat,count.int$pair),]
juan.genes.dt_nw = juan.genes.dt[,match(count.int.mat$count.int,colnames(juan.genes.dt))]
snv_count = as.matrix(rowSums(juan.genes.dt_nw))
juan_dt1 = juan.genes.dt[,1:6]
juan_dt1$snv.count = snv_count
sikk.cent = readRDS("~/Projects/COO/sikk_cent.rds")
sikk.cent_nw = sikk.cent[match(juan_dt1$gene_name,rownames(sikk.cent)),]
juan_dt2 = cbind(juan_dt1,sikk.cent_nw,juan.genes.dt[,6558:6562])
juan_dt2 = as.data.table(juan_dt2)
juan.genes.gr.nw = dt2gr(juan_dt2)

#kofi_centroids = sikk.cent
centroids= epicells
juan.genes.gr.nw$avgexpr = apply(mcols(juan.genes.gr.nw)[, colnames(mcols(juan.genes.gr.nw)) %in% centroids],1,mean)
juan.genes.dt1 = gr2dt(juan.genes.gr.nw)
juan.genes.dt1 = juan.genes.dt1[!is.na(gc),]
juan.genes.dt1 = juan.genes.dt1[!is.na(rept),]
#juan.genes.dt = dt2gr(juan.genes.dt)

## include only lung cell types for individual samples
cmut = grep('snv.count', names(mcols(juan.genes.gr.nw)), value = TRUE)

cts = epicells
combos = expand.grid(cmut, cts) %>% as.data.table


res = mclapply(1:nrow(combos), function(i)
{
  tt = combos[i,]$Var1 %>% as.character;
  celltype = combos[i,]$Var2 %>% as.character
  message(paste(tt, celltype));
  dat = data.table(snv = mcols(juan.genes.gr.nw)[[tt]],
                   gc = juan.genes.gr.nw$gc,
                   rept = juan.genes.gr.nw$rept,
                   avg = juan.genes.gr.nw$avgexpr,
                   cov =  mcols(juan.genes.gr.nw)[[celltype]],
                   exon.frac = log(juan.genes.gr.nw$exon.frac),
                   wid = log(width(juan.genes.gr.nw)))
  if(sum(dat$snv) < 1)
  {
    message(" -- No mutation data -- ")
    res = data.frame(name = c(NA,NA), method = c(NA,NA), p = c(1,1), estimate = c(1,1), ci.lower = c(1,1), ci.upper = c(1,1), effect = c(1, 1))
  } else {
    #res = glm(snv ~ cov + exon.frac + avg + rept + gc + offset(wid), family ="poisson",data = dat) %>% dflm
    res = glm.nb(snv ~ cov + avg + exon.frac + rept + gc + offset(wid), data = dat) %>% dflm
  }
  res$tt = tt
  res$celltype = celltype
  return(res)
}, mc.cores = 35, mc.preschedule = TRUE)

res = rbindlist(res)
res[, fdr := p.adjust(p), by = name] ## include BH correction
res = readRDS("~/data/scLung/Sikkema/resforfirsquartile.rds")
#res$tt = paste0(res$tt,".LUSC")


## create grouped naming annotation
res$plot.label = "LUSC"
res$plot.label = gsub("snv.count.","",res$plot.label)
#res$plot.label = paste0(toupper(gsub("_num_"," (",res$plot.label)),')')

res$combination = paste0(toupper(gsub("_num_"," (",res$tt)), ')')
res$combination = gsub(")","",res$combination)
res$combination = paste0(res$combination," - ", res$celltype)
#res$tt = gsub("_num.*","",res$tt)
unique(res$plot.label)


## set estimates value for each celltype below fdr threshold to null value of 1
dt = dcast.data.table(res[grep('cov', name), ][, estimate2 := ifelse(fdr<0.1, estimate, 1)], tt ~ celltype, value.var = 'estimate2')

dt.plot = melt(dt, id.vars = c("tt"))
colnames(dt.plot)[2] = "celltype"
dt.plot = as.data.table(dt.plot)
dt.plot[, combination := paste0(tt,"_",celltype)]

tmp.res = res[grep('cov', name), ][,.(combination,tt,celltype,estimate,ci.lower,ci.upper,p,fdr,plot.label)]
colnames(tmp.res)[1] = "combination.formal"
tmp.res[, combination := paste0(tt,"_",celltype)]

tmp.res$tt = NULL
tmp.res$celltype = NULL

dt.plot = merge(dt.plot, tmp.res, by = "combination", all.x = TRUE)
unique(dt.plot$plot.label)

##dt.plot$celltype = factor(dt.plot$celltype, levels = c(lung.cluster, stromal.cluster))
dt.plot$tumor.type = "LUSC"
dt.plot$tumor.type = toupper(dt.plot$tumor.type)
dt.plot$combo.label = paste0(dt.plot$tumor.type,"_",dt.plot$celltype)

## center values at 1
dt.plot$estimate = dt.plot$estimate - 1
dt.plot$ci.lower = dt.plot$ci.lower - 1
dt.plot$ci.upper = dt.plot$ci.upper - 1


## add annotations for distal and proximal
dt.plot$Cell_Type = dt.plot$celltype
#dt.plot$Cell_Type = ifelse(dt.plot$celltype %in% c('Alveolar.Type.I.cells','Alveolar.Type.II.cells'),'Distal Lung Cell Types', 'Proximal Lung Cell Types')
#dt.plot$Cell_Type = factor(dt.plot$Cell_Type, levels = c('Distal Lung Cell Types', 'Proximal Lung Cell Types'))


## Relative risk barplots
in.pattern.tt = c('LUSC')
in.pattern.snv = "snv.count"
##in.pattern.snv = "snv.tob.count"

ppdf(
  for (ik in in.pattern.tt)
  {
    #ik = "HMF.SCLC"
    message('\nPlotting mutational density for: ', ik, ' samples')
    in.ct.sub = epicells
    in.pattern = in.pattern.snv
    in.count.sub = grep(paste0(in.pattern),dt.plot$tt, value = TRUE)
    dt.plot.sub = dt.plot[tumor.type %in% ik & celltype %in% in.ct.sub & tt %in% in.count.sub,]
    setorder(dt.plot.sub,celltype)
    print(dt.plot.sub)
    in.title = unique(gsub(" - .*","",dt.plot.sub$combination.formal))
    ## set upper and lower limits
    in.axis.breaks = 0.05
    in.axis.min = -0.25
    in.axis.max = 10.05
    
    p = ggplot(dt.plot.sub, aes(x = reorder(celltype,estimate), y = estimate)) +
      geom_bar(stat = "identity", alpha = 0.9) +
      scale_fill_manual(values = c("grey")) +
      geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.85, width = 0.25, color = 'black') +
      #scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1, limits = c(in.axis.min, in.axis.max)) +
      scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1) +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      xlab("") +
      ylab("Relative Risk") +
      ggtitle(paste0(in.title)) +
      theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) +
      coord_flip()
    print(p)
  }, file= "new_non_apobec_LUSC.pdf",cex = c(0.95,1.25))

## Supplementary Fig 10 A

tmp.plot = readRDS("~/Projects/COO/tmp_plot_for_umap_fig5_knn6.rds")

#Adjust carcinoma annotations. To be 100% sure we're doing it right, I'm using a cycle.
for(i in 1:nrow(suk.annot)){
  tmp.plot$celltype.ident[tmp.plot$cell==rownames(suk.annot)[i]] = suk.annot$ann_finest_lev_transferred_label[i]
}

highlight_color <- "darkslategray3"  # Color for WCM-1
grey_color <- "grey"        # Color for other patients

tmp.plot$color <- ifelse(tmp.plot$patient == "WCM-1", highlight_color, grey_color)

p <- ggplot(tmp.plot, aes(x, y, colour = color), size = 3) + 
  geom_point() + 
  theme_bw() + 
  theme(legend.position = "top") +
  scale_colour_identity()

ppdf(print(p),filename="UMAP with WCM-1 patient knn6_nw1.pdf")

## Supplementary Fig 10 B
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(tmp.plot$cluster)))

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

## Fig 5B Left Panel ##
#Load libraries
library(skitools)
library(parallel)
library(MASS)

#Load csv with relative risk information per centroid for every sample.
rel_risk = read.table("/gpfs/commons/home/jandrademartínez/luad_rel_risk_SP.csv",sep=",",header=T,row.names=1)

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
wcm = readRDS("/gpfs/commons/home/jandrademartínez/cov_res.rds")
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
#Results are saved in the distances data.frame object (figure was manually maded based on this object).
distances = as.data.frame(distances)
distances$Vector = rownames(distances)

## Fig 5B Right Panel ##
#Load libraries
library(skitools)

#Load GLM output for WCM-1.
dt.plot = readRDS("/gpfs/commons/home/jandrademartínez/cov_res.rds")
dt.plot$Cell_Type = dt.plot$celltype
dt.plot$tumor.type = "LUAD"

#Center relative risk and confidence interval estimate values at 1 for plot.
dt.plot$estimate = dt.plot$estimate - 1
dt.plot$ci.lower = dt.plot$ci.lower - 1
dt.plot$ci.upper = dt.plot$ci.upper - 1

# Relative risk barplots. The code allows for filtering for specific cancer type, snv columns, and celltypes via in.pattern.tt, in.pattern.snv, and in.ct.sub respectively. Here we add the values needed for plotting all GLM outputs defined above with no restriction. Note that multiple cancer subtype grpahs can be done by adding more values to in.pattern.tt. 
in.pattern.tt = c('LUAD')
in.pattern.snv = "snv.count"

ppdf(
  for (ik in in.pattern.tt)
  {
    message('\nPlotting mutational density for: ', ik, ' samples')
#Filter for any cancer type, celltype, or snv column if needed.     
    in.ct.sub = dt.plot$celltype
    in.pattern = in.pattern.snv
    in.count.sub = grep(paste0(in.pattern),dt.plot$tt, value = TRUE)
    dt.plot.sub = dt.plot[tumor.type %in% ik & celltype %in% in.ct.sub & tt %in% in.count.sub,]
#Order results by relative risk (lower to higher).         
    setorder(dt.plot.sub,celltype)
    print(dt.plot.sub)
    in.title = unique(gsub(" - .*","",dt.plot.sub$combination.formal))
# Set upper and lower limits for plot (x axis).
    in.axis.breaks = 0.05
    in.axis.min = -0.8
    in.axis.max = 10.05
#Construct ggplot2 object and plot.        
    p = ggplot(dt.plot.sub, aes(x = reorder(celltype,estimate), y = estimate)) +
      geom_bar(stat = "identity", alpha = 0.9) +
      scale_fill_manual(values = c("grey")) +
      geom_errorbar(aes(ymax = ci.upper, ymin = ci.lower), size = 0.85, width = 0.25, color = 'black') +
      scale_y_continuous(breaks = seq(in.axis.min-0.05, in.axis.max+0.05, in.axis.breaks), labels = function(y) y + 1) +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      xlab("") +
      ylab("Relative Risk") +
      ggtitle(paste0(in.title)) +
      theme(axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1)) +
      coord_flip()
    print(p)
  }, file= "LUAD_COO_NB_FDR_0_1_all_celltype_WCM1_nw.pdf",cex = c(0.95,1.25))

## Fig 5C ##
#Load libraries. 
library(skitools)
#Load file with umap of cancer cells output. This allows for gathering WCM-1A and B cells. 
tmp.plot = readRDS("~/../spanja/Projects/COO/tmp_plot_for_umap_fig5_knn6.rds")

#Filter for WCM-1 cells, split in WCM-1 A and B based on their cluster numbers.
tmp_plot_wcm = tmp.plot[which(tmp.plot$patient=="WCM-1"),]
wcm_clusters = tmp_plot_wcm[which(tmp_plot_wcm$cluster!=1),]
wcm_clusters[which(wcm_clusters$cluster==2),"type"]="WCM-1-A"
wcm_clusters[which(wcm_clusters$cluster==3),"type"]="WCM-1-B"
wcm_clusters = as.data.frame(wcm_clusters)
wcm1a = wcm_clusters[which(wcm_clusters$type=="WCM-1-A"),]
wcm1b = wcm_clusters[which(wcm_clusters$type=="WCM-1-B"),]

#Load aneuploidy score output from CONICSmat.
tmp_norm = readRDS("~/../spanja/Projects/COO/anneuploidy zscore.rds")
tmp_norm = as.data.frame(tmp_norm)

#Get aneuploidy score for WCM-1 cells.
tmp_norm_wcm = tmp_norm[match(wcm_clusters$cell,rownames(tmp_norm)),]
tmp_norm_wcm = as.matrix(tmp_norm_wcm)
tmp_norm_wcm1a = tmp_norm_wcm[match(wcm1a$cell,rownames(tmp_norm_wcm)),]
tmp_norm_wcm1b = tmp_norm_wcm[match(wcm1b$cell,rownames(tmp_norm_wcm)),]
tmp_norm_wcm = rbind(tmp_norm_wcm1a,tmp_norm_wcm1b)
tmp_norm_wcm = as.matrix(tmp_norm_wcm)

#Generate heatmap annotation and then plot heatmap.
row.ha = HeatmapAnnotation(Patient = wcm_clusters$type, name = 'Tissue', width = unit(2, "cm"), show_legend = TRUE, show_annotation_name = FALSE, col = list(Patient = c("WCM-1-A"= "#00BA38", "WCM-1-B" = "#619CFF")), annotation_legend_param = list("WCM-1-A", "WCM-1-B"), which = "row")
p = Heatmap(tmp_norm_wcm, name = 'z-score', cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 8), column_names_side = c("bottom"), show_column_names = TRUE, show_row_names = FALSE, column_title_gp = gpar(fontsize = 20, fontface = "bold"))
ppdf(print(p),filename="Fig5C.pdf")

## Fig 5D ##
#Load libraries 
library(skitools)

#Load per cell identity labels. Remove cells with no identity.
labels = read.table("~/../spanja/Projects/COO/combined_emb_for_carc_for_homogenous_sikk_ep_500_nw.csv",sep=",",header=TRUE)
labels$ann_finest_lev_transferred_label_filtered = ifelse(labels$ann_finest_lev_transfer_uncert<0.5,labels$ann_finest_lev_transferred_label_unfiltered,"Unknown")

#Define proximal and distal cell groupings.
distal_sikk = c("AT1","AT2","AT2 proliferating","AT0")
rare_sikk = c("Neuroendocrine","Tuft","Ionocyte")
cells = c(distal_sikk,rare_sikk)
proximal_sikk = setdiff(labels$ann_finest_lev_transferred_label_filtered,cells)

#Get WCM-1 A and B cells from the carcinoma cell UMAP  (same as for Fig 5C).
tmp.plot = readRDS("~/Projects/COO/tmp_plot_for_umap_fig5_knn6.rds")
tmp_plot_wcm = tmp.plot[which(tmp.plot$patient=="WCM-1"),]
wcm_clusters = tmp_plot_wcm[which(tmp_plot_wcm$cluster!=1),]
wcm_clusters[which(wcm_clusters$cluster==2),"type"]="WCM-1-A"
wcm_clusters[which(wcm_clusters$cluster==3),"type"]="WCM-1-B"
wcm_clusters = as.data.frame(wcm_clusters)
wcm1a = wcm_clusters[which(wcm_clusters$type=="WCM-1-A"),]
wcm1b = wcm_clusters[which(wcm_clusters$type=="WCM-1-B"),]

#For WCM-1 A and B, get the cell labels and assign Distal/Non-Distal to the corresponding celltypes.
wcm1a_labels = labels[match(wcm1a$cell,labels$CellBarcode),]
wcm1a_labels = wcm1a_labels[which(wcm1a_labels$ann_finest_lev_transferred_label_filtered!="Unknown"),]
wcm1a_DND = as.data.frame(table(wcm1a_labels$ann_finest_lev_transferred_label_filtered))
wcm1a_DND$pat = "WCM-1-A"
wcm1a_DND$cat = "Non-Distal"
wcm1a_DND[which(wcm1a_DND$Var1 %in% distal_sikk),"cat"] = "Distal"

wcm1b_labels = labels[match(wcm1b$cell,labels$CellBarcode),]
wcm1b_labels = wcm1b_labels[which(wcm1b_labels$ann_finest_lev_transferred_label_filtered!="Unknown"),]
wcm1b_DND = as.data.frame(table(wcm1b_labels$ann_finest_lev_transferred_label_filtered))
wcm1b_DND$pat = "WCM-1-B"
wcm1b_DND$cat = "Non-Distal"
wcm1b_DND[which(wcm1b_DND$Var1 %in% distal_sikk),"cat"] = "Distal"

#Combine and compute percentage per celltype per subpopulation..
wcm1_DND = rbind(wcm1a_DND,wcm1b_DND)
result <- wcm1_DND %>%
  group_by(pat) %>%
  mutate(percentage = Freq / sum(Freq) * 100)

#Define color palette and plot.
library(RColorBrewer)
palette1 = brewer.pal(12, "Set3")
palette2 = brewer.pal(5, "Set2")
combined_palette = c(palette1, palette2)

p = ggplot(result, aes(x = factor(cat), y = percentage, fill = Var1)) +
  geom_bar(stat = "identity") +                # Create barplot
  theme_bw() +                                 # Use a clean theme
#  theme(legend.position = "none") +            # Remove legend
  facet_wrap(~ pat) +         # Create separate barplots for each 'pat'
  labs(x = "Frequency", y = "Percentage",      # Label axes
       title = "Barplot of Percentage by Frequency and Pat") + 
  scale_fill_manual(values = combined_palette)+ylim(0,100)

ppdf(print(p),filename="barplot for WCM1 A and B distal non-distal_v1.pdf")

##Fig 5E ##
#load libraries
library(skitools)
library(EnhancedVolcano)

#Load file with marker genes.
seu.markers = readRDS("~/../spanja/Projects/COO/seu_markers.rds")

#Filter for markers from WCM-1 A and B (based on their cluster number in the carcinoma cell umap).
markers.1 = seu.markers[seu.markers$cluster == 2,]
markers.1$wcm = 'WCM-1-A'
markers.1$avg_log2FC = markers.1$avg_log2FC * -1
markers.2 = seu.markers[seu.markers$cluster == 3,]
markers.2$wcm = 'WCM-1-B'

#Combine the markers and plot.
seu.markers.mast.plot = rbind(markers.1, markers.2)
p = EnhancedVolcano(seu.markers.mast.plot,
                    lab = seu.markers.mast.plot$gene,
                    x = 'avg_log2FC', y = 'p_val',
                    xlim = c(-3, 3),
                    title = 'WCM-1-A (Alveolar-like) v. WCM-1-B (Basal-like)',
                    subtitle = "",
                    pCutoff = 10e-10,
                    FCcutoff = 0.5,
                    labSize = 3, colAlpha = 1,
                    pointSize = 2)
ppdf(print(p),filename="Fig5E.pdf" ,cex = 0.8)


## EDF 7 -- Univariate ##
#Load libraries.
library(skitools)
library(MASS)
library(parallel)

#Load GRanges file with genomic covariates and centroids.
juan.genes.gr.nw = readRDS("~/../spanja/Projects/COO/genes_gr_LUAD.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)
#Define cells for simulation (epithelial cells).
range_Cells = c("AT1","AT2","AT2 proliferating","AT0","Suprabasal","Basal resting","Hillock-like","Multiciliated (non-nasal)","Multiciliated (nasal)","Deuterosomal","Neuroendocrine","Ionocyte","Tuft","Goblet (nasal)","Club (nasal)","Club (non-nasal)","pre-TB secretory","Goblet (bronchial)","Goblet (subsegmental)","SMG serous (bronchial)","SMG mucous","SMG duct","SMG serous (nasal)")


#Simulate COO Centroids.
simCent.dt = data.table()
#Define number of transcripts per cell.
txpercell = 4500
#For each centroid to be used (epithelial cell centroids), compute the fraction of total counts associated to each gene. The product of this quantity times txpercell is the expected number of counts for said gene. 
for(i in 1:length(range_Cells)){
  thisCell = range_Cells[i]
  centroid.col = gr2dt(juan.genes.gr.nw[,thisCell])[,6]
  centroid.col = exp(centroid.col)
  rel_exp_per_gene = centroid.col/sum(centroid.col)
  lambda_per_gene = as.matrix(rel_exp_per_gene*txpercell)
#To get the counts per gene for the simulated centroid, sample from a poisson distribution per gene, where the lambda is the expected number of counts for each gene. Log this quantity for the GLM and save.  
  simCent.dt[, newcol := log(rpois(nrow(centroid.col), lambda_per_gene) + 0.001)]
  colnames(simCent.dt)[i] = thisCell
}

#Join simulated centroid values per gene with the GRanges with covariate information per gene.
juan.dt =  gr2dt(juan.genes.gr.nw)
juan.dt = juan.dt[,1:6]
simCent.dt = as.data.table(simCent.dt)
simCent.dt_nw = cbind(juan.dt,simCent.dt)
simCent.dt.gr = dt2gr(simCent.dt_nw)
#Compute average simulated centroid expression (needed as covariate for univariable GLM) .
epicells = range_Cells
juan.genes.gr.nw$simavgexpr = apply(mcols(simCent.dt.gr)[, colnames(mcols(simCent.dt.gr)) %in% epicells],1,mean)

#Set range of TMB values for simulation.
rval = c(seq(1000, 3000, by = 100),seq(10, 90, by = 10),seq(1,9,1))

#For each value in rval (i.e.: for each TMB value), we run 100 iterations of the simulation.
for(i in 1:length(rval))
{
 rval1 = rval[i]
 print(rval1)
 sim_COO_mat = mclapply(1:100, function(x)
  {
#Per simulation, set the corresponding TMB.      
    print(x)
    randTMB = rval1
#Randomly sample one of the epithelial cells to act as the true celltype.      
  thisCell = sample(range_Cells, 1)
  expVec = gr2dt(juan.genes.gr.nw[,thisCell])[,6]
  juan.genes.gr.nw.dt$expr1 = exp(expVec)
#Set values for fraction of lineage specific mutations (thislfrac) and correlation strength between expression and TMB (thisalpha).      
  thislfrac = 0.5
  thisalpha = -0.65
  juan.genes.gr.nw.dt[, paste0("exprSim") := expr1]

#Simulate TMB per gene for this simulation iteration. First, we shift the replication timing vector values so that the minimum value is 0 (negative or 0 values lead to errors in TMB simulation).        
  juan.genes.gr.nw.dt[["reptAdj"]] = juan.genes.gr.nw.dt[["rept"]] + abs(min(juan.genes.gr.nw.dt[["rept"]]))+1
#Generate a base vector representing the combined effect of genomic covariates (gc content, replication timing, and exon fraction) on the TMB of each gene).
  base_lambda = (nrow(juan.genes.gr.nw.dt)*juan.genes.gr.nw.dt$gc^(-1.0)/sum((juan.genes.gr.nw.dt$gc)^(-1.0)))*(nrow(juan.genes.gr.nw.dt)*juan.genes.gr.nw.dt$reptAdj^(-0.06)/sum((juan.genes.gr.nw.dt$reptAdj)^(-0.06)))*(nrow(juan.genes.gr.nw.dt)*juan.genes.gr.nw.dt$exon.frac^(-0.16)/sum((juan.genes.gr.nw.dt$exon.frac)^(-0.16)))
#Generate a vector representing the effect of gene expression on the TMB of each gene.     
  exp_sim = (nrow(juan.genes.gr.nw.dt)*((juan.genes.gr.nw.dt$exprSim)^(thisalpha))/sum((juan.genes.gr.nw.dt$exprSim)^(thisalpha)))
#Get the expected TMB value of each gene based on target TMB and expression/covariate effects.     
  tmb_nw = randTMB*juan.genes.gr.nw.dt$width*10^-6
#Sample TMB per gene using a poisson distribution where lambda is the expected value of TMB per gene.     
  juan.genes.gr.nw.dt$sim.muts = rpois(nrow(juan.genes.gr.nw.dt), (tmb_nw*((1-thislfrac)+thislfrac*exp_sim))*base_lambda)

#Run GLM. Set the combinations of centroid expression and TMB.      
  cmut = grep('sim.muts', names(juan.genes.gr.nw.dt), value = TRUE)
  cts = epicells
  combos = expand.grid(cmut, cts) %>% as.data.table
  res = mclapply(1:nrow(combos), function(i)
  {
    celltype = combos[i,]$Var2 %>% as.character
#Generate a data.frame with all values needed for the model (simulated TMB, centroid expression, and covariates).
    sim.dat = data.table(snv = juan.genes.gr.nw.dt$sim.muts,
                         gc = juan.genes.gr.nw.dt$gc,
                         rept = juan.genes.gr.nw.dt$rept,
                         avg = mcols(juan.genes.gr.nw)[["simavgexpr"]],
                         cov =  mcols(simCent.dt.gr)[[celltype]],
                         exon.frac = log(juan.genes.gr.nw.dt$exon.frac),
                         wid = log(width(juan.genes.gr.nw)))
#Run model.    
    if(sum(sim.dat$snv) < 1)
    {
      message(" -- No mutation data -- ")
      res = data.frame(name = c(NA,NA), method = c(NA,NA), p = c(1,1), estimate = c(1,1), ci.lower = c(1,1), ci.upper = c(1,1), effect = c(1, 1))
    } else {
      res = glm.nb(snv ~ cov + exon.frac + avg + rept + gc + offset(wid), data = sim.dat) %>% dflm
    }
    res$tt = "snv.sim.count"
    res$celltype = celltype
    return(res)
  }, mc.cores = 2, mc.preschedule = TRUE)

#Combine results of all GLMs for all centroids.
  res = rbindlist(res)
#Adjust p-value and filter covariates with fdr higher than 0.1.
  res[, fdr := p.adjust(p,method="BH"), by = name]
  res = res[which(res$fdr<0.1),]
  res$plot.label = "LUAD"

#Gather top celltype and those overlapping it's relative risk confidence interval. 
  allOverlap = "No_Overlap"
#Filter results for only expression rows and sort by relative risk estimate.      
  cov_res = res[grep("cov",name),]
  sorted_cov_res = cov_res[order(cov_res$estimate),]
#Compute confidence interval length for top cell and check which other confidence intervals for other cells, if any, overlap it by more than 60% of thee latter's length.
    calc.num = sorted_cov_res$ci.upper[1] - sorted_cov_res$ci.lower
    calc.denom = sorted_cov_res$ci.upper[1] - sorted_cov_res$ci.lower[1]
    overlap.calc =  calc.num/ calc.denom
    overlap.calc = as.matrix(overlap.calc)
    rownames(overlap.calc)= sorted_cov_res$celltype
    overlap.calc.cells = overlap.calc[which(overlap.calc[,1]>0.6),]
#Save top and overlapping cells separated by commas. Join with true celltype in data.frame and return result for this iteration.     
    allOverlap = paste0(names(overlap.calc.cells),collapse=",")
    overlap.calc.cells = as.matrix(overlap.calc.cells)
    sim_COO_mat = as.data.frame(cbind(thisCell,allOverlap))
    return(sim_COO_mat)
   }, mc.cores = 8, mc.preschedule = TRUE)

 
#Define function to check for errors in GLM output. This is needed for low TMB sims in case the model doesn't converge.
replace_errors <- function(x, replacement) {
  if (inherits(x, "try-error")) {
    return(replacement)
  } else {
    return(x)
  }
 }

#Set a replacement row for non-converging results. Since no convergence is reached, there's no overlap. Arbitrarily assign a celltype for this cases (since the allOverlap column is only equal to No_Overlap in this case, we can account this non-convergent results when computing accuracy).
replacement_value = data.table(thisCell="Deuterosomal",allOverlap="No_Overlap")

#Apply the function to each element in the list of results.
  fixed_list <- lapply(sim_COO_mat, replace_errors, replacement = replacement_value)

#Combine the results of the runs for this TMB and save, add TMB value to rows, and save.  
  sim_COO_mat_v1 = rbindlist(fixed_list)
  toSave = sim_COO_mat_v1 
  toSave$TMB = rval1
  thisRDS = readRDS("~/projects/scLung/db/testRunEDF6_withOverlap_Cov.rds")
  thisRDS = rbind(thisRDS,toSave)
  saveRDS(thisRDS,"~/projects/scLung/db/testRunEDF6_withOverlap_Cov.rds")

}

## EDF 7 -- Multivariate ##
#Load libraries.
library(skitools)
library(MASS)
library(parallel)

#Most of the code is the same as for univariate. Comments are left as is for those parts, and changed when neccesary.
#Load GRanges file with genomic covariates and centroids.
juan.genes.gr.nw = readRDS("~/../spanja/Projects/COO/genes_gr_LUAD.rds")
juan.genes.gr.nw.dt = gr2dt(juan.genes.gr.nw)
#Define cells for simulation (epithelial cells).
range_Cells = c("AT1","AT2","AT2 proliferating","AT0","Suprabasal","Basal resting","Hillock-like","Multiciliated (non-nasal)","Multiciliated (nasal)","Deuterosomal","Neuroendocrine","Ionocyte","Tuft","Goblet (nasal)","Club (nasal)","Club (non-nasal)","pre-TB secretory","Goblet (bronchial)","Goblet (subsegmental)","SMG serous (bronchial)","SMG mucous","SMG duct","SMG serous (nasal)")


#Simulate COO Centroids.
simCent.dt = data.table()
#Define number of transcripts per cell.
txpercell = 4500
#For each centroid to be used (epithelial cell centroids), compute the fraction of total counts associated to each gene. The product of this quantity times txpercell is the expected number of counts for said gene. 
for(i in 1:length(range_Cells)){
  thisCell = range_Cells[i]
  centroid.col = gr2dt(juan.genes.gr.nw[,thisCell])[,6]
  centroid.col = exp(centroid.col)
  rel_exp_per_gene = centroid.col/sum(centroid.col)
  lambda_per_gene = as.matrix(rel_exp_per_gene*txpercell)
#To get the counts per gene for the simulated centroid, sample from a poisson distribution per gene, where the lambda is the expected number of counts for each gene. Log this quantity for the GLM and save.  
  simCent.dt[, newcol := log(rpois(nrow(centroid.col), lambda_per_gene) + 0.001)]
  colnames(simCent.dt)[i] = thisCell
}

#Join simulated centroid values per gene with the GRanges with covariate information per gene.
juan.dt =  gr2dt(juan.genes.gr.nw)
juan.dt = juan.dt[,1:6]
simCent.dt = as.data.table(simCent.dt)
simCent.dt_nw = cbind(juan.dt,simCent.dt)
simCent.dt.gr = dt2gr(simCent.dt_nw)
#Compute average simulated centroid expression (needed as covariate for univariable GLM) .
epicells = range_Cells
juan.genes.gr.nw$simavgexpr = apply(mcols(simCent.dt.gr)[, colnames(mcols(simCent.dt.gr)) %in% epicells],1,mean)

#Set range of TMB values for simulation.
rval = c(seq(1000, 3000, by = 100),seq(10, 90, by = 10),seq(1,9,1))

#For each value in rval (i.e.: for each TMB value), we run 100 iterations of the simulation.
for(i in 1:length(rval))
{
 rval1 = rval[i]
 print(rval1)
 sim_COO_mat = mclapply(1:100, function(x)
  {
#Per simulation, set the corresponding TMB.      
    print(x)
    randTMB = rval1
#Randomly sample one of the epithelial cells to act as the true celltype.      
  thisCell = sample(range_Cells, 1)
  expVec = gr2dt(juan.genes.gr.nw[,thisCell])[,6]
  juan.genes.gr.nw.dt$expr1 = exp(expVec)
#Set values for fraction of lineage specific mutations (thislfrac) and correlation strength between expression and TMB (thisalpha).      
  thislfrac = 0.5
  thisalpha = -0.65
  juan.genes.gr.nw.dt[, paste0("exprSim") := expr1]

#Simulate TMB per gene for this simulation iteration. First, we shift the replication timing vector values so that the minimum value is 0 (negative or 0 values lead to errors in TMB simulation).        
  juan.genes.gr.nw.dt[["reptAdj"]] = juan.genes.gr.nw.dt[["rept"]] + abs(min(juan.genes.gr.nw.dt[["rept"]]))+1
#Generate a base vector representing the combined effect of genomic covariates (gc content, replication timing, and exon fraction) on the TMB of each gene).
  base_lambda = (nrow(juan.genes.gr.nw.dt)*juan.genes.gr.nw.dt$gc^(-1.0)/sum((juan.genes.gr.nw.dt$gc)^(-1.0)))*(nrow(juan.genes.gr.nw.dt)*juan.genes.gr.nw.dt$reptAdj^(-0.06)/sum((juan.genes.gr.nw.dt$reptAdj)^(-0.06)))*(nrow(juan.genes.gr.nw.dt)*juan.genes.gr.nw.dt$exon.frac^(-0.16)/sum((juan.genes.gr.nw.dt$exon.frac)^(-0.16)))
#Generate a vector representing the effect of gene expression on the TMB of each gene.     
  exp_sim = (nrow(juan.genes.gr.nw.dt)*((juan.genes.gr.nw.dt$exprSim)^(thisalpha))/sum((juan.genes.gr.nw.dt$exprSim)^(thisalpha)))
#Get the expected TMB value of each gene based on target TMB and expression/covariate effects.     
  tmb_nw = randTMB*juan.genes.gr.nw.dt$width*10^-6
#Sample TMB per gene using a poisson distribution where lambda is the expected value of TMB per gene.     
  juan.genes.gr.nw.dt$sim.muts = rpois(nrow(juan.genes.gr.nw.dt), (tmb_nw*((1-thislfrac)+thislfrac*exp_sim))*base_lambda)

#Run the GLM. Since only one model is needed no cycle is done here. Generate a data.frame with simulated TMB and covariates.
         dat = data.table(snv = juan.genes.gr.nw.dt$sim.muts,
                     gc = juan.genes.gr.nw$gc,
                     rept = juan.genes.gr.nw$rept,
                     avg = mcols(juan.genes.gr.nw)[["simavgexpr"]],
                     exon.frac = log(juan.genes.gr.nw$exon.frac),
                     wid = log(width(juan.genes.gr.nw)))
#Add simulated centroid expression for epithelial cells.
    dat = cbind(dat,gr2dt(simCent.dt.gr)[,7:29])
#Define the GLM formula. Note that average expression is not included in this formula since all centroids are already included (if it were to be added, then one of the centroids would be excluded due to being a linear combination of the average and the other centroids).
    multi_formula <- as.formula(paste0("snv ~ ", paste0(colnames(dat[, -c(1,4,6)]), collapse = " + "), " + offset(wid)") )
    if(sum(dat$snv) < 1)
    {
      message(" -- No mutation data -- ")
      res = data.frame(name = c(NA,NA), method = c(NA,NA), p = c(1,1), estimate = c(1,1), ci.lower = c(1,1), ci.upper = c(1,1), effect = c(1, 1))
    } else {
     res = glm.nb(multi_formula,  data = dat) %>% dflm
    }
#Format output columns to be able to identify expression rows below.
    res$tt = "snv.sim.counts"
    res$name[5:27] = epicells
    res$celltype = res$name
#Adjust p-value and filter covariates with fdr higher than 0.1.     
  res[, fdr := p.adjust(p,method="BH"), by = name]
  res = res[which(res$fdr<0.1),]
  res$plot.label = "LUAD"
      
#Gather top celltype and those overlapping it's relative risk confidence interval. 
  allOverlap = "No_Overlap"
#Filter results for only expression rows and sort by relative risk estimate.      
  cov_res = res[grep("cov",name),]
  sorted_cov_res = cov_res[order(cov_res$estimate),]
#Compute confidence interval length for top cell and check which other confidence intervals for other cells, if any, overlap it by more than 60% of thee latter's length.
    calc.num = sorted_cov_res$ci.upper[1] - sorted_cov_res$ci.lower
    calc.denom = sorted_cov_res$ci.upper[1] - sorted_cov_res$ci.lower[1]
    overlap.calc =  calc.num/ calc.denom
    overlap.calc = as.matrix(overlap.calc)
    rownames(overlap.calc)= sorted_cov_res$celltype
    overlap.calc.cells = overlap.calc[which(overlap.calc[,1]>0.6),]
#Save top and overlapping cells separated by commas. Join with true celltype in data.frame and return result for this iteration.     
    allOverlap = paste0(names(overlap.calc.cells),collapse=",")
    overlap.calc.cells = as.matrix(overlap.calc.cells)
    sim_COO_mat = as.data.frame(cbind(thisCell,allOverlap))
    return(sim_COO_mat)
   }, mc.cores = 8, mc.preschedule = TRUE)

 
#Define function to check for errors in GLM output. This is needed for low TMB sims in case the model doesn't converge.
replace_errors <- function(x, replacement) {
  if (inherits(x, "try-error")) {
    return(replacement)
  } else {
    return(x)
  }
 }

#Set a replacement row for non-converging results. Since no convergence is reached, there's no overlap. Arbitrarily assign a celltype for this cases (since the allOverlap column is only equal to No_Overlap in this case, we can account this non-convergent results when computing accuracy).
replacement_value = data.table(thisCell="Deuterosomal",allOverlap="No_Overlap")

#Apply the function to each element in the list of results.
  fixed_list <- lapply(sim_COO_mat, replace_errors, replacement = replacement_value)

#Combine the results of the runs for this TMB and save, add TMB value to rows, and save.  
  sim_COO_mat_v1 = rbindlist(fixed_list)
  toSave = sim_COO_mat_v1 
  toSave$TMB = rval1
  thisRDS = readRDS("~/projects/scLung/db/testRunEDF6_withOverlap_Cov.rds")
  thisRDS = rbind(thisRDS,toSave)
  saveRDS(thisRDS,"~/projects/scLung/db/testRunEDF6_withOverlap_Cov.rds")

}


## EDF 7 --Accuracy Plots ##
#Load libraries.
library(skitools)

#Plots can be made for univariate and multivariate simultaneously or individually. Comments are made when relevant to do one or the other.

#Load simulation result files. Add Method column to identify if results are uni or multivariate.
resDT = readRDS("~/projects/scLung/db/testRunEDF6_withOverlapMultivariate.rds")
resDT$Method = "Univariate"

#If results are plotted in the same graph, uncomment the following lines, adjusting the Method column value accordingly. 
#MulDT = readRDS("~/projects/scLung/db/testRunEDF6_withOverlap.rds")
# MulDT$Method = "Univariate"
#resDT = rbind(resDT,MulDT,fill=TRUE)

#Load csv with Sikkema annotations.
annotSikk = read.csv("/gpfs/commons/home/pmantri/projects/scLung/db/sikk_meta_epi3.csv")
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



rmarkdown::render(
        input = normalizePath("/gpfs/commons/groups/imielinski_lab/projects/scLung/draft_scLung.rmd"),
        output_format = "html_document",
        output_file = normalizePath("/gpfs/commons/groups/imielinski_lab/projects/scLung/draft_scLung.html"),
        knit_root_dir = normalizePath("/gpfs/commons/groups/imielinski_lab/projects/scLung/"),
        params = list(set_title = "Test",fileLoc = normalizePath("/gpfs/commons/groups/imielinski_lab/projects/scLung/db/scLung_RMDFiles/")),quiet = FALSE)

