## Figure 6B ##
```{r, fig6B, results = "asis", echo = TRUE, out.width= "100%", fig.align = "center"}
#Load libraries
library(skitools)

#Load GLM output for WCM-1.
dt.plot = readRDS("../data/Fig6/6B_cov_res.rds")
dt.plot$Cell_Type = dt.plot$celltype
dt.plot$tumor.type = "LUAD"

#Center relative risk and confidence interval estimate values at 1 for plot.
dt.plot$estimate = dt.plot$estimate - 1
dt.plot$ci.lower = dt.plot$ci.lower - 1
dt.plot$ci.upper = dt.plot$ci.upper - 1

# Relative risk barplots. The code allows for filtering for specific cancer type, snv columns, and celltypes via in.pattern.tt, in.pattern.snv, and in.ct.sub respectively. Here we add the values needed for plotting all GLM outputs defined above with no restriction. Note that multiple cancer subtype grpahs can be done by adding more values to in.pattern.tt. 
in.pattern.tt = c('LUAD')
in.pattern.snv = "snv.count"

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
}

## Figure 6C ##
#Load libraries. 
library(skitools)
#Load file with umap of cancer cells output. This allows for gathering WCM-1A and B cells. 
tmp.plot = readRDS("../data/Fig6/6CD_tmp_plot_for_umap_fig6_knn6.rds")

#Filter for WCM-1 cells, split in WCM-1 A and B based on their cluster numbers.
tmp_plot_wcm = tmp.plot[which(tmp.plot$patient=="WCM-1"),]
wcm_clusters = tmp_plot_wcm[which(tmp_plot_wcm$cluster!=1),]
wcm_clusters[which(wcm_clusters$cluster==2),"type"]="WCM-1-A"
wcm_clusters[which(wcm_clusters$cluster==3),"type"]="WCM-1-B"
wcm_clusters = as.data.frame(wcm_clusters)
wcm1a = wcm_clusters[which(wcm_clusters$type=="WCM-1-A"),]
wcm1b = wcm_clusters[which(wcm_clusters$type=="WCM-1-B"),]

#Load aneuploidy score output from CONICSmat.
tmp_norm = readRDS("../data/Fig6/Fig6C_anneuploidy_zscore.rds")
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
print(p)

## Figure 6D ##
#Load libraries 
library(skitools)

#Load per cell identity labels. Remove cells with no identity.
labels = read.table("../data/Fig6/6D_combined_emb_for_carc_for_homogenous_sikk_ep_500_nw.csv",sep=",",header=TRUE)
labels$ann_finest_lev_transferred_label_filtered = ifelse(labels$ann_finest_lev_transfer_uncert<0.5,labels$ann_finest_lev_transferred_label_unfiltered,"Unknown")

#Define proximal and distal cell groupings.
distal_sikk = c("AT1","AT2","AT2 proliferating","AT0")
rare_sikk = c("Neuroendocrine","Tuft","Ionocyte")
cells = c(distal_sikk,rare_sikk)
proximal_sikk = setdiff(labels$ann_finest_lev_transferred_label_filtered,cells)

#Get WCM-1 A and B cells from the carcinoma cell UMAP  (same as for Fig 6C).
tmp.plot = readRDS("../data/Fig6/6CD_tmp_plot_for_umap_fig6_knn6.rds")
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

print(p)

## Figure 6E ##
#Load libraries
library(skitools)
library(EnhancedVolcano)

#Load file with marker genes.
seu.markers = readRDS("../dataFig6/6E_seu_markers.rds")

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
print(p)
