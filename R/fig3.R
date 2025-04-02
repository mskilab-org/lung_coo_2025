library(skitools)

mat.in.2 <- readRDS('../data/mat.in.2.rds')
# fig 3A
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


