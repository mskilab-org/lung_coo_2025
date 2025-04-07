## Figure 1C ##
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
proximal = setdiff(colnames(epicells),distal_rare)
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
dft2 <- LUAD_query %>%
  group_by(Patient, cat1) %>%
  summarise(count = n()) %>%             # Count the rows per group
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
