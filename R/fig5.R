## Fig 5B Left Panel ##
#Load libraries
library(skitools)
library(parallel)
library(MASS)

#Load csv with relative risk information per centroid for every sample.
rel_risk = read.table("../Fig5B_luad_rel_risk.csv",sep=",",header=T,row.names=1)

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
wcm = readRDS("../data/Fig5B_cov_res.rds")
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
dt.plot = readRDS("../data/Fig5B_cov_res.rds")
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
  }, file= "Fig_5BRight.pdf",cex = c(0.95,1.25))
