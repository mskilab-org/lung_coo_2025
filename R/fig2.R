## FIGURE 2B ##
library(skitools)
library(MASS)
library(parallel)
res = readRDS("../data/Fig2B.rds")

#Set estimates value for each celltype below fdr threshold (0.1) to null value of 1. This effectively discards said celltypes from COO calling. 
dt = dcast.data.table(res[grep('cov', name), ][, estimate2 := ifelse(fdr<0.1, estimate, 1)], tt ~ celltype, value.var = 'estimate2')

#Take the estimates and add back GLM information required for plotting (confidence interval values, celltype, fdr, etc.).
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
dt.plot$tumor.type = "LUAD"
dt.plot$tumor.type = toupper(dt.plot$tumor.type)
dt.plot$combo.label = paste0(dt.plot$tumor.type,"_",dt.plot$celltype)
dt.plot$Cell_Type = dt.plot$celltype


#Center relative risk and confidence interval estimate values at 1 for plot.
dt.plot$estimate = dt.plot$estimate - 1
dt.plot$ci.lower = dt.plot$ci.lower - 1
dt.plot$ci.upper = dt.plot$ci.upper - 1

# Relative risk barplots. The code allows for filtering for specific cancer type, snv columns, and celltypes via in.pattern.tt, in.pattern.snv, and in.ct.sub respectively. Here we add the values needed for plotting all GLM outputs defined above with no restriction. Note that multiple cancer subtype grpahs can be done by adding more values to in.pattern.tt. 
in.pattern.tt = c('LUAD')
in.pattern.snv = "snv.count"
in.ct.sub = epicells
ppdf(
  for (ik in in.pattern.tt)
  {
    message('\nPlotting mutational density for: ', ik, ' samples')
#Filter for any cancer type, celltype, or snv column if needed. 
    in.pattern = in.pattern.snv
    in.count.sub = grep(paste0(in.pattern),dt.plot$tt, value = TRUE)
    dt.plot.sub = dt.plot[tumor.type %in% ik & celltype %in% in.ct.sub & tt %in% in.count.sub,]
#Order results by relative risk (lower to higher).     
    setorder(dt.plot.sub,celltype)
    in.title = unique(gsub(" - .*","",dt.plot.sub$combination.formal))
# Set upper and lower limits for plot (x axis).
    in.axis.breaks = 0.05
    in.axis.min = -0.25
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
  }, file= "Fig2B_LUAD.pdf",cex = c(0.95,1.25))

## FIGURE 2C ##
library(skitools)
library(MASS)
library(parallel)
res = readRDS("../data/Fig2C.rds")

#Now set estimates value for each celltype below fdr threshold (0.1) to null value of 1. This effectively discards said celltypes from COO calling. 
dt = dcast.data.table(res[grep('cov', name), ][, estimate2 := ifelse(fdr<0.1, estimate, 1)], tt ~ celltype, value.var = 'estimate2')

#Take the estimates and add back GLM information required for plotting (confidence interval values, celltype, fdr, etc.).
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
dt.plot$tumor.type = "LUSC"
dt.plot$tumor.type = toupper(dt.plot$tumor.type)
dt.plot$combo.label = paste0(dt.plot$tumor.type,"_",dt.plot$celltype)
dt.plot$Cell_Type = dt.plot$celltype

#Center relative risk and confidence interval estimate values at 1 for plot.
dt.plot$estimate = dt.plot$estimate - 1
dt.plot$ci.lower = dt.plot$ci.lower - 1
dt.plot$ci.upper = dt.plot$ci.upper - 1

# Relative risk barplots. The code allows for filtering for specific cancer type, snv columns, and celltypes via in.pattern.tt, in.pattern.snv, and in.ct.sub respectively. Here we add the values needed for plotting all GLM outputs defined above with no restriction. Note that multiple cancer subtype grpahs can be done by adding more values to in.pattern.tt. 
in.pattern.tt = c('LUSC')
in.pattern.snv = "snv.count"
ppdf(
  for (ik in in.pattern.tt)
  {
    message('\nPlotting mutational density for: ', ik, ' samples')
#Filter for any cancer type, celltype, or snv column if needed.     
    in.ct.sub = epicells
    in.pattern = in.pattern.snv
    in.count.sub = grep(paste0(in.pattern),dt.plot$tt, value = TRUE)
    dt.plot.sub = dt.plot[tumor.type %in% ik & celltype %in% in.ct.sub & tt %in% in.count.sub,]
#Order results by relative risk (lower to higher).         
    setorder(dt.plot.sub,celltype)
    print(dt.plot.sub)
    in.title = unique(gsub(" - .*","",dt.plot.sub$combination.formal))
# Set upper and lower limits for plot (x axis).
    in.axis.breaks = 0.05
    in.axis.min = -0.25
    in.axis.max = 10.05
#Construct ggplot2 object and plot.        
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
  }, file= "Fig2C_LUSC.pdf",cex = c(0.95,1.25))
    
