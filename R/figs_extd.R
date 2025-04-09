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

## EDF 7 -- Univariate ##
#Load libraries.
library(skitools)
library(MASS)
library(parallel)

#Load GRanges file with genomic covariates and centroids.
juan.genes.gr.nw = readRDS("../data/edf7_genes_gr_LUAD.rds")
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
#IF NEW SIMS ARE TO BE DONE GENERATE AN EMPTY DATA.TABLE FOR SAVING RESULTS BEFORE RUNNING.
  thisRDS = readRDS("../data/testRun_withOverlap_Cov.rds")
  thisRDS = rbind(thisRDS,toSave)
  saveRDS(thisRDS,"../data/testRun_withOverlap_Cov.rds")

}

## EDF 7 -- Multivariate ##
#Load libraries.
library(skitools)
library(MASS)
library(parallel)

#Most of the code is the same as for univariate. Comments are left as is for those parts, and changed when neccesary.
#Load GRanges file with genomic covariates and centroids.
juan.genes.gr.nw = readRDS("../data/edf7_genes_gr_LUAD.rds")
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
#IF NEW SIMS ARE TO BE DONE GENERATE AN EMPTY DATA.TABLE FOR SAVING RESULTS BEFORE RUNNING.  
  thisRDS = readRDS("../data/testRun_withOverlapMultivariatewithCovs.rds")
  thisRDS = rbind(thisRDS,toSave)
  saveRDS(thisRDS,"../data/testRun_withOverlapMultivariatewithCovs.rds")

}


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
# EDF 9A
# ------------------------------------------------------------------------------------------------



pts_coo_id_4B <- readRDS( '../data/pts_coo_id_4B.rds')
pts_coo_id_4C <- copy(pts_coo_id_4B)
pts_coo_id_4C$Lineage_plasticity <- ''
pts_coo_id_4C[Identity == 'Distal Lung' & Origin == 'Distal Lung' ]$Lineage_plasticity <- 'Lineage conserved'
pts_coo_id_4C[Identity == 'Proximal Lung' & Origin == 'Distal Lung' ]$Lineage_plasticity <- 'Lineage plasticity'
pts_coo_id_4C[Identity == 'Distal Lung' & Origin == 'Non-Distal Lung' ]$Lineage_plasticity <- 'Lineage plasticity'
pts_coo_id_4C[Identity == 'Proximal Lung' & Origin == 'Non-Distal Lung' ]$Lineage_plasticity <- 'Lineage conserved'

res.plot =  pts_coo_id_4C[, prop.test(sum(Lineage_plasticity == 'Lineage plasticity'), .N) %>% dflm %>% cbind(nprox = sum(Lineage_plasticity == 'Lineage plasticity'), tot = .N), by = .(TP53_mut = ifelse(TP53_mut, 'TP53 MUT', 'WT'), Origin)][, fracprox := estimate]
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
