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
