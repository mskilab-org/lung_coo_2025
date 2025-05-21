# Introduction

This github serves as a repository for the data and code associated to the manuscript entitled: _Passenger mutations link cellular origin and transcriptional identity in human lung adenocarcinomas_, currently in revision.

The repository is composed of 2 folders: `data` and `R`, which contain the files and code required to reproduce the figures within the manuscript. Note that certain figures were further modified via Illustrator to improve legibility and/or presentation, so in some instance the specific styling (e.g.: colors) don't correspond completely between the manuscript figures and the ones shown here.

# Repository Structure

## Code Files

- The `R` folder contains a total of 8 files. Out of these, the 6 that begin with the prefix `fig` contain code for either 1 of the 5 main figures (`fig1.R` to `fig5.R`), or the extended data figures (`figs_extd.R`). Within these files, the code is divided by panel, is commented, and can be run individually for single panels without running the whole file. The required libraries and files for generating each panel (see below) are loaded individually. 
- The `RMDKnit.R` file contains the code used to knit the html notebook with all figures (see below). It's designed to run as is from within it's current location and generate the html in the `data` folder. More details can be found in the initial comments of said file.
- The `COO_Simulations.R` file contains the code used for simulations for extended data figures 6 and 7. The code as is will generate simulations for extended data figure 7, but instructions are included in the code's comments to adjust for generating simulations for extended data figure 6.
- The `scLung_Figures.rmd` file contains the code required to generate an html file with the manuscript figures, via the knitting function ran from `RMDKnit.R`.

## Data Files

- The `data` folder contains all files required for generating the manuscript figures. For specifics of which folder is required for which figure panel, please refer to the code for said panel in the `R` folder.
- This folder also contains the html file with the code and output of all panels, which can be generated again via `RMDKnit.R` (see above).

# Running Instructions

As mentioned above, all code panels can be ran independently of one another. _To run the code as is without having to deal with path issues, please set the `R` as working directory_, as all paths within the code files use the syntaxis to `../data/` when looking for their associated data files to load. However, note that individual paths can be modified within a given panel code with no problem (as long as the data files are loaded the panel code runs correctly).

## Session Information

This code was created and tested using R version 4.3.3. 

- The following base packages are loaded for at least one panel: `parallel`, `grid`, `stats4`, `stats`, `graphics`, `grDevices`, `utils`, `datasets`, `methods`, `base`.

- Non-base packages, with their associated versions, are listed below:

`DiagrammeR_1.0.11`           `data.tree_1.1.0`        `circlize_0.4.16`            `Seurat_5.0.0`               
`SeuratObject_5.0.2`          `sp_2.1-1`               `EnhancedVolcano_1.20.0`      `ggrepel_0.9.4`              
`RColorBrewer_1.1-3`          `ggforce_0.4.2`          `ggpubr_0.6.0`                `forcats_1.0.0`              
`effects_4.2-2`               `carData_3.0-5`          `ggsankey_0.0.99999`          `ggalluvial_0.12.5`          
`readxl_1.4.3`                `wesanderson_0.3.7`      `plyr_1.8.9`                  `skidb_0.1`                  
`MASS_7.3-60`                 `dplyr_1.1.4`            `skitools_0.0.0.9000`         `igraph_2.1.4`               
`reshape2_1.4.4`              `plotly_4.10.4`          `ggplot2_3.5.1`               `stringr_1.5.1`              
`gUtils_0.2.0`                `data.table_1.17.0`      `devtools_2.4.5`              `usethis_3.1.0`              
`htmlwidgets_1.6.4`           `VariantAnnotation_1.46.0`      `Rsamtools_2.16.0`            `Biostrings_2.68.1`          
`XVector_0.42.0`              `SummarizedExperiment_1.32.0`      `Biobase_2.62.0`              `MatrixGenerics_1.14.0`      
`matrixStats_1.5.0`           `ComplexHeatmap_2.18.0`      `GenomicRanges_1.54.1`        `GenomeInfoDb_1.38.1`        
`IRanges_2.36.0`              `S4Vectors_0.40.2`           `BiocGenerics_0.48.1` 

- To check which libraries are loaded for which panels, refer to the code for said panel in the `R` folder.
