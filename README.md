# Introduction

This github serves as a repository for the data and code associated to the manuscript entitled: _Passenger mutations link cellular origin and transcriptional identity in human lung adenocarcinomas_, currently in revision.

The repository is composed of 2 folders: `data` and `R`, which contain the files and code required to reproduce the figures within the manuscript. Note that certain figures were further modified via Illustrator to improve legibility and/or presentation, so in some instance the specific styling (e.g.: colors) don't correspond completely between the manuscript figures and the ones shown here.

# Repository Structure

## Code Files

- The `R` folder contains a total of 8 files. Out of these, the 6 that begin with the prefix `fig` contain code for either 1 of the 5 main figures (`fig1.R` to `fig5.R`), or the extended data figures (`figs_extd.R`). Within these files, the code is divided by panel, is commented, and can be run individually for single panels without running the whole file. The required libraries and files for generating each panel (see below) are loaded individually. 
- The `RMDKnit.R` file contains the code used to generate the html notebook with all figures (see below). It's designed to run as is from within it's current location and generate the html in the `data` folder. More details can be found in the initial comments of said file.
- The `COO_Simulations.R` file contains the code used for simulations for extended data figures 6 and 7. The code as is will generate simulations for extended data figure 7, but instructions are included in the code's comments to adjust for generating simulations for extended data figure 6.

