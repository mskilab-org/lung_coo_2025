## Short script for knitting the RMD file into an html. By default this code is ran assuming that the RMD file is in the data folder, but this can be changed below.
## The variables of relevance used in this script are as follows:
#' @param workingDir is the main working folder, which contains both the folder with R code (R) and the one with data files (data). It is assumed by default that the script is ran from it's location within the repository directory.
#' @param input is the location of the RMD to be knitted.
#' @param output_file is the location in which the html will be saved (by default the data folder) and it's name (by default scLungFigures.html).
#' @param knit_root_dit is the directory used during the knitting process (by default the data folder).
#' @param params is a list with the parameters to be passed to the RMD. There're only 2 parameters: set_title, which indicates the title to be given to the html, and fileLoc, which indicates where to find the data files used in the generation of the figures (by default the data directory).
workingDir = dirname(getwd())
rmarkdown::render(
        input = normalizePath(file.path(workingDir,"data","scLungFigures.rmd")),
        output_format = "html_document",
        output_file = normalizePath(file.path(workingDir,"data","scLungFigures.html")),
        knit_root_dir = normalizePath(file.path(workingDir,"data"),
        params = list(set_title = "Paper Figures",fileLoc = normalizePath(file.path(workingDir,"data")),quiet = FALSE)
