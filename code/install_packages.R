#'
#' Installation of the required R packages
#'


# required packages
packages <- c("fda",
              "fdakma",
              "roahd",
              "coda",
              "devtools",
              "fastmatrix",
              "R.matlab",
              "Matrix",
              "progress",
              "pbmcapply",
              "latex2exp",
              "calculus",
              "ReconstPoFD",
              "tidyverse",
              "ggplot2",
              "ggpubr",
              "snowfall",
              "psych",
              "wesanderson")

# new packages
new_packages  <- packages[!(packages %in% installed.packages()[,"Package"])] 

# install required packages
if(length(new_packages))
{
  install.packages(new_packages)
}

library(devtools)
install_github("lidom/ReconstPoFD/ReconstPoFD")
