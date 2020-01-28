#!/usr/bin/env Rscript

list.of.packages <- c("ggplot2", "reshape", "grid", "gridExtra", "gtable", "scales", "argparse", "assertthat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = T, repos="http://cran.ma.imperial.ac.uk/")
