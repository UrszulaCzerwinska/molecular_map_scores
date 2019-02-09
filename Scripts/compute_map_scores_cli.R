#!/usr/bin/env Rscript

#+++++++++++++
# command line interface for computing module scores
# expression frame as the first argument, cells in columns,genes in rows
# ranks file as the second argument, cells in rows, ranks in the column (row names need to match), only one column needed
# map.gmt - either a table with pathways in colums or a gmt as the third argument
# opional:
  # path to the folder where the files should be saved
  # extreme_thr - threshold for cells to be selected from ranks dist
  # target - name of the folder with result to be created 
  # th_var_genes - thereshold on the percentage of the top variant genes
  # thr.module - min number of genes in the module 

# you need to provide either only the three first params or all params
  
# result files are saved in the 'target' repo or "res" by default
# R object is savec in res_computue_map_scores.Rdata in the target repo
#++++++++++++++++
  
# this script need to be run from the same repo as the scripts "getPathways.R" and  "compute_map_scores.R"
  
# examples

# only compulsory params
# Rscript compute_map_scores_cli.R  ../Data/preprocessed2/macro_dc_outrem.txt ../Data/ica/ranks_macro.txt ../Data/maps/innate_immune_master_14_05_2017_MODULES_FOR\ STAINING.gmt

# all params 
# Rscript compute_map_scores_cli.R  ../Data/preprocessed2/macro_dc_outrem.txt ../Data/ica/ranks_macro.txt ../Data/maps/innate_immune_master_14_05_2017_MODULES_FOR\ STAINING.gmt 0.25 0.5 ./ myres/ 4


# get args
  
args = commandArgs(trailingOnly = TRUE)
# print(length(args))
# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("provide args", call. = FALSE)
} else if (length(args) > 8) {
  stop("too many args", call. = FALSE)
} else{
  df.file = args[1]
  rnk.file = args[2]
  gmt.file = args[3]
  if (length(args) == 8) {
    extreme_thr = as.numeric(args[4])
    th_var_genes = as.numeric(args[5])
    path = args[6]
    target = args[7]
    thr.module = as.numeric(args[8])
  }
  else{
    extreme_thr = 0.2
    th_var_genes = 0.5
    path="./"
    target = "res/"
    thr.module = 4
  }
 }


# +++++++++++++++++++++++
# used for testinh

# df.file =  "../Data/preprocessed2/macro_dc_outrem.txt"
# rnk.file = "../Data/ica/ranks_macro.txt"
# gmt.file = "../Data/maps/innate_immune_master_14_05_2017_MODULES_FOR STAINING.gmt"
# extreme_thr = as.numeric("0.25")
# th_var_genes = as.numeric("0.5")
# path = "./"
# target = "myres/"
# thr.module = as.numeric("4")

# +++++++++++++++++++++++++

# import files

if (length(grep(".gmt", gmt.file, fixed = T, value = TRUE)) > 0) {
  source("getPathways.R")
  gmt.t <- gmtPathways(gmt.file)
  map.gmt <- gmt_to_table(gmt.t)

} else{
  map.gmt = read.table(
    gmt.file,
    sep = "\t",
    row.names = NULL,
    header = TRUE,
    stringsAsFactors = FALSE
  )
}



df = read.table(
  df.file,
  sep = "\t",
  row.names = 1,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


ranks = read.table(
  rnk.file,
  sep = "\t",
  row.names = 1,
  header = TRUE,
  stringsAsFactors = FALSE
)[,1, drop=FALSE]



#source function

source("compute_map_scores.R")

# run

res <- compute_map_scores(
  df,
  ranks,
  map.gmt,
  extreme_thr = extreme_thr,
  th_var_genes = th_var_genes,
  path = path,
  target = target,
  thr.module = thr.module
)

save(res, file=paste0(path, target,"res_compute_map_scores.RData"))
