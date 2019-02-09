#!/bin/bash

Rscript compute_map_scores_cli.R  ../Data/preprocessed2/macro_dc_outrem.txt ../Data/ica/ranks_macro.txt ../Data/maps/innate_immune_master_14_05_2017_MODULES_FOR\ STAINING.gmt 0.25 0.5 ./ myres/ 4

Rscript compute_map_scores_cli.R  ../Data/preprocessed2/macro_dc_outrem.txt ../Data/ica/ranks_macro.txt ../Data/maps/innate_immune_master_14_05_2017_MODULES_FOR\ STAINING.gmt
