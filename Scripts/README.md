# Intro

This scripts aim to compute activity scores in groups of *single cells* using scRNAseq genes expression in functional modules of molecular map created in Sysbio group at Institute Curie.

# How it works

## Description

The function **compute_map_scores** takes the gene expression, ranks vector and molecular map modules table and computes module scores.
Cells are divided into groups based on the quantile threshold of the ranks distribution. For each map module genes present in the module are selected from gene expression table. The mean of the x% most variant genes is computed for each module

## Params 

- df - gene expression file, tab deliminated, genes in rows and calls in columns
- ranks - ranks file forcell grouping, here IC1 projections, cells in rownames,one column only needed
- map - table with modules as columns, to transform gmt file into a gmt table use getPathways.R script
- extreme_thr - threshold used to define groups of cells, i.e. 0.25 => 25% of top and 25% of down ranked genes
- th_var_genes - threshold applied to select the percentage of the top variant genes, 0.5 by default
- path - the path where the results will be written
- target - name of the directory with results that will be created
- thr.module - min number of genes in the module to be consided, default 4

## Output

### In R

- genes_in_modules - list of matrices with genes selected by map module
- scores - matrix of computed scores
- t_test - t-test results
- neg_cells - list of cells in the group of negative tail of ranks
- pos_cells - list of cells in the positive tail of ranks

### Saved in the filesystem 

As a result files are saved in the path/target directory:
  aggregated scores :
        - folder with the heatmap of scores with t-test stars
        - scores matrix used map for staining
  heatmap_extreme_cells :
         heatmaps of only selected cell in groups with threshold applied on genes
  heatmap_most_var_genes :
        heatmap of all cells with threshold applied on genes
and files :
   - _modules_size.txt - number of genes present in the data and in the module
   - _t_test.txt - t-test results
   - most_neg_pos_genes.txt - for each module three genes with highest (up) and three genes lowest expression (down) are listed

# examples of use

see/execute

**running_scores_function.R**

Can be run directly from command line with input of .txt and .gmt files

**compute_map_scores_cli.R**

and 

see/execute 

**command_line_example.sh**

