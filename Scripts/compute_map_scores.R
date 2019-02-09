

#++++++++++++++++++++++++++++++++++
# Compute map scores
#
# This function takes the gene expression, ranks vector and molecular map modules table and computes module scores.
# Cells are divided into groups based on the quantile threshold of the ranks distribution
# For each map module genes present in the module are selected from gene expression table
# The mean of the x% most variant genes is computed for each module
# As a result files are saved in the path/target directory:
#   aggregated scores :
#         - folder with the heatmap of scores with t-test stars
#         - scores matrix used map for staining
#   heatmap_extreme_cells :
#         heatmaps of only selected cell in groups with threshold applied on genes
#   heatmap_most_var_genes :
#         heatmap of all cells with threshold applied on genes
# and files :
#   - _modules_size.txt - number of genes present in the data and in the module
#   - _t_test.txt - t-test results
#   - most_neg_pos_genes.txt - for each module three genes with highest (up) and three genes lowest expression (down) are listed

# @param df - gene expression file, tab deliminated, genes in rows and calls in columns
# @param ranks - ranks file forcell grouping, here IC1 projections, cells in rownames,one column only needed
# @param map - table with modules as columns, to transform gmt file into a gmt table use getPathways.R script
# @param extreme_thr - threshold used to define groups of cells, i.e. 0.25 => 25% of top and 25% of down ranked genes
# @param th_var_genes - threshold applied to select the percentage of the top variant genes, 0.5 by default
# @param path - the path where the results will be written
# @param target - name of the directory with results that will be created
# @param thr.module - min number of genes in the module to be consided, default 4
#
# @return
#  genes_in_modules - list of matrices with genes selected by map module
#  scores - matrix of computed scores
# t_test - t-test results
# neg_cells - list of cells in the group of negative tail of ranks
# pos_cells - list of cells in the positive tail of ranks
# @export
#
# @examples
# for example with parameters used in the article see
# running_scores_function.R
#
#+++++++++++++++++++++++++++++++++++
# NOTE
# this function can be run directly from command line with input of .txt and .gmt files
# see:
# compute_map_scores_cli.R

compute_map_scores <-
  function(df,
           ranks,
           map,
           extreme_thr = 0.2,
           th_var_genes = 0.5,
           path = "./",
           target = "res",
           thr.module = 4) {
    require(pheatmap)
    require(ggplot2)
    dir.create(paste(path, target, sep = ""), showWarnings = FALSE)
    diff_list <- NULL
    diff_list_short <- NULL
    if (nrow(ranks) != ncol(df)) {
      stop("ranks need to have the same length as number of rows in df")
    } else if (!identical(rownames(ranks), colnames(df))) {
      stop("ranks and df need to have the same rownames")
    }
    
    get_samples_from_quantiles <- function(ranks, quantiles) {
      #quantile(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7, ...)
      quantiles.dataset <-
        quantile(
          ranks[, 1],
          probs = quantiles,
          na.rm = FALSE,
          names = TRUE,
          type = 7
        )
      
      
      extreme.dataset.neg  <-
        ranks[ranks < quantiles.dataset [1], , drop = FALSE]
      extreme.dataset.pos  <-
        ranks[ranks[, 1] > quantiles.dataset [2], , drop = FALSE]
      
      
      return(list(extreme.dataset.neg, extreme.dataset.pos))
    }
    
    
    module_size <-
      data.frame(colnames(map), rep(NA, length(colnames(map))))
    quant = get_samples_from_quantiles(ranks, c(extreme_thr, 1 - extreme_thr))
    ic1neg <- data.frame(quant[1])
    ic1plus <- data.frame(quant[2])
    
    ic1neg$group <-
      rep("IC1neg", nrow(ic1neg))
    
    ic1plus$group <-
      rep("IC1plus", nrow(ic1plus))
    
    ic1 <- rbind(ic1neg, ic1plus)
    # df = t(df)
    mean.class_map_innate = apply(map, 2, function(x) {
      na.omit(df[which(row.names(df) %in% as.character(x)),])
    })
    
    
    
    sum.df <-  data.frame(
      MODULE = character(),
      IC1neg = numeric(),
      IC1pos = numeric(),
      stringsAsFactors = FALSE
    )
    
    sample.mean.df <-
      data.frame(matrix(
        ,
        nrow = ncol(map)  ,
        ncol = length(ic1plus$group) + length(ic1neg$group)
      ),
      stringsAsFactors = FALSE)
    
    row.names(sample.mean.df) <- colnames(map)
    colnames(sample.mean.df) <-
      c(rownames(ic1neg), rownames(ic1plus))
    
    RowVar <- function (x, na.rm = TRUE)
    {
      sqr = function(x)
        x * x
      n = rowSums(!is.na(x))
      n[n <= 1] = NA
      return(rowSums(sqr(x - rowMeans(x, na.rm = na.rm)), na.rm = na.rm) /
               (n - 1))
    }
    compute_t_test_for_many <- function (df, array1, array2) {
      pvalue = NULL # Empty list for the p-values
      tstat = NULL # Empty list of the t test statistics
      
      for (i in 1:nrow(df)) {
        # For each gene :
        x = df[i, array1] # condition1 of gene number i
        y = df[i, array2] # condition2 of gene number i
        
        # Compute t-test between the two conditions
        t = t.test(x, y)
        
        # Put the current p-value in the pvalues list
        pvalue[i] = t$p.value
        
        # Put the current t-statistic in the tstats list
        tstat[i] = t$statistic
        
      }
      
      pvalue.df = data.frame(pvalue)
      tstat.df = data.frame(tstat)
      
      row.names(pvalue.df) = row.names(df)
      row.names(tstat.df)[i] = row.names(df)[i]
      
      return(cbind(pvalue.df, tstat))
    }
    
    add_significance <- function(df, col = "pvalue") {
      df$significance <- rep("", nrow(df))
      for (i in 1:nrow(df)) {
        if (df[[col]][i] < 0.001) {
          df$significance[i] <- "***"
        } else if (df[[col]][i] > 0.001 &  df[[col]][i] < 0.01) {
          df$significance[i] <- "**"
        } else if (df[[col]][i] > 0.01 &  df[[col]][i] < 0.05) {
          df$significance[i] <- "*"
        } else if (df[[col]][i] > 0.05 &  df[[col]][i] < 0.1) {
          df$significance[i] <- "."
        } else {
          
        }
      }
      df$legend = paste(row.names(df), df$significance, sep = "  ")
      return(df)
    }
    #
    pheatmap_palette_zero_center <- function(dataframe) {
      paletteLength <- 50
      myColor <-
        colorRampPalette(c("green", "white", "red"))(paletteLength)
      # use floor and ceiling to deal with even/odd length pallettelengths
      myBreaks <-
        c(
          seq(min(dataframe), 0, length.out = ceiling(paletteLength / 2) + 1),
          seq(
            max(dataframe) / paletteLength,
            max(dataframe),
            length.out = floor(paletteLength / 2)
          )
        )
      
      return(list(myColor, myBreaks))
    }
    for (module in 1:length(mean.class_map_innate)) {
      if (is.null(nrow(mean.class_map_innate[[module]]))) {
        module_size[module, 2] <-  1
        
      } else if (nrow(mean.class_map_innate[[module]]) <= thr.module) {
        module_size[module, 2] <- nrow(mean.class_map_innate[[module]])
        
      } else{
        #get 50% of most variable genes per module
        mean.class_map_innate[[module]] <-
          data.frame(mean.class_map_innate[[module]], var = RowVar(mean.class_map_innate[[module]]))
        mean.class_map_innate[[module]] <-
          mean.class_map_innate[[module]][order(-mean.class_map_innate[[module]][, ncol(df) +
                                                                                   1]), ]
        mean.class_map_innate[[module]] <-
          mean.class_map_innate[[module]][mean.class_map_innate[[module]][, ncol(df) +
                                                                            1] >
                                            quantile(mean.class_map_innate[[module]][, ncol(df) +
                                                                                       1], th_var_genes) ,]
        main <-
          paste(names(mean.class_map_innate)[module], "top_var_genes", sep = "_")
        
        dir.create(paste(path, target, "/heatmap_most_var_genes/", sep = ""),
                   showWarnings = FALSE)
        #pheatmap for all cells per module
        pheatmap(
          mean.class_map_innate[[module]][, 1:ncol(df)],
          main = paste(
            "heatmap of top",
            th_var_genes * 100,
            "% of most variant genes, all cells in",
            names(mean.class_map_innate)[module],
            sep = " "
          ),
          width = 10,
          height = 12,
          filename = paste(
            path,
            paste(path, target, "/heatmap_most_var_genes/", sep = ""),
            main,
            ".pdf" ,
            sep = ""
          )
        )
        
        module_size[module, 2] <-
          nrow(mean.class_map_innate[[module]][, rownames(ic1)])
        
        #pheatmap for extremes only per module
        annotation_cols <- ic1[, 2, drop = FALSE]
        
        
        mean_gene_neg <-
          rowMeans(mean.class_map_innate[[module]][, rownames(ic1neg)])
        mean_gene_pos <-
          rowMeans(mean.class_map_innate[[module]][, rownames(ic1plus)])
        
        diff_neg_pos <- mean_gene_neg - mean_gene_pos
        diff_list[[module]] <-
          diff_neg_pos[order(mean_gene_neg - mean_gene_pos)]
        
        if (tail(diff_list[[module]], 1) > 0 &
            head(diff_list[[module]], 1) < 0) {
          diff_list_short[[module]] <-
            c(head(diff_neg_pos[order(mean_gene_neg - mean_gene_pos)], 3), tail(diff_neg_pos[order(mean_gene_neg - mean_gene_pos)], 3))
          diff_genes <- paste(
            names(mean.class_map_innate)[module],
            paste("UP_PLUS", paste(
              names(diff_list_short[[module]][1:3]), collapse = ", "
            ), sep = " : "),
            paste("UP_NEG", paste(
              names(diff_list_short[[module]][4:6]), collapse = ", "
            ), sep = " : "),
            sep = " : "
          )
          
        } else if (tail(diff_list[[module]], 1) > 0) {
          diff_list_short[[module]] <-
            tail(diff_neg_pos[order(mean_gene_neg - mean_gene_pos)], 3)
          diff_genes <- paste(names(mean.class_map_innate)[module],
                              paste("UP_NEG", paste(
                                names(diff_list_short[[module]][4:6]), collapse = ", "
                              ), sep = " : "),
                              sep = " : ")
          
        } else {
          diff_list_short[[module]] <-
            head(diff_neg_pos[order(mean_gene_neg - mean_gene_pos)], 3)
          diff_genes <- paste(names(mean.class_map_innate)[module],
                              paste("UP_PLUS", paste(
                                names(diff_list_short[[module]]), collapse = ", "
                              ), sep = " : "),
                              sep = " : ")
          
        }
        
        
        diff_genes <- paste(
          names(mean.class_map_innate)[module],
          paste("UP_PLUS", paste(
            names(diff_list_short[[module]][1:3]), collapse = ", "
          ), sep = " : "),
          paste("UP_NEG", paste(
            names(diff_list_short[[module]][4:6]), collapse = ", "
          ), sep = " : "),
          sep = " : "
        )
        
        write.table(
          diff_genes,
          file = paste(
            path,
            target,
            "/",
            paste("most_neg_pos_genes.txt" , sep = ""),
            sep = ""
          ),
          append = TRUE,
          quote = FALSE,
          sep = " ",
          row.names = FALSE,
          col.names = FALSE
        )
        
        dir.create(paste(path, target, "/heatmap_extreme_cells/", sep = ""),
                   showWarnings = FALSE)
        
        pheatmap(
          mean.class_map_innate[[module]][, rownames(ic1)],
          main = paste(
            "heatmap of top",
            th_var_genes * 100,
            "% of most variant genes,  IC1 extremes",
            extreme_thr,
            " in ",
            names(mean.class_map_innate)[module],
            sep = " "
          ),
          cluster_cols = FALSE,
          width = 10,
          height = 12,
          annotation = annotation_cols,
          filename = paste(
            paste(path, target, "/heatmap_extreme_cells/", sep = ""),
            paste(
              names(mean.class_map_innate)[module],
              "var_genes_extremes",
              sep = "_"
            ),
            ".pdf" ,
            sep = ""
          )
        )
        
        top_genes <-
          paste(names(mean.class_map_innate)[module],
                paste(rownames(mean.class_map_innate[[module]])[1:3], collapse = ", "),
                sep = " : ")
        
        sum.df[module, 1] <- names(mean.class_map_innate)[module]
        sum.df[module, 2] <-
          mean(rowMeans(mean.class_map_innate[[module]][, rownames(ic1neg)]))
        sum.df[module, 3] <-
          mean(rowMeans(mean.class_map_innate[[module]][, rownames(ic1plus)]))
        
        sample.mean.df[module,] <-
          cbind(t(colMeans(mean.class_map_innate[[module]][, rownames(ic1neg)])),
                t(colMeans(mean.class_map_innate[[module]][, rownames(ic1plus)])))
        
        
        
      }
    }
    
    t.test.res <-
      compute_t_test_for_many(na.omit(sample.mean.df),
                              rownames(ic1neg),
                              rownames(ic1plus))
    t.test.res.stars <- add_significance(t.test.res)
    
    
    write.table(
      t.test.res.stars,
      file = paste(path,  target, "/" , deparse(substitute(map)), "_t_test.txt", sep =
                     ""),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    
    colnames(module_size) = c("MODULE", "N_GENES")
    
    write.table(
      module_size,
      file = paste(
        path,
        target,
        "/",
        deparse(substitute(map)),
        "_modules_size.txt",
        sep = ""
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    
    sum.df_rows <- na.omit(sum.df)
    row.names(sum.df_rows) <- t.test.res.stars$legend
    sum.df_rows$MODULE <- NULL
    col.breaks.3 = pheatmap_palette_zero_center(sum.df_rows)
    
    dir.create(paste(path, target, "/aggregated_scores/", sep = ""),
               showWarnings = FALSE)
    
    pheatmap(
      sum.df_rows,
      color = col.breaks.3[[1]],
      breaks = col.breaks.3[[2]],
      cluster_cols = FALSE,
      main = paste(
        "heatmap_of_ica_extrems",
        extreme_thr,
        "gene_means",
        deparse(substitute(map)),
        sep = "_"
      ),
      filename = paste0(
        paste(path, target, "/aggregated_scores/", sep = ""),
        paste(
          "heatmap_of_ica_extrems",
          extreme_thr,
          "gene_means",
          deparse(substitute(map)),
          sep = "_"
        ),
        ".pdf"
      )
    )
    
    
    sum.df <-
      data.frame(sum.df, DIFF_NEG_POS = sum.df[, 2] - sum.df[, 3])
    write.table(
      sum.df,
      file = paste(
        paste(path, target, "/aggregated_scores/", sep = ""),
        deparse(substitute(map)),
        "_heatmap_modules.txt",
        sep = ""
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    scores <- sum.df
    genes_in_modules <- mean.class_map_innate
    return(
      list(
        genes_in_modules = genes_in_modules,
        scores= scores,
        t_test= t.test.res.stars,
        neg_cells = rownames(ic1neg),
        pos_cells = rownames(ic1plus)
      )
    )
  }
