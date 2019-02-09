# import maps
load('../Data/maps/maps.RData')

# import preprocessed data resulting from "../preprcessing/Macrophages_extraaction.r", "../preprcessing/NK_extraaction.r"
load('../Data/preprocessed2/macro_nk_expression.RData')

# impor weigths after running ICA from '../ica/ICA.r"
load("../Data/ica/ica_res_nk_macro.RData")

# import function
source("compute_map_scores.R")


res_master_map_macro <- compute_map_scores(Macrophages.dc.t,
           A.macro[, 1, drop=FALSE],
           MAP_GMT,
           extreme_thr = 0.25,
           th_var_genes = 0.5,
           path = "../Results/",
           target = "res_macrophages_master_map",
           thr.module = 0)

res_master_map_nk <- compute_map_scores(NK.dc.t,
                                           A.nk[, 1, drop=FALSE],
                                           MAP_GMT,
                                           extreme_thr = 0.5,
                                           th_var_genes = 0.5,
                                           path = "../Results/",
                                           target = "res_nk_master_map",
                                           thr.module = 4)

res_macro_map_macro <- compute_map_scores(Macrophages.dc.t,
                                           A.macro[, 1, drop=FALSE],
                                           MAP_MACRO,
                                           extreme_thr = 0.25,
                                           th_var_genes = 0.5,
                                           path = "../Results/",
                                           target = "res_macrophages_macrophages_map",
                                           thr.module = 0)

res_nk_map_nk <- compute_map_scores(NK.dc.t,
                                        A.nk[, 1, drop=FALSE],
                                        MAP_NK,
                                        extreme_thr = 0.5,
                                        th_var_genes = 0.5,
                                        path = "../Results/",
                                        target = "res_nk_nk_map",
                                        thr.module = 4)
