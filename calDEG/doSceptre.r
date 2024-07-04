library("Seurat")
library("SCREE")
library("anndata")
library("reticulate")
library("sceptre")
path_to_python <- "/home/wzt/anaconda3/bin/python"
use_python(path_to_python)

Args = commandArgs(T)

dirName = Args[1]
species = Args[2]
setwd(dirName)

data <- read_h5ad("filterRaw.h5ad")
response_matrix_lowmoi = t(as.matrix(data$to_df()))  ### 用count数据*****

grna_matrix_lowmoi = read.table('sgRNA_matrix.txt', sep='\t', header = T, row.names = 1, check.names = F)
grna_matrix_lowmoi = as.matrix(grna_matrix_lowmoi)
response_matrix_lowmoi = response_matrix_lowmoi[, colnames(grna_matrix_lowmoi)]  ### 保证顺序一致

grna_group_data_frame_lowmoi = read.table('grna_target_data_frame.tsv', sep='\t', header=T)

sceptre_object_lowmoi <- import_data(
  response_matrix = response_matrix_lowmoi, 
  grna_matrix = grna_matrix_lowmoi, 
  grna_target_data_frame = grna_group_data_frame_lowmoi,
  moi = "low"
)
positive_control_pairs_lowmoi <- construct_positive_control_pairs(
  sceptre_object = sceptre_object_lowmoi
)
discovery_pairs_lowmoi <- construct_trans_pairs(
  sceptre_object = sceptre_object_lowmoi,
  positive_control_pairs = positive_control_pairs_lowmoi,
  pairs_to_exclude = "pc_pairs"
)

set_analysis_parameters(discovery_pairs_lowmoi, positive_control_pairs_lowmoi)







sceptre_object <- sceptre_object_lowmoi  |> # |> is R's base pipe, similar to %>%
  set_analysis_parameters(discovery_pairs_lowmoi, positive_control_pairs_lowmoi) |>
  run_calibration_check() |>
  run_power_check() |>
  run_discovery_analysis()








# tmp = data1$obs
# covariate_data_frame_lowmoi = tmp[colnames(grna_matrix_lowmoi), c('total_counts', 'n_genes_by_counts', 'pct_counts_mt')]
# colnames(covariate_data_frame_lowmoi) = c('response_n_umis', 'response_n_nonzero', 'p_mito')


# response_grna_group_pairs <- generate_all_pairs(response_matrix_lowmoi,
#                                                 grna_group_data_frame_lowmoi)

# if (sum(covariate_data_frame_lowmoi$p_mito) == 0){
# formula_object <- formula(~log(response_n_umis) +  log(response_n_nonzero))
# }else{
# formula_object <- formula(~log(response_n_umis) +  log(response_n_nonzero) + p_mito)
# }


# calibration_result <- run_sceptre_lowmoi(
#   response_matrix = response_matrix_lowmoi,
#   grna_matrix = grna_matrix_lowmoi,
#   covariate_data_frame = covariate_data_frame_lowmoi,
#   grna_group_data_frame = grna_group_data_frame_lowmoi,
#   formula_object = formula_object,
#   response_grna_group_pairs = response_grna_group_pairs,
#   calibration_check = TRUE
# )



# discovery_result <- run_sceptre_lowmoi(
#   response_matrix = response_matrix_lowmoi,
#   grna_matrix = grna_matrix_lowmoi,
#   covariate_data_frame = covariate_data_frame_lowmoi,
#   grna_group_data_frame = grna_group_data_frame_lowmoi,
#   formula_object = formula_object,
#   response_grna_group_pairs = response_grna_group_pairs,
#   calibration_check = FALSE
# )

# write.table(discovery_result, 'sceptre/rawResult.tsv', sep='\t', col.names = T, quote = F, row.names = F)
