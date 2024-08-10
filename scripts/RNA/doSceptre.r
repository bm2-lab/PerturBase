library("Seurat")
library("SCREE")
library("anndata")
library("reticulate")
library("sceptre")
Args = commandArgs(T)

species = Args[1]
path_to_python = Args[2] 
use_python(path_to_python)
myUtil_r = Args[3] 
# source('/NFS_home/NFS_home_2/wzt/project/HC-CrisprScreen/myUtil.r')
# source(myUtil_r)



data <- read_h5ad("filterRaw.h5ad")
response_matrix_lowmoi = t(as.matrix(data$to_df()))  ### count data

grna_matrix_lowmoi = read.table('barcode_sceptre.txt', sep='\t', header = T, row.names = 1, check.names = F)
grna_matrix_lowmoi = as.matrix(grna_matrix_lowmoi)
response_matrix_lowmoi = response_matrix_lowmoi[, colnames(grna_matrix_lowmoi)]

grna_group_data_frame_lowmoi = read.table('grnaGroup.tsv', sep='\t', header=T)

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

#set_analysis_parameters(discovery_pairs_lowmoi, positive_control_pairs_lowmoi)



sceptre_object <- sceptre_object_lowmoi  |> # |> is R's base pipe, similar to %>%
  set_analysis_parameters(discovery_pairs_lowmoi, positive_control_pairs_lowmoi) |>
  run_calibration_check() |>
  run_power_check() |>
  run_discovery_analysis(parallel = TRUE)
result <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)
write.table(result, file = "sceptre/rawResult.tsv", sep = "\t",row.names = FALSE, col.names = TRUE)
