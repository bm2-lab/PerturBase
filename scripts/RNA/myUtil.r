Input_preprocess_10X<-function(directory){
  require(hash)
  require(stringr)
  perturb_seq <- Read10X(directory)
  perturb_seq <- as.matrix(perturb_seq)
  cbc_gbc <- read.table(paste(directory, "cbc_gbc_dict.tsv",sep = "/"), stringsAsFactors = FALSE)
  cbc_gbc <- unique(cbc_gbc)
  data_preprocess <- function(perturb_seq, cbc_gbc) {
    cell_KO_hash = hash()
    for (i in 1:nrow(cbc_gbc)) {
      cbc = cbc_gbc[i, 1]
      gbc = cbc_gbc[i, 2]
      if (has.key(cbc, cell_KO_hash)) {
        cell_KO_hash[cbc] = paste(cell_KO_hash[[cbc]],gbc, sep = ",")
      }
      else {
        cell_KO_hash[cbc] = gbc
      }
    }
    perturb_information = c()
    j = 1
    k = 1
    nogbc_col = c()
    for (i in 1:ncol(perturb_seq)) {
      if (!is.null(cell_KO_hash[[colnames(perturb_seq)[i]]])) {
        perturb_information[k] = cell_KO_hash[[colnames(perturb_seq)[i]]]
        k = k + 1
      }
      else {
        nogbc_col[j] <- i
        j = j + 1
      }
    }

    if (!is.null(nogbc_col)){
        perturb_seq <- perturb_seq[, -nogbc_col]
    }
    names(perturb_information) = colnames(perturb_seq)
    for (i in 1:length(perturb_information)) {
      sample_info_arr <- unlist(strsplit(perturb_information[i],","))
      if (length(sample_info_arr) > 1) {
        sortedMultiKO<-c()
        sample_info_arr <- sort(sample_info_arr)
        if(sample_info_arr[1]!="CTRL"){
          sortedMultiKO=sample_info_arr[1]
          for (j in 2:length(sample_info_arr)) {
            if(sample_info_arr[j]!="CTRL"){
              sortedMultiKO = paste(sortedMultiKO, sample_info_arr[j],sep = ",")
            }
          }
        }
        else{
          sortedMultiKO=sample_info_arr[2]
          if(length(sample_info_arr)>=3){
            for (j in 3:length(sample_info_arr)) {
              if(sample_info_arr[j]!="CTRL"){
                sortedMultiKO = paste(sortedMultiKO, sample_info_arr[j],sep = ",")
              }
            }
          }
        }
        perturb_information[i] = sortedMultiKO
      }
    }
    perturb_seq <- perturb_seq[!str_detect(row.names(perturb_seq),"^MRP"), ]
    perturb_seq <- perturb_seq[!str_detect(row.names(perturb_seq),"^RP"), ]
    return(list("perturb_data" = perturb_seq, "perturb_information" = perturb_information))
  }
  perturb_seq_list <- data_preprocess(perturb_seq, cbc_gbc)
  perturb_list = list("expression_profile" = perturb_seq_list$perturb_data, "perturb_information" = perturb_seq_list$perturb_information)
  return(perturb_list)
}




scmageck_rra = function (BARCODE, RDS, GENE, RRAPATH = NULL, LABEL = NULL, NEGCTRL = NULL, 
    SIGNATURE = NULL, KEEPTMP = FALSE, PATHWAY = FALSE, SAVEPATH = "./") 
{
    if (is.null(RRAPATH)) {
        RRAPATH = system.file("bin", "RRA", package = "scMAGeCK")
    }
    message("Checking RRA...")
    if (!file.exists(RRAPATH)) {
        if (system("RRA", ignore.stdout = TRUE, ignore.stderr = TRUE) != 
            0) {
            message("RRA does not exist! Please check RRA executable file path")
            return(NULL)
        }
        else {
            RRAPATH = NULL
        }
    }
    if (!is.null(LABEL)) {
        data_label = LABEL
    }
    else {
        data_label = "sample1"
    }
    bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE, sep='\t')
    if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 
        3) {
        stop("cell, barcode, or gene column names not found in barcode file.")
    }
    keep_tmp = KEEPTMP
    message(paste("keep_tmp:", keep_tmp))
    if (!is.null(SIGNATURE)) {
        message(paste("ispathway: TRUE"))
        message(paste("run_signature: TRUE"))
    }
    else {
        ispathway = PATHWAY
        message(paste("ispathway:", ispathway))
    }
    if (!is.null(NEGCTRL)) {
        negctrl_gene = NEGCTRL
    }
    else {
        negctrl_gene = NULL
    }
    guide_count = table(bc_dox$cell)
    ncnt = table(table(bc_dox$cell))
    dupsq = bc_dox[duplicated(bc_dox$cell), 1]
    bc_dox_uq = bc_dox[!bc_dox[, 1] %in% dupsq, ]
    rownames(bc_dox_uq) = bc_dox_uq[, 1]
    message(paste("Total barcode records:", nrow(bc_dox)))
    message(paste("Unique barcode records:", nrow(bc_dox_uq)))
    if (is.character(RDS)) {
        message(paste("Reading RDS file:", RDS))
        targetobj = readRDS(RDS)
    }
    else {
        targetobj = RDS
    }
    nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
    if (nmatch == 0) {
        message("Cell names in expression matrix and barcode file do not match. Try to remove possible trailing \"-1\"s...")
        if (length(grep("-\\d$", bc_dox[, 1])) > 0) {
            bc_dox[, 1] = sub("-\\d$", "", bc_dox[, 1])
            bc_dox_uq[, 1] = sub("-\\d$", "", bc_dox_uq[, 1])
            rownames(bc_dox_uq) = bc_dox_uq[, 1]
        }
        nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
        if (nmatch == 0) {
            stop("No cell names match in expression matrix and barcode file.")
        }
    }
    if ("scale.data" %in% names(attributes(targetobj))) {
        scalef = targetobj@scale.data
    }
    else {
        scalef = GetAssayData(object = targetobj, slot = "scale.data")
    }
    if (!is.null(SIGNATURE)) {
        newdir <- paste0("GENE_SET")
        dir.create(file.path(SAVEPATH, newdir))
        cwd <- getwd()
        setwd(file.path(SAVEPATH, newdir))
        gmt <- read.delim(SIGNATURE, header = FALSE)
        gmt <- t(as.matrix(gmt))
        colnames(gmt) <- gmt[1, ]
        gmt <- gmt[-1:-2, ]
        message(paste("Total signature records:", ncol(gmt)))
        for (num in (1:ncol(gmt))) {
            GENE <- gmt[, num]
            GENE <- as.character(subset(GENE, GENE != ""))
            message(paste("Target gene_signature:", colnames(gmt)[num]))
            if (!any(GENE %in% rownames(scalef))) {
                GENE <- capitalize(tolower(GENE))
                if (!any(GENE %in% rownames(scalef))) {
                  message(paste("This gene signature is not found in expression list."))
                  next
                }
            }
            GENE <- GENE[GENE %in% rownames(scalef)]
            if (length(GENE) < 2) {
                message("This gene signature is not found in expression list.")
                next
            }
            else {
                texp = colMeans(scalef[GENE, ])
                texp = sort(texp)
                texp_withg = texp[names(texp) %in% rownames(bc_dox_uq) & 
                  !is.na(bc_dox_uq[names(texp), "barcode"])]
                other_table = get_rank_tables_from_rra(texp_withg, 
                  bc_dox_uq, tmpprefix = paste("sample_", runif(1, 
                    1, 10000), sep = ""), rrapath = RRAPATH, 
                  keeptmp = keep_tmp, negctrlgenelist = negctrl_gene)
            }
            if (!is.null(SAVEPATH)) {
                write.table(other_table, file = file.path(paste(colnames(gmt)[num], 
                  "_RRA.txt", sep = "")), sep = "\t", quote = FALSE, 
                  row.names = FALSE)
            }
        }
        setwd(cwd)
    }
    else {
        target_gene_list = strsplit(GENE, ",")[[1]]
        message(paste("Target gene:", paste(target_gene_list, 
            collapse = ";")))
        if (ispathway == TRUE) {
            for (target_gene in target_gene_list) {
                if (!target_gene %in% rownames(scalef)) {
                  message(paste("Error: gene ", target_gene, 
                    " not in expression list."))
                  quit()
                }
            }
            texp = colMeans(scalef[target_gene_list, ])
            texp = sort(texp)
            texp_withg = texp[names(texp) %in% rownames(bc_dox_uq) & 
                !is.na(bc_dox_uq[names(texp), "barcode"])]
            other_table = get_rank_tables_from_rra(texp_withg, 
                bc_dox_uq, tmpprefix = paste("sample_", runif(1, 
                  1, 10000), sep = ""), rrapath = RRAPATH, keeptmp = keep_tmp, 
                negctrlgenelist = negctrl_gene)
            if (!is.null(SAVEPATH)) {
                write.table(other_table, file = file.path(SAVEPATH, 
                  paste(data_label, "_PATHWAY", "_RRA.txt", sep = "")), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
            }
            return(other_table)
        }
        else {
            for (target_gene in target_gene_list) {
                if (!target_gene %in% rownames(scalef)) {
                  message(paste("Warning: gene ", target_gene, 
                    " not in expression list."))
                  next
                }
                else {
                  message(paste("Testing gene ", target_gene, 
                    "..."))
                }
                texp = scalef[target_gene, ]
                texp = sort(texp)
                texp_withg = texp[names(texp) %in% rownames(bc_dox_uq) & 
                  !is.na(bc_dox_uq[names(texp), "barcode"])]
                other_table = get_rank_tables_from_rra(texp_withg, 
                  bc_dox_uq, tmpprefix = paste("sample_", runif(1, 
                    1, 10000), sep = ""), rrapath = RRAPATH, 
                  keeptmp = keep_tmp, negctrlgenelist = negctrl_gene)
                if (!is.null(SAVEPATH)) {
                  write.table(other_table, file = paste(SAVEPATH, 
                    data_label, "_", target_gene, "_RRA.txt", 
                    sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
                }
                return(other_table)
            }
        }
    }
}


get_rank_tables_from_rra = function (rankexp, bc_dox_u, rrapath = NULL, pcutoff = 0.1, tmpprefix = paste("sample_", 
    runif(1, 1, 10000), sep = ""), negctrlgenelist = "NonTargetingControlGuideForHuman", 
    more_rra = "", negsel = T, possel = T, keeptmp = F) 
{
    rankexp = rankexp[names(rankexp) %in% rownames(bc_dox_u) & 
        !is.na(bc_dox_u[names(rankexp), "barcode"])]
    if (length(rankexp) < 3) {
        print("Error: cannot find enough cells.")
        return(NULL)
    }
    rankexp = sort(rankexp)
    texp_guide_ass = bc_dox_u[names(rankexp), "barcode"]
    texp_gene_ass = bc_dox_u[names(rankexp), "gene"]
    texp_guide_ass1 = paste(texp_guide_ass, 1:length(texp_guide_ass), 
        sep = "_r")
    rra_oframe = data.frame(guide = texp_guide_ass1, gene = texp_gene_ass, 
        list = rep("list", length(texp_guide_ass)), value = rankexp, 
        prob = rep(1, length(texp_guide_ass)), chosen = rep(1, 
            length(texp_guide_ass)))
    low_file = paste(tmpprefix, "_rra_low.txt", sep = "")
    write.table(rra_oframe, file = low_file, row.names = F, quote = F, 
        sep = "\t")
    rra_oframe_h = rra_oframe
    rra_oframe_h[, "value"] = -1 * rra_oframe_h[, "value"]
    rra_oframe_h = rra_oframe_h[order(rra_oframe_h[, "value"]), 
        ]
    high_file = paste(tmpprefix, "_rra_high.txt", sep = "")
    write.table(rra_oframe_h, file = high_file, row.names = F, 
        quote = F, sep = "\t")
    ngguidefile = paste(tmpprefix, "_negctrl.txt", sep = "")
    if (!is.null(negctrlgenelist)) {
        ngguidelist = texp_guide_ass1[texp_gene_ass %in% negctrlgenelist]
        write.table(ngguidelist, file = ngguidefile, sep = "\t", 
            row.names = F, col.names = F, quote = F)
        ngguidecommand = paste("--control", ngguidefile)
    }
    else {
        ngguidecommand = ""
    }
    if (is.null(rrapath)) {
        rracommand = "RRA"
    }
    else {
        rracommand = rrapath
    }
    rra_low_out = paste(tmpprefix, "_rra_low.out", sep = "")
    rra_c = paste(rracommand, "-i", low_file, "-o", rra_low_out, 
        ngguidecommand, "-p", pcutoff, "--max-sgrnapergene-permutation 100 ", 
        more_rra)
    if (negsel) {
        print(rra_c)
        system(rra_c, ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    rra_high_out = paste(tmpprefix, "_rra_high.out", sep = "")
    rra_c = paste(rracommand, "-i", high_file, "-o", rra_high_out, 
        ngguidecommand, "-p", pcutoff, "--max-sgrnapergene-permutation 100 ", 
        more_rra)
    if (possel) {
        print(rra_c)
        system(rra_c, ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    if (negsel) {
        frame_l = read.table(rra_low_out, header = T, as.is = T, 
            row.names = 1, na.strings = "")
    }
    if (possel) {
        frame_h = read.table(rra_high_out, header = T, as.is = T, 
            row.names = 1, na.strings = "")
    }
    if (negsel & !possel) {
        system(paste("rm", low_file, rra_low_out))
        if (!is.null(negctrlgenelist)) {
            system(paste("rm", ngguidefile))
        }
        return(frame_l)
    }
    if (!negsel & possel) {
        system(paste("rm", high_file, rra_high_out))
        if (!is.null(negctrlgenelist)) {
            system(paste("rm", ngguidefile))
        }
        return(frame_h)
    }
    report_f = merge(frame_l, frame_h, by = 0, suffixes = c(".low", 
        ".high"))
    if (!keeptmp) {
        system(paste("rm", low_file, high_file, rra_low_out, 
            rra_high_out))
        if (!is.null(negctrlgenelist)) {
            system(paste("rm", ngguidefile))
        }
    }
    return(report_f)
}





dotplot_beta_PIP <- function(fit,
                             target_names = NULL,
                             reorder_targets = target_names,
                             reorder_factors = NULL,
                             exclude_offset = TRUE){
  if (!inherits(fit, "gsfa_fit")){
    stop("Input argument \"fit\" should be an object of class ",
         "\"gsfa_fit\", such as an output of fit_gsfa_multivar().")
  }
  beta_pip <- t(fit$posterior_means$Gamma_pm) # factor by target association PIP matrix
  beta_pm <- t(fit$posterior_means$beta_pm) # factor by target association effect matrix
  if (exclude_offset){
    beta_pip <- beta_pip[, -ncol(beta_pip)]
    beta_pm <- beta_pm[, -ncol(beta_pm)]
  }
  if (is.null(target_names)){
    target_names <- colnames(beta_pm)
  } else {
    colnames(beta_pip) <- target_names
    colnames(beta_pm) <- target_names
  }
  if (is.null(reorder_targets)){
    reorder_targets <- target_names
  }

  beta_pip <- beta_pip[, reorder_targets]
  beta_pip_df <- as.data.frame(beta_pip)
  beta_pip_df$Factor <- paste0("Factor ", 1:nrow(beta_pip_df))
  if (!is.null(reorder_factors)){
    beta_pip_df <- beta_pip_df[reorder_factors, ]
  }
  beta_pip_plot_df <- melt(beta_pip_df, value.name = "PIP")

  beta_pm <- beta_pm[, reorder_targets]
  beta_pm_df <- as.data.frame(beta_pm)
  beta_pm_df$Factor <- paste0("Factor ", 1:nrow(beta_pm_df))
  if (!is.null(reorder_factors)){
    beta_pm_df <- beta_pm_df[reorder_factors, ]
  }
  beta_pm_plot_df <- melt(beta_pm_df, id.var = "Factor",
                                    variable.name = "Perturbation",
                                    value.name = "Effect size")
  beta_pm_plot_df <- beta_pm_plot_df %>%
    mutate(PIP = beta_pip_plot_df$PIP,
           Perturbation = factor(Perturbation, levels = reorder_targets))
  beta_pm_plot_df <- beta_pm_plot_df %>%
    mutate(Factor = factor(Factor, levels = beta_pip_df$Factor[1:nrow(beta_pip_df)]))

  plot_out <- ggplot(beta_pm_plot_df) +
    geom_point(aes(x = Factor, y = Perturbation,
                   size = PIP, color = `Effect size`)) +
    scale_color_gradient2(low = "purple3", mid = "grey90", high = "darkorange1") +
    # color scale can be changed externally
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13, hjust = 1),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12))
  return(plot_out)
}

dotplot_total_effect <- function(fit, gene_indices, gene_names = gene_indices,
                                 target_names = NULL, reorder_targets = NULL,
                                 plot_max_score = NULL){
  # Both inputs should be gene by guide/marker matrices,
  # Dots will be colored by effect size and sized according to LFSR bins.
  require(dplyr)
  require(ggplot2)
  lfsr_binning <- function(lfsr){
    if (lfsr <= 0.05){
      return("0 - 0.05")
    } else if (lfsr <= 0.25){
      return("0.05 - 0.25")
    } else {
      return("> 0.25")
    }
  }
  lfsr_matrix <- fit$lfsr[, -ncol(fit$lfsr)] # remove offset
  if (is.null(target_names)){
    target_names <- colnames(lfsr_matrix)
  } else {
    colnames(lfsr_matrix) <- target_names
  }
  if (is.null(reorder_targets)){
    reorder_targets <- target_names
  }
  selected_lfsr_mat <- lfsr_matrix[gene_indices, reorder_targets]
  rownames(selected_lfsr_mat) <- gene_names

  effect_matrix <- fit$posterior_means$W_pm %*%
    t(fit$posterior_means$beta_pm[-nrow(fit$posterior_means$beta_pm), ])
  colnames(effect_matrix) <- target_names
  selected_effect_mat <- effect_matrix[gene_indices, reorder_targets]
  rownames(selected_effect_mat) <- gene_names

  pip_df <- as.data.frame(selected_lfsr_mat)
  pip_df$gene <- rownames(selected_lfsr_mat)
  pip_plot_df <- melt(pip_df, variable.name = "Perturbation",
                      value.name = "LFSR")

  effect_df <- as.data.frame(selected_effect_mat)
  effect_df$gene <- rownames(selected_effect_mat)
  effect_plot_df <- melt(effect_df, variable.name = "Perturbation",
                         value.name = "Effect_size")

  combined_plot_df <- pip_plot_df %>%
    mutate(Effect_size = effect_plot_df$Effect_size,
           gene = factor(gene, levels = gene_names),
           Perturbation = factor(Perturbation, levels = reorder_targets))
  ## Stratify LSFR values into discrete bins
  combined_plot_df <- combined_plot_df %>%
    rowwise() %>%
    mutate(LFSR_bin = lfsr_binning(LFSR)) %>%
    mutate(LFSR_bin = factor(LFSR_bin, levels = c("> 0.25", "0.05 - 0.25", "0 - 0.05")))

  if (!is.null(plot_max_score)){
    ## Capping effect size values on both ends for more friendly visualization
    ## useful when the data contain only a few extreme values
    plot_min_score <- plot_max_score * (-1)
    plot_by_score <- plot_max_score/2
    combined_plot_df$Effect_size[combined_plot_df$Effect_size > plot_max_score] <- plot_max_score
    combined_plot_df$Effect_size[combined_plot_df$Effect_size < plot_min_score] <- plot_min_score
  }

  plot_out <- ggplot(combined_plot_df) +
    geom_point(aes(x = Perturbation, y = gene,
                   size = LFSR_bin, color = Effect_size)) +
    theme_void() +
    theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 13),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12)) +
    labs(color = "Summarized effect", size = "LFSR")
  if (!is.null(plot_max_score)){
    plot_out <- plot_out +
      scale_color_gradientn(limits = c(plot_min_score, plot_max_score),
                            colours = c("blue3", "blue", "grey90", "red", "red3"),
                            breaks = seq(plot_min_score, plot_max_score, plot_by_score))
  } else {
    plot_out <- plot_out +
      scale_color_gradient2(low = "blue3", mid = "grey90", high = "red3")
  }
  return(plot_out)
}



sgRNA_quality_plot <- function(mtx, sg_lib, title.size = 16, legend.text.size = 8, legend.title.size = 8, x.text.size = 8, x.title.size = 8, y.text.size = 8, y.title.size = 8, label.size = 6, bar_width = NULL, plot.save = TRUE, width = 8, height = 8) {
    

    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T, sep='\t')
    }

    #remove cells in sgRNA library that are not included in matrix
    
    sg_lib_filtered <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))
    mtx = mtx[, colnames(mtx) %in% intersect(sg_lib$cell, colnames(mtx))]
    mtx <- Add_meta_data(sg_lib = sg_lib_filtered,  mtx = mtx)



    #count sgRNA in each cell;count cells of each sgRNA
    sg_count <- plyr::count(subset(sg_lib_filtered, cell %in% colnames(mtx))$barcode)
    colnames(sg_count) <- c("sgRNA", "freq")
    sg_count <- sg_count[order(-sg_count$freq), ]
    sg_count$order <- seq(1, nrow(sg_count))

    #prepare for plot each gene
    
    sg_gene <- unique(sg_lib[ ,c("barcode", "gene")])
    rownames(sg_gene) <- sg_gene$barcode
    sg_count$gene <- sg_gene[sg_count$sgRNA, "gene"]
    sg_count['Perturbation'] = sg_count['sgRNA']

    g1 <- ggplot(data = sg_count, mapping = aes(x = order, y = freq)) +
    geom_line(size = 1) + theme_test() +
    geom_point(data = head(sg_count, 30),
               mapping = aes(x = order, y = freq, color = Perturbation), size = 5) +
    labs(x = "Rank", y = "Cell Count", title = "Cell Count of Perturbation") +
    theme(plot.title = element_text(hjust = 0.5, size = title.size), 
          legend.text = element_text(size = legend.text.size),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          legend.title = element_text(size = legend.title.size),
          axis.title.x = element_text(size = x.title.size), 
          axis.title.y = element_text(size = y.title.size),
          axis.text.y = element_text(hjust = 0.5, size = y.text.size)) +
    guides(colour = guide_legend(title.hjust = 0.5))
  
    sg_num_count <- plyr::count(mtx[["sgRNA_num"]])
    if (max(mtx$sgRNA_num) > 10) {
        over_10 <- subset(sg_num_count, sgRNA_num > 10)
        sg_num_count <- subset(sg_num_count, sgRNA_num <= 10)
        over_10 <- sum(over_10$freq)
        sg_num_count <- rbind(sg_num_count, c(">10", over_10))
    }
    colnames(sg_num_count) <- c("sgRNA_num", "freq")
    sg_num_count$freq <- as.numeric(sg_num_count$freq)
    sg_num_count$sgRNA_num <- factor(sg_num_count$sgRNA_num, levels = sg_num_count$sgRNA_num, ordered = T)

    g2 <- ggplot(sg_num_count,mapping = aes(x = sgRNA_num, y = freq)) +
    geom_bar(stat = "identity",color = "black", fill = "#8491B4FF", width = bar_width) + theme_test() + 
    labs(x = "Perturbation Count", y = "Cell Count", title = "Perturbation in Each Cell") +
    theme(plot.title = element_text(hjust = 0.5, size = title.size),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = x.text.size),
          axis.title.x = element_text(size = x.title.size), 
          axis.title.y = element_text(size = y.title.size),
          axis.text.y = element_text(hjust = 0.5, size = y.text.size)) +
    geom_text(aes(label = freq), stat = "identity", vjust = -0.5, size = label.size)
    
    #save plot
    
    if (plot.save == TRUE) {
        
        dir <- file.path('figures')
        if (!(dir.exists(dir))) {
            dir.create(dir)
        }
        
        combined_plots <- arrangeGrob(g1, g2, ncol = 2)
        ggsave("figures/sgRNA_quality.pdf", combined_plots, width = width * 2, height = height, dpi = 300)
    }
}



improved_scmageck_lr_edit <- function(BARCODE, RDS, NEGCTRL = "NTC", SELECT_GENE = NULL, LABEL = NULL, PERMUTATION = NULL, SAVEPATH = ".", LAMBDA = 0.01, NTC_baseline = TRUE) {
    if (!is.null(LABEL)) {
        data_label = LABEL
    } else {
        data_label = "sample1"
    }

    if (!is.null(PERMUTATION)) {
        n_permutation = as.integer(PERMUTATION)
    } else {
        n_permutation = 10000
    }

    # read cell assignment and libray file ####
    
    if (is.character(BARCODE)) {
        bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
    } else {
        bc_dox = BARCODE
    }
    
    

    if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
        stop("cell, barcode, or gene column names not found in barcode file.")
    }

    message(paste("Total barcode records:", nrow(bc_dox)))

    # load neg control guides ####
    
    ngctrlgenelist = strsplit(NEGCTRL, ",")[[1]] #NTC need to be a string split by ","
    message(paste("Neg Ctrl guide:", paste(ngctrlgenelist, collapse = ";")))

    # read Seurat RDS file ####
    
    if (is.character(RDS)) {
        message(paste("Reading RDS file:", RDS))
        targetobj = readRDS(RDS)
    } else {
        targetobj = RDS
    }
    
    if (isS4(targetobj)) {
        targetobj <- SeuratObject::GetAssayData(object = targetobj, slot = "scale.data")
    }
    
    # check if names are consistent
    
    nmatch = sum(bc_dox[, "cell"] %in% colnames(x = targetobj))
    if (nmatch == 0) {
        stop("No cell names match in expression matrix and barcode file.")
    }
    bc_dox <- subset(bc_dox, cell %in% colnames(x = targetobj))

    # convert to ind_matrix ####
    
    ind_matrix <- frame2indmatrix(bc_dox, targetobj) #return TRUE and FALSE matrix
    message(paste("Index matrix dimension:", nrow(ind_matrix), ",", ncol(ind_matrix)))

    # try to perform matrix regresson on single genes ####
    
    mat_for_single_reg = single_gene_matrix_regression(targetobj, selected_genes_list = SELECT_GENE,
                                                       ngctrlgene = ngctrlgenelist, indmatrix = ind_matrix, 
                                                       NTC_baseline = NTC_baseline)
    Xmat = mat_for_single_reg[[1]]

    # Xmat[,which(colnames(Xmat)%in%ngctrlgenelist)[1]]=1 # already integrated into function
    
    Ymat = mat_for_single_reg[[2]]
    
    Ymat = Ymat[, rownames(Xmat)]
    rm(list = c("mat_for_single_reg", "targetobj", "ind_matrix"))

    # remove values in Y mat
    
    Amat_pm_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, Ymat, lambda = LAMBDA, 
                                                              npermutation = n_permutation)
    Amat = Amat_pm_lst[[1]]
    Amat_pval = Amat_pm_lst[[2]]
   
    if (!is.null(SAVEPATH)) {
        write.table(data.frame(Perturbedgene = rownames(Amat), Amat, check.names=FALSE), 
                    file = file.path(SAVEPATH, paste(data_label, "_score.txt", sep = "")), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(data.frame(Perturbedgene = rownames(Amat), Amat_pval, check.names=FALSE), 
                    file = file.path(SAVEPATH, paste(data_label, "_score_pval.txt", sep = "")), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(list(data.frame(Perturbedgene = rownames(Amat), Amat, check.names=FALSE), 
                data.frame(Perturbedgene = rownames(Amat), Amat_pval, check.names=FALSE)))
}

frame2indmatrix <- function(bc_d, targetobj) {
    #rnm = unique(bc_d$cell)
    rnm = colnames(targetobj)
    cnm = unique(bc_d$gene)
    scalef = targetobj
    message(paste(length(unique(bc_d$cell)), "..."))
    message(paste(ncol(scalef), "..."))
    #rnm = rnm[!is.na(rnm)] #remove NA
    #rnm = rnm[rnm %in% colnames(scalef)]
    if (length(rnm) == 0) {
        stop("Cell names do not match in expression matrix and barcode.")
    }
    cnm = cnm[!is.na(cnm)]#remove NA
    ind_matrix = matrix(rep(0, length(rnm) * length(cnm)), nrow = length(rnm))
    rownames(ind_matrix) = rnm
    colnames(ind_matrix) = cnm
    row <- bc_d[, 'cell']
    col <- bc_d[, 'gene']
    test <- (row %in% rnm) & (col %in% cnm)
    idx <- as.matrix(data.frame(row[test], col[test]))
    #idx <- cbind(row[test], col[test])
    ind_matrix[idx] <- 1
    return(ind_matrix)
}

single_gene_matrix_regression <- function(targetobj, ngctrlgene = c("NonTargetingControlGuideForHuman"),
                                          indmatrix = NULL, selected_genes_list = NULL, NTC_baseline = TRUE) {
    
    # return X matrix and Y matrix for regression note that all the ngctrlgene are merged into one
    # column, 'NegCtrl' if indmatrix is provided, the Xmat will be constructed from indmatrix
    
    outlier_threshold = 0.95

    #because select cells and genes before input,in this step, there is no need to select
    
    YmatT = targetobj
    select_genes <- rownames(YmatT)
    if (!is.null(selected_genes_list)) {
        select_genes = select_genes[select_genes %in% selected_genes_list]
        if (length(select_genes) == 0) {
            stop("No genes left for regression. Check your selected gene list.")
        }
    }
    message(paste("Selected genes:", length(select_genes)))

    if (NTC_baseline == TRUE) {
        if (length(ngctrlgene) == 1) {
            colnames(indmatrix)[colnames(indmatrix) == ngctrlgene]<- "NegCtrl"
            indmatrix[, "NegCtrl"] <- 1
        } else if (length(ngctrlgene) > 1) {
            indmatrix <- indmatrix[, -which(colnames(indmatrix) %in% ngctrlgene[-1])] #only remain one columns of NTC
            colnames(indmatrix)[colnames(indmatrix) == ngctrlgene[1]] <- "NegCtrl"
            indmatrix[, "NegCtrl"] <- 1
        }
    }

    YmatT = YmatT[select_genes, ]
    
    # remove outliers
        
    Ymat_outlier = apply(YmatT, 1, function(X) {
        return(quantile(X, probs = outlier_threshold))
    })
    YmatT <- YmatT - Ymat_outlier
    YmatT[YmatT > 0] <- 0
    YmatT <- YmatT + Ymat_outlier

    return(list(indmatrix, YmatT))
}


getsolvedmatrix_with_permutation_cell_label <- function(Xm, Ym, lambda = 0.01, npermutation = 1000) {
    Amat_ret = getsolvedmatrix(Xm, Ym, lambda = lambda)
    Amat_ret_higher = Amat_ret * 0
    
    # permute N times randomly shuffle cell labels
    
    for(npm in 1:npermutation){
        if (npm%%100 == 0) {
            message(paste("Permutation:", npm, "/", npermutation, "..."))
    }
    cells_shu = sample(colnames(Ym), ncol(Ym))
    Xm_s = Xm[cells_shu, ]
    Amat_random = getsolvedmatrix(Xm_s, Ym, lambda = lambda)

    Amat_ret_higher = Amat_ret_higher + (abs(Amat_random) > abs(Amat_ret)) * 1
    }
    Amat_ret_higher = Amat_ret_higher/npermutation
    return(list(Amat_ret, Amat_ret_higher))
}



getsolvedmatrix <- function(Xm, Ym, lambda = 0.01) {
    # Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
    TMmat_g = crossprod(Xm, Xm) + lambda * diag(ncol(Xm))

    Amat_g = Ym %*% Xm %*% solve(TMmat_g)
    return(Amat_g)
}