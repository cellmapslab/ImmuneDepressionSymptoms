library(dplyr)
library(purrr)

source("03_cytokine_clustering/01_code/plotting_functions.R")

generate_report <- function(
    type,
    cluster_assignments,
    general_data,
    add_healthy_controls = FALSE,
    ...  
) {
  dots <- list(...)
  
  clean_cluster_assignments <-
    get_cluster_assignments(
      cluster_assignments = cluster_assignments,
      general_data = general_data,
      add_healthy_controls = add_healthy_controls
    )
  
  do.call(paste0("report_", type), c(list(data = general_data,
                                          cluster_assignments_df = clean_cluster_assignments$cluster_assignments_df),
                                     dots))
}

report_variable_importance <- function(
    data,
    cluster_assignments_df,
    omic_type = c("genesets", "cytokines_RNA", "RNA"),
    cytokine_version = c("standard", "additional", "corrected_1", "corrected_2",
                         "additional_unimputed"),
    rna_version = c("standard", "corrected_1"),
    ...
) {
  
  if (omic_type == "genesets") {
    report_variable_importance_genesets(
      data = data,
      cluster_assignments_df = cluster_assignments_df,
      ...
    )
  } else if (omic_type == "cytokines_RNA") {
    report_variable_importance_cytokines_rna(
      data = data,
      cluster_assignments_df = cluster_assignments_df,
      cytokine_version = cytokine_version,
      rna_version = rna_version,
      ...
    )
  } else if (omic_type == "RNA") {
    report_variable_importance_rna(
      data = data,
      cluster_assignments_df = cluster_assignments_df,
      cytokine_version = cytokine_version,
      rna_version = rna_version,
      ...
    )
  } else {
    stop("This omic type is not implemented yet")
  }
}

report_variable_importance_cytokines_rna <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    cytokine_version = c("standard", "additional", "corrected_1", "corrected_2",
                         "additional_unimputed"),
    rna_version = c("standard", "corrected_1"),
    show_gene_names = FALSE,
    include_age_bmi = FALSE,
    include_celltypes = FALSE,
    ...
) {
  variable_importance <- calculate_variable_importance_cytokines_rna(
    data = data,
    cluster_assignments_df = cluster_assignments_df,
    cluster_colours = cluster_colours,
    cytokine_version = cytokine_version,
    rna_version = rna_version,
    show_gene_names = show_gene_names,
    include_age_bmi = include_age_bmi,
    include_celltypes = include_celltypes
  )
  
  variable_importance %>% 
    arrange(desc(F_value))
}

report_variable_importance_genesets <- function(
    data,
    cluster_assignments_df,
    geneset_results,
    include_p_value = FALSE,
    ...
) {
  
  variable_importance_data <- geneset_results %>% 
    right_join(cluster_assignments_df, by = "ID") %>%
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    dplyr::select(-ID)
  
  geneset_names <- colnames(variable_importance_data)
  geneset_names <- geneset_names[geneset_names != "cluster"]
  
  variable_importance <- set_names(geneset_names) %>% 
    map_dfr(function(x) {
      lm_res <- lm(variable_importance_data[[x]] ~ variable_importance_data[["cluster"]])
      anova_res <- anova(lm_res)
      res <- data.frame(F_value = anova_res[["F value"]][1])
      if (include_p_value) {
        res$p_value <- anova_res[["Pr(>F)"]][1]
      }
      res
    }, .id = "geneset")
  
  variable_importance %>% 
    arrange(desc(F_value))
}

report_variable_importance_rna <- function(
    data,
    cluster_assignments_df,
    rna_version = c("standard", "corrected_1"),
    ...
) {
  rename_correction <- function(name) {
    stringr::str_extract(name, ".+?(?=_corrected)")
  }
  if (rna_version == "standard") {
    variable_importance_data <- data$rna %>% 
      pivot_wider(
        id_cols = ID,
        names_from = "gene",
        values_from = "expression"
      ) %>% 
      right_join(cluster_assignments_df, by = "ID") %>%
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      dplyr::select(-ID)
  } else if (rna_version == "corrected_1"){
    variable_importance_data <- data$rna_corrected_1 %>% 
      pivot_wider(
        id_cols = ID,
        names_from = "gene",
        values_from = "expression"
      ) %>% 
      rename_with(rename_correction, starts_with("ENSG")) %>% 
      right_join(cluster_assignments_df, by = "ID") %>%
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      dplyr::select(-ID)
  } else {
    stop("Not implemented yet.")
  }
  
  gene_names <- colnames(variable_importance_data)
  gene_names <- gene_names[gene_names != "cluster"]
  
  variable_importance <- set_names(gene_names) %>% 
    map_dfr(function(x) {
      lm_res <- lm(variable_importance_data[[x]] ~ variable_importance_data[["cluster"]])
      anova_res <- anova(lm_res)
      data.frame(F_value = anova_res[["F value"]][1])
    }, .id = "gene")
  
  variable_importance %>% 
    arrange(desc(F_value))
}

report_cytokine_RNA_correlation <- function(
    data,
    cluster_assignments_df,
    cytokine_version = c("standard", "additional", "corrected_1", "additional_unimputed"),
    rna_type = c("standard", "corrected_1"),
    sample_type = c("cases", "controls", "lifetime cases",  "all"),
    pretty_names = FALSE,
    ...
) {
  corr_data <- calculate_cytokine_RNA_correlation(
    data = data,
    cluster_assignments_df = cluster_assignments_df,
    cytokine_version = cytokine_version,
    rna_type = rna_type,
    sample_type = sample_type,
    pretty_names = pretty_names,
    ...
  )
  
  res <- as.data.frame(corr_data)
  res$variable_1 <- rownames(res)
  res <- res %>% 
    pivot_longer(
      cols = -variable_1,
      names_to = "variable_2",
      values_to = "correlation"
    )
  
  res
}

report_RNA_HIRDS_distribution <- function(
    data,
    cluster_assignments_df,
    rna_type = c("standard", "corrected_1"),
    show_gene_id = TRUE,
    ...
) {
  if (rna_type == "standard") {
    plot_data <- data$rna
  } else if (rna_type == "corrected_1") {
    plot_data <- data$rna_corrected_1 %>% 
      dplyr::mutate(gene = str_extract(gene, ".+?(?=_corrected)"))
  } else {
    stop("This rna_type is not supported.")
  }
  
  if ("gene_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    rna_index <- plot_data$gene %in% match.call()[["gene_list"]]
    plot_data <- plot_data[rna_index, ]
  } else {
    stop("Please provide a `gene_list` for which genes the expression should be shown.")
  }
  
  # check for which genes the median of the HIRDS clusters are either lower
  # or higher than MIDS/LIRDS
  res <- plot_data %>% 
    inner_join(cluster_assignments_df, by = "ID") %>% 
    filter(cluster != "control") %>% 
    group_by(gene, cluster) %>% 
    summarise(median_expression = median(expression, na.rm = TRUE)) %>% 
    pivot_wider(
      id_cols = gene,
      names_from = cluster,
      values_from = median_expression
    ) %>% 
    filter((`cluster 1` < `cluster 2` & `cluster 1` < `cluster 4` & `cluster 3` < `cluster 2` & `cluster 3` < `cluster 4`) |
             (`cluster 1` > `cluster 2` & `cluster 1` > `cluster 4` & `cluster 3` > `cluster 2` & `cluster 3` > `cluster 4`)) %>% 
    # add the gene name info 
    left_join(data$rna_info %>% select(ensembl_gene_id, external_gene_name),
              by = c("gene" = "ensembl_gene_id"))
  res
}
