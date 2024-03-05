library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(purrr)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(forcats)
library(ggcorrplot)

read_in_general_data <- function(
  proband_type = c("current_cases_new_cidi"),
  # additional_correct also includes the updated BMI data
  cytokine_version = c("additional_correct",
                       # corrected for age, sex, BMI
                       "additional_corrected_1"),
  add_healthy_controls = TRUE,
  healthy_control_version = "strict",
  add_lifetime_cases = FALSE
) {
  all_data <- list()
  if (cytokine_version == "additional_correct") {
    all_data$imputed_data <- readRDS("02_imputed_data_add_cytokines_2022_10_28_3.rds")
  } else if (cytokine_version == "additional_corrected_1"){
    imputed_data <- readRDS("cytokine_data_imputed_corrected_age_sex_bmi.rds")
    rename_correction <- function(name) {
      stringr::str_extract(name, ".+?(?=_corrected)")
    }
    all_data$imputed_data <- imputed_data %>% 
      rename_with(rename_correction, bFGF_corrected_1:sIL_6R_corrected_1)
    rm(imputed_data)
  } else {
    stop("This cytokine_version is not available.")
  }
  
  # cytokine data uncorrected
  all_data$cytokine_uncorrected <- readRDS("01_cytokine_data_for_imputation_3.rds")
  # cytokine data corrected for age/sex/BMI
  all_data$cytokine_corrected_1 <- readRDS("cytokines_corrected_for_age_sex_bmi_norm_1.rds")
  
  all_data$cidi_data <- readRDS("cidi_has_diagnosis.Rds")
  all_data$cidi_single <- readRDS("cidi_single_tidy_processed_jonas.rds")
  # case control status according to CIDI
  all_data$case_control <- readRDS("ids_cases_and_controls_for_clustering_new_cidi_case_control_status_2022_11_23.Rds")
  all_data$cidi_depression <- readRDS("cidi_depression_diagnosis.Rds")
  all_data$ids_analysis <- readRDS("02_data/ids_for_analysis.rds")
  
  if (proband_type == "current_cases_new_cidi") {
    all_data$ids_joined <- readRDS("ids_cases_for_clustering_new_cidi_2022_11_23.Rds")
    all_data$proband_status <- data.frame(
      ID = all_data$ids_joined$ID,
      status = "current case"
    )
  }
  
  if (add_healthy_controls) {
    if (healthy_control_version == "strict") {
      all_data$control_ids <- readRDS("cidi_no_diagnosis_strict_2022_11_25.Rds")
    }
    
    ids_for_filtering <- bind_rows(
      all_data$ids_joined,
      all_data$control_ids
    )
  } else {
    ids_for_filtering <- all_data$ids_joined
  }
  
  if (add_lifetime_cases) {
    all_data$lifetime_ids <- readRDS("cidi_only_lifetime_diagnosis_2022_12_30.rds")
    
    ids_for_filtering <- bind_rows(
      ids_for_filtering,
      all_data$lifetime_ids
    )
  }
  
  all_data$imputed_data <- all_data$imputed_data %>% 
    filter(ID %in% ids_for_filtering$ID)
  
  # add updated age info
  all_data$imputed_data <- all_data$imputed_data %>% 
    dplyr::select(-age)
  all_data$age <- readRDS("age_data.Rds")

  all_data$celltype_data_rna_updated_14 <- readRDS("deconvolution_featurecounts_raw_dtangle_lm22_14_celltypes.rds")
  
  
  all_data$medication_data <- readRDS("medication_data_both.Rds")
  all_data$medication_data <- all_data$medication_data %>% 
    right_join(ids_for_filtering, by = "ID")
  
  all_data$medication_become_categories <- readRDS("medication_data_become_categories_tanja_2022_11_17.Rds")
  all_data$medication_become_categories <- all_data$medication_become_categories %>% 
    right_join(ids_for_filtering, by = "ID")
  
  all_data$somatic_data <- readRDS("somatic_disease_overview.Rds")
  all_data$somatic_data_single <- readRDS("somatic_disease_single_items.Rds")
  all_data$somatic_data_single <- all_data$somatic_data_single %>% 
    right_join(ids_for_filtering, by = "ID")
  
  # BDI data
  all_data$bdi_data <- readRDS("bdi_single_items.Rds")
  all_data$bdi_data <- all_data$bdi_data %>% 
    right_join(ids_for_filtering, by = "ID")
  
  # additional BDI sumscores
  all_data$bdi_sumscores <- readRDS("bdi_sumscores.Rds")
  all_data$bdi_sumscores <- all_data$bdi_sumscores %>% 
    right_join(ids_for_filtering, by = "ID")
  
  # RNA data
  rna <- readRDS("rna_filtered_batch_corrected_normalised_featurecounts.rds")
  rna <- as.data.frame(t(rna))
  rna$combined_id <- rownames(rna)
  id_matching <- readRDS("matched_ids.Rds")
  all_data$rna <- rna %>% 
    left_join(id_matching %>% dplyr::select(ID, combined_id, library_rna_amount), by = "combined_id") %>% 
    filter(library_rna_amount == 200) %>% 
    dplyr::select(-c(combined_id, library_rna_amount)) %>% 
    right_join(ids_for_filtering, by = "ID") %>% 
    pivot_longer(
      cols = -ID,
      names_to = "gene",
      values_to = "expression"
    )
  
  all_data$rna_info <- readRDS("all_rna_info.rds")
  
  # RNA corrected for age/sex/BMI
  rna_corrected_1 <- readRDS("rna_corrected_for_age_sex_bmi_norm_1.rds")
  all_data$rna_corrected_1 <- rna_corrected_1 %>% 
    right_join(ids_for_filtering, by = "ID") %>% 
    pivot_longer(
      cols = -ID,
      names_to = "gene",
      values_to = "expression"
    )
  
  # MADRS
  madrs <- readRDS("madrs_updated.Rds")
  all_data$madrs <- madrs %>% 
    right_join(ids_for_filtering, by = "ID")
  
  # cytokine panel information
  data_msd_raw <- readRDS("01_dataMSD_all_cols_Jonas.rds")
  all_data$cytokine_panel_info <- data_msd_raw %>% 
    dplyr::select(Assay, panel_name) %>% 
    distinct(Assay, .keep_all = TRUE) %>% 
    mutate(Assay = str_replace_all(Assay, "[ -/]", "_")) %>% 
    # add info for additional markers
    bind_rows(
      data.frame(Assay = c("hsCRP", "IL_6HS"),
                 panel_name = c("Vascular Injury Panel", "Proinflammatory Panel"))
    ) %>% 
    dplyr::rename(cytokine = Assay,
           `panel name` = panel_name)
  
  # ECG data
  all_data$ecg <- readRDS("ecg_data_2023_08_31.rds")
  all_data$ecg_wo_excl <- readRDS("ecg_data_wo_excl_2023_08_31.rds")
  
  # Imaging data
  thickness_data_revised <- readRDS("thickness_data_clean_2.rds")
  all_data$imaging_thickness_revised <- thickness_data_revised %>% 
    right_join(ids_for_filtering, by = "ID")
  
  volumina_data <- readRDS("volumina_data_clean.rds")
  all_data$imaging_volumina <- volumina_data %>% 
    right_join(ids_for_filtering, by = "ID")
  
  residuals_ventricles <- readRDS("residuals_ventricles.rds")
  all_data$imaging_residuals_ventricles <- residuals_ventricles %>% 
    right_join(ids_for_filtering, by = "ID")
  
  # BMI data
  # the same data as contained in "cytokine data uncorrected"
  bmi_data <- readRDS("01_cytokine_data_for_imputation_3.rds")
  all_data$bmi <- bmi_data %>% 
    select(ID, BMI) %>% 
    right_join(ids_for_filtering, by = "ID")
  
  all_data
}

set_defaults <- function(
    compare_means_test = "wilcox.test",
    add_lifetime_cases = FALSE,
    add_healthy_controls = TRUE
) {
  list(
    compare_means_test = compare_means_test,
    add_lifetime_cases = add_lifetime_cases,
    add_healthy_controls = add_healthy_controls
  )
}

add_dashes <- function(
    p,
    plot_data
) {
  if ("controls" %in% plot_data$cluster && "lifetime cases" %in% plot_data$cluster) {
    intercept_position <- length(unique(plot_data$cluster)) - 1.5
    p <- p +
      geom_vline(aes(xintercept = intercept_position),
                 linetype = "longdash",
                 colour = "#555555")
  }
  
  if (xor("controls" %in% plot_data$cluster, "lifetime cases" %in% plot_data$cluster)) {
    intercept_position <- length(unique(plot_data$cluster)) - 0.5
    p <- p +
      geom_vline(aes(xintercept = intercept_position),
                 linetype = "longdash",
                 colour = "#555555")
  }
  
  p
}

pretty_cytokine_names <- function(plot_data) {
  plot_data %>% 
    dplyr::mutate(cytokine = str_replace(cytokine, "IL_8_HS", "IL-8 HS"),
                  cytokine = str_replace(cytokine, "IL_6_HS", "IL-6 HS"),
                  cytokine = str_replace(cytokine, "VEGF_A_HS", "VEGF-A HS"),
                  cytokine = str_replace(cytokine, "VEGF_A_LS", "VEGF-A LS"),
                  cytokine = str_replace(cytokine, "IL_12_IL_23p40", "IL-12/IL-23p40"),
                  cytokine = str_replace(cytokine, "_", "-"))
}

correct_cytokine_names <- function(plot_data, var_name = cytokine) {
  var_name <- rlang::enquo(var_name)
  plot_data %>% 
    dplyr::mutate({{var_name}} := str_replace({{var_name}}, "bFGF", "FGF2"),
                  {{var_name}} := str_replace({{var_name}}, "Eotaxin$", "CCL11"),
                  {{var_name}} := str_replace({{var_name}}, "Eotaxin_3", "CCL26"),
                  {{var_name}} := str_replace({{var_name}}, "Flt_1", "VEGFR1"),
                  {{var_name}} := str_replace({{var_name}}, "IP_10", "CXCL10"),
                  {{var_name}} := str_replace({{var_name}}, "MCP_1", "CCL2"),
                  {{var_name}} := str_replace({{var_name}}, "MCP_4", "CCL13"),
                  {{var_name}} := str_replace({{var_name}}, "MDC", "CCL22"),
                  {{var_name}} := str_replace({{var_name}}, "MIP_1alpha", "CCL3"),
                  {{var_name}} := str_replace({{var_name}}, "MIP_1beta", "CCL4"),
                  {{var_name}} := str_replace({{var_name}}, "TARC", "CCL17"),
                  {{var_name}} := str_replace({{var_name}}, "TNF_alpha", "TNF"),
                  {{var_name}} := str_replace({{var_name}}, "TNF_alpha", "TNF"),
                  {{var_name}} := str_replace({{var_name}}, "TNF_beta", "LT-alpha"),
                  {{var_name}} := str_replace({{var_name}}, "_", "-"))
}

get_cluster_assignments <- function(
    cluster_assignments,
    general_data = NULL,
    add_healthy_controls = TRUE,
    add_lifetime_cases = FALSE,
    cluster_names = NULL
) {
  
  if (inherits(cluster_assignments, "data.frame")) {
    if (any(!c("ID", "cluster") %in% colnames(cluster_assignments))) {
      stop("cluster_assignments needs the columns ID and cluster")
    }
    
    id_order <- cluster_assignments %>% 
      pull(ID)
    
    cluster_assignments_df <- cluster_assignments %>% 
      dplyr::select(ID, cluster) %>% 
      arrange(cluster, ID)
  } else {
    # to get the ID order as in the clustering:
    id_order <- general_data$imputed_data %>% 
      filter(.imp == 1) %>% 
      right_join(general_data$ids_joined, by = "ID") %>% 
      dplyr::select(ID, bFGF:Cortisol) %>% 
      pull(ID)
    
    if (length(id_order) != length(cluster_assignments)) {
      stop("The length of the cluster assignments and the IDs is not equal.")
    }
    
    cluster_assignments_df <- data.frame(
      ID = id_order,
      cluster = cluster_assignments
    ) %>% 
      arrange(cluster, ID)
  }
  
  if (is.null(cluster_names)) {
    cluster_assignments_df <- cluster_assignments_df %>% 
      dplyr::mutate(cluster = paste0("cluster ", cluster))
  } else {
    cluster_assignments_df <- cluster_assignments_df %>% 
      left_join(cluster_names, by = "cluster") %>% 
      rename(cluster_old = cluster) %>% 
      rename(cluster = cluster_new) %>% 
      dplyr::select(-cluster_old)
  }
  
  
  if (add_healthy_controls) {
    control_info <- general_data$control_ids %>% 
      dplyr::mutate(cluster = "controls") %>% 
      arrange(ID)
    
    cluster_assignments_df <- bind_rows(
      cluster_assignments_df,
      control_info
    )
  }
  
  if (add_lifetime_cases) {
    lifetime_info <- general_data$lifetime_ids %>% 
      dplyr::mutate(cluster = "lifetime cases") %>% 
      arrange(ID)
    
    cluster_assignments_df <- bind_rows(
      cluster_assignments_df,
      lifetime_info
    )
  }
  
  # if I change the cluster names, then I need to specify a factor with a
  # correct order. This can't be done before adding controls or lifetime cases
  if (!is.null(cluster_names)) {
    if  (add_healthy_controls || add_lifetime_cases) {
      if (add_healthy_controls & !add_lifetime_cases) additional_factors <- "controls"
      if (!add_healthy_controls & add_lifetime_cases) additional_factors <- "lifetime cases"
      if (add_healthy_controls & add_lifetime_cases) additional_factors <- c("controls", "lifetime cases")
      correct_order_cluster_names <- c(cluster_names$cluster_new, additional_factors)
    } else {
      correct_order_cluster_names <- cluster_names$cluster_new
    }
    cluster_assignments_df <- cluster_assignments_df %>% 
      mutate(cluster = factor(cluster, levels = correct_order_cluster_names))
  }
  
  correct_order_factor_levels <- cluster_assignments_df$ID
  cluster_assignments_df <- cluster_assignments_df %>% 
    dplyr::mutate(ID = factor(ID, levels = correct_order_factor_levels))
  
  list(
    id_order = id_order,
    cluster_assignments_df = cluster_assignments_df
  )
}

generate_plot <- function(
    type = c("heatmap", "BMI", "age", "sex"),
    cluster_assignments,
    general_data,
    cluster_colours = define_cluster_colours(),
    previous_cluster_colours = define_previous_cluster_colours(),
    cytokine_version = c("standard", "additional"),
    add_healthy_controls = set_defaults()$add_healthy_controls,
    add_lifetime_cases = set_defaults()$add_lifetime_cases,
    cluster_names = NULL,
    ...
) {
  
  dots <- list(...)
  
  clean_cluster_assignments <-
    get_cluster_assignments(
      cluster_assignments = cluster_assignments,
      general_data = general_data,
      add_healthy_controls = add_healthy_controls,
      add_lifetime_cases = add_lifetime_cases,
      cluster_names = cluster_names
    )
  
  # add white for controls
  number_additional_groups <- sum(add_healthy_controls, add_lifetime_cases)
  number_clusters <- length(unique(clean_cluster_assignments[["cluster_assignments_df"]][["cluster"]])) - 
    number_additional_groups
  cluster_colour_selection <- cluster_colours[seq_len(number_clusters)]
  if (add_healthy_controls) {
    cluster_colour_selection <- c(cluster_colour_selection, "controls" = "#FFFFFF")
  }
  if (add_lifetime_cases) {
    cluster_colour_selection <- c(cluster_colour_selection, "lifetime cases" = "#525252")
  }
  
  do.call(paste0("plot_", type),
          c(list(data = general_data,
                 cluster_assignments_df = clean_cluster_assignments$cluster_assignments_df,
                 cluster_colours = cluster_colour_selection,
                 previous_cluster_colours = previous_cluster_colours,
                 id_order = clean_cluster_assignments$id_order,
                 cytokine_version = cytokine_version),
            dots))
}

plot_heatmap <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    previous_cluster_colours,
    id_order,
    cytokine_version = c("additional_unimputed_corrected",
                         "additional_unimputed"),
    previous_assignments = NULL,
    previous_assignments_order = NULL,
    annotation_legend = TRUE,
    legend = TRUE,
    cluster_cols = FALSE,
    show_cases_controls = FALSE,
    show_psychopharmaca = TRUE,
    show_BMI = FALSE,
    var_name_study_instead_med = FALSE,
    included_annotations = NULL,
    pretty_names = FALSE,
    sumo_cytokine_names = FALSE,
    show_cytokine_panel_names = FALSE,
    sort_by_panel_names = FALSE,
    ...
) {
  
  if (!is.null(previous_assignments)) {
    # have the possibility to use a data.frame for previous_assignments
    # (which then can also have a different order) than using cluster and
    # id_order vectors separately
    if (inherits(previous_assignments, "data.frame")) {
      previous_assignments_clean <- get_cluster_assignments(
        cluster_assignments = previous_assignments,
        add_healthy_controls = FALSE,
        add_lifetime_cases = FALSE
      )
      
      previous_assignments_df <- previous_assignments_clean[["cluster_assignments_df"]]
      previous_assignments_df <- previous_assignments_df %>% 
        dplyr::rename(previous_cluster = cluster)
    } else {
      # have the ability to pass different ID order for previous assignments
      # because the data with standard cytokines was differently ordered than
      # the data with additional cytokines
      if (is.null(previous_assignments_order)) {
        old_id_order <- id_order
      } else {
        old_id_order <- previous_assignments_order
      }
      
      previous_assignments_df <- data.frame(
        ID = old_id_order,
        previous_cluster = previous_assignments
      )
    }
    
    cluster_assignments_df <- cluster_assignments_df %>% 
      # only use ids for clustering, no healthy controls or lifetime cases
      filter(cluster != "controls", cluster != "lifetime cases") %>% 
      left_join(previous_assignments_df, by = "ID")
  } else {
    # exclude healthy controls from cluster assignments
    # (needs to be done additionally to removing them from the cytokine values,
    # because otherwise these IDs are added by the right_join with the
    # cluster assignments later)
    cluster_assignments_df <- cluster_assignments_df %>% 
      # only use ids for clustering, no healthy controls
      filter(cluster != "controls", cluster != "lifetime cases")
  }
  
if (cytokine_version == "additional_unimputed_corrected") {
    cytokine_data_median <- data$cytokine_corrected_1 %>% 
      left_join(data$ids_joined, by = "ID") %>% 
      pivot_longer(
        cols = bFGF_corrected_1:sIL_6R_corrected_1,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      ) %>% 
      dplyr::mutate(cytokine = str_extract(cytokine, ".+?(?=_corrected_1)"))
  } else if (cytokine_version == "additional_unimputed") {
    cytokine_data_median <- data$cytokine_uncorrected %>% 
      left_join(data$ids_joined, by = "ID") %>% 
      dplyr::select(-c(Barcode_Rack:BDI_total)) %>% 
      pivot_longer(
        cols = bFGF:sIL_6R,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      )
  } else {
    stop("This cytokine_version is not available.")
  }
  
  cytokine_data_median <- cytokine_data_median %>% 
    group_by(cytokine, ID) %>% 
    summarise(
      median_normalised_concentration = median(normalised_concentration)
    ) %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    arrange(cluster)
  
  panel_name_annotation <- data$cytokine_panel_info
  
  if (pretty_names) {
    cytokine_data_median <- cytokine_data_median %>% 
      dplyr::mutate(cytokine = str_replace(cytokine, "IL_8_HS", "IL-8 HS"),
                    cytokine = str_replace(cytokine, "VEGF_A_HS", "VEGF-A HS"),
                    cytokine = str_replace(cytokine, "VEGF_A_LS", "VEGF-A LS"),
                    cytokine = str_replace(cytokine, "IL_12_IL_23p40", "IL-12/IL-23p40"),
                    cytokine = str_replace(cytokine, "_", "-"))
    
    panel_name_annotation <- panel_name_annotation %>% 
      dplyr::mutate(cytokine = str_replace(cytokine, "IL_8_HS", "IL-8 HS"),
                    cytokine = str_replace(cytokine, "VEGF_A_HS", "VEGF-A HS"),
                    cytokine = str_replace(cytokine, "VEGF_A_LS", "VEGF-A LS"),
                    cytokine = str_replace(cytokine, "IL_12_IL_23p40", "IL-12/IL-23p40"),
                    cytokine = str_replace(cytokine, "_", "-"))
  }
  
  cytokine_data_median_wide <- cytokine_data_median %>% 
    pivot_wider(
      id_cols = ID,
      names_from = "cytokine",
      values_from = "median_normalised_concentration"
    ) %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    arrange(cluster)
  
  cluster_rows = TRUE
  
  if (sort_by_panel_names) {
    cluster_rows = FALSE
  }
  
  # always calculate the cytokine order and apply it (because the order isn't
  # respected if the cytokines/rows should be clustered anyway, so the code
  # is easier)
  cytokine_order <- data.frame(
    cytokine = colnames(cytokine_data_median_wide)
  ) %>% 
    filter(!cytokine %in% c("ID", "cluster", "previous_cluster")) %>% 
    left_join(panel_name_annotation, by = "cytokine") %>% 
    arrange(`panel name`, cytokine) %>% 
    pull(cytokine)
  
  
  if ("previous_cluster" %in% colnames(cytokine_data_median_wide)) {
    heatmap_data <- t(as.matrix(cytokine_data_median_wide %>%
                                  dplyr::select(-c(ID, cluster, previous_cluster)) %>% 
                                  relocate(all_of(cytokine_order))))
  } else {
    heatmap_data <- t(as.matrix(cytokine_data_median_wide %>%
                                  dplyr::select(-c(ID, cluster)) %>% 
                                  relocate(all_of(cytokine_order))))
  }
  
  colnames(heatmap_data) <- cytokine_data_median_wide$ID
  
  # remove the "controls" and "lifetime" colour from the cluster colours
  if ("controls" %in% names(cluster_colours)) {
    cluster_colours <- cluster_colours[names(cluster_colours) != "controls"]
  }
  if ("lifetime cases" %in% names(cluster_colours)) {
    cluster_colours <- cluster_colours[names(cluster_colours) != "lifetime cases"]
  }
  
  # if there are NAs in the clustering assignment, set the last cluster colour
  # for NA to grey
  if (any(is.na(cluster_assignments_df$cluster))) {
    names(cluster_colours)[length(cluster_colours)] <- NA
    cluster_colours[length(cluster_colours)] <- "grey"
  }
  
  # prepare additional data for margin visualisation
  bmi_data <- data$imputed_data %>% 
    group_by(ID) %>% 
    summarise(BMI = median(BMI)) %>% 
    dplyr::select(ID, BMI) %>% 
    dplyr::mutate(BMI = cut(BMI, breaks = c(16, 18.5, 25, 30, 35, 40, Inf), right = FALSE))
    
  col_add_data <- data$imputed_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(.imp == 1) %>% 
    left_join(data$age, by = "ID") %>% 
    dplyr::mutate(age = cut(age, breaks = c(18, 25, 40, 55, Inf))) %>% 
    dplyr::select(ID, age, sex)
  col_info <- cytokine_data_median_wide %>% 
    dplyr::mutate(# check for NAs in cluster assignment
      cluster = if_else(cluster == "cluster NA", NA_character_, cluster),
      medication = if_else(grepl("^F", ID), "mostly not medicated", "mostly medicated")) %>%
    left_join(col_add_data, by = "ID") %>% 
    left_join(data$proband_status, by = "ID") %>% 
    dplyr::rename(`diagnosis status` = status) %>% 
    left_join(data$cidi_depression, by = "ID") %>% 
    dplyr::mutate(`current dep.` = if_else(depression_current, "yes", "no")) %>% 
    left_join(data$medication_data, by = "ID") %>% 
    dplyr::mutate(psychopharmaca = if_else(t0_any_psych_medication, "yes", "no")) %>% 
    left_join(bmi_data, by = "ID")
  
  diagnosis_status_colours <- c(`current case` = "#1B9E77")
  if ("control" %in% col_info$`diagnosis status`) {
    diagnosis_status_colours <- c(diagnosis_status_colours,
                                  `control` = "#D95F02")
  }
  if ("lifetime case" %in% col_info$`diagnosis status`) {
    diagnosis_status_colours <- c(diagnosis_status_colours,
                                  `lifetime case` = "#E6AB02")
  }
  
  heatmap_ann_colors <- list(
    `diagnosis status` = diagnosis_status_colours,
    `current dep.` = c(no = "#FEAF16", yes = "#B00068"),
    sex = c(male = "#7570B3", female = "#E7298A"),
    medication = c(`mostly medicated` = "#595959", `mostly not medicated` = "#bbbbbb"),
    study = c(Optima = "#595959", Become = "#bbbbbb"),
    psychopharmaca = c(yes = "#E41A1C", no = "#377EB8"),
    cluster = cluster_colours,
    `cluster cytokines + RNA` = cluster_colours,
    age = c(`(18,25]` = "#b3d1ff", `(25,40]` = "#66a3ff", `(40,55]` ="#0066ff",
            `(55,Inf]` = "#0047b3"),
    BMI = c(`[16,18.5)` = "#D8F0D8", `[18.5,25)` = "#92D192", `[25,30)` = "#59B359",
            `[30,35)` = "#2C942C", `[35,40)` = "#0C750C", `[40,Inf)` = "#005C00"),
    `panel name` = c(`Angiogenesis Panel` = viridis(7)[7],
                     `Chemokine Panel` =  viridis(7)[6],
                     `Cytokine Panel 1` = viridis(7)[5],
                     `Cytokine Panel 2` = viridis(7)[4],
                     `Proinflammatory Panel` = viridis(7)[3],
                     `TH17 Panel` = viridis(7)[2],
                     `Vascular Injury Panel` = viridis(7)[1])
  )
  
  if ("previous_cluster" %in% colnames(cytokine_data_median_wide)) {
    col_info <- col_info %>% 
      dplyr::mutate(`prev. cluster` = previous_cluster,
                    # check for NAs in cluster assignment
                    `prev. cluster` = if_else(`prev. cluster` == "cluster NA", NA_character_, `prev. cluster`)) %>% 
      dplyr::select(cluster, `prev. cluster`, medication, psychopharmaca,
                    age, sex, BMI, `current dep.`, `diagnosis status`, ID)
    
    temp_previous_cluster_colours <- previous_cluster_colours[seq_len(length(unique(previous_assignments_df$previous_cluster)))]
    
    if (any(is.na(previous_assignments_df$previous_cluster))) {
      names(temp_previous_cluster_colours)[length(temp_previous_cluster_colours)] <- NA
      temp_previous_cluster_colours[length(temp_previous_cluster_colours)] <- "grey"
    }
    
    heatmap_ann_colors[["prev. cluster"]] <- temp_previous_cluster_colours
    heatmap_ann_colors[["cluster cytokines"]] <- temp_previous_cluster_colours
    
  } else {
    col_info <- col_info %>% 
      dplyr::select(cluster, medication, psychopharmaca,
                    age, sex, BMI, `current dep.`, `diagnosis status`, ID)
  }
  
  if (!show_cases_controls) {
    col_info <- col_info %>% 
      dplyr::select(-`diagnosis status`)
  }
  
  if (!show_BMI) {
    col_info <- col_info %>% 
      dplyr::select(-BMI)
  }
  
  if (!show_psychopharmaca) {
    col_info <- col_info %>% 
      dplyr::select(-psychopharmaca)
  }
  
  if (var_name_study_instead_med) {
    col_info <- col_info %>% 
      dplyr::rename(study = medication) %>% 
      dplyr::mutate(study = case_when(
        study == "mostly medicated" ~ "Optima",
        study == "mostly not medicated" ~ "Become",
        TRUE ~ study
      ))
  }
  
  if (sumo_cytokine_names) {
    col_info <- col_info %>% 
      dplyr::rename(`cluster cytokines` = `prev. cluster`,
                    `cluster cytokines + RNA` = cluster)
  }
  
  if (!is.null(included_annotations)) {
    col_info <- col_info %>% 
      dplyr::select(all_of(included_annotations), ID)
  }
  
  col_id <- col_info$ID
  col_info <- col_info %>% 
    dplyr::select(-ID)
  col_info <- as.data.frame(col_info)
  
  rownames(col_info) <- col_id
  
  if (show_cytokine_panel_names) {
    cytokine_id <- panel_name_annotation$cytokine
    panel_name_annotation <- panel_name_annotation %>% 
      dplyr::select(-cytokine)
    row_info <- as.data.frame(panel_name_annotation)
    rownames(row_info) <- cytokine_id
  } else {
    row_info <- NA
  }
  
  pheatmap(heatmap_data,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           show_colnames = FALSE,
           annotation_col = col_info,
           annotation_row = row_info,
           annotation_colors = heatmap_ann_colors,
           annotation_legend = annotation_legend,
           annotation_names_row = FALSE,
           legend = legend,
           na_col = "green",
           ...)
  
}

plot_age <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    compare_means_test = set_defaults()$compare_means_test,
    compare_means_vjust = 1,
    compare_means_hjust = 0,
    x_axis_legend = element_text(),
    ...
) {
  plot_data <- data$age %>% 
    right_join(cluster_assignments_df, by = "ID")
  compare_data <- plot_data %>% 
    filter(cluster != "controls", cluster != "lifetime cases")
  p <- plot_data %>% 
    ggplot(aes(x = cluster, y = age)) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(aes(group = cluster), width = 0.2) +
    stat_compare_means(
      data = compare_data,
      method = compare_means_test,
      vjust = compare_means_vjust,
      hjust = compare_means_hjust) + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title.x = x_axis_legend) +
    scale_fill_manual(values = cluster_colours)
  
  add_dashes(p, plot_data)
}

plot_BMI <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    compare_means_test = set_defaults()$compare_means_test,
    compare_means_vjust = 1,
    compare_means_hjust = 0,
    median_BMI_annotation = TRUE,
    x_axis_legend = element_text(),
    use_unimputed_data = FALSE,
    ...
) {
  bmi_var_name <- "median_BMI"
  if (use_unimputed_data) {
    plot_data <- data$cytokine_uncorrected %>% 
      right_join(cluster_assignments_df, by = "ID")
    bmi_var_name <- "BMI"
  } else {
    plot_data <- data$imputed_data %>% 
      group_by(ID) %>% 
      summarise(median_BMI = median(BMI)) %>% 
      right_join(cluster_assignments_df, by = "ID")
    
    if(!median_BMI_annotation) {
      plot_data <- plot_data %>% 
        dplyr::rename(BMI = median_BMI)
      bmi_var_name <- "BMI"
    }
  }
  
  
  compare_data <- plot_data %>% 
    filter(cluster != "controls", cluster != "lifetime cases")
  p <- plot_data %>% 
    ggplot(aes(x = cluster, y = .data[[bmi_var_name]])) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(aes(group = cluster), width = 0.2) +
    stat_compare_means(
      data = compare_data,
      method = compare_means_test,
      vjust = compare_means_vjust,
      hjust = compare_means_hjust) + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title.x = x_axis_legend
      ) +
    scale_fill_manual(values = cluster_colours)
  
  if(median_BMI_annotation) {
    p <- p +
      ylab("median BMI")
  }
  
  add_dashes(p, plot_data)
}

plot_sex <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    ...
) {
  data$imputed_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(.imp == 1) %>%
    ggplot(aes(x = cluster, fill = sex)) +
    geom_bar() +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("#E7298A", "#7570B3"))
}

plot_BDI <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    compare_means_test = set_defaults()$compare_means_test,
    compare_means_vjust = 1,
    compare_means_hjust = 0,
    x_axis_legend = element_text(),
    ...
) {
  plot_data <- data$bdi_sumscores %>% 
    group_by(ID) %>% 
    right_join(cluster_assignments_df, by = "ID")
  compare_data <- plot_data %>% 
    filter(cluster != "controls", cluster != "lifetime cases")
  p <- plot_data %>% 
    ggplot(aes(x = cluster, y = BDI_total)) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(aes(group = cluster), width = 0.2) +
    stat_compare_means(
      data = compare_data,
      method = compare_means_test,
      vjust = compare_means_vjust,
      hjust = compare_means_hjust) + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title.x = x_axis_legend
      ) +
    ylab("median BDI") +
    scale_fill_manual(values = cluster_colours)
  
  add_dashes(p, plot_data)
}

plot_celltypes <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    compare_means_vjust = 1,
    compare_means_hjust = 0,
    compare_means_test = set_defaults()$compare_means_test,
    celltype_type = c("RNA_updated_14"),
    ...
) {
  if (celltype_type == "RNA_updated_14") {
    plot_data <- data$celltype_data_rna_updated_14 %>% 
      pivot_longer(
        cols = -ID,
        names_to = "celltype",
        values_to = "percentage"
      ) %>% 
      dplyr::mutate(percentage = percentage / 100,
                    celltype = case_when(
                      celltype == "Dendritic cells activated" ~ "Dendritic cells\nactivated",
                      celltype == "T cells CD4 memory activated" ~ "T cells CD4 memory\nactivated",
                      celltype == "T cells CD4 memory resting" ~ "T cells CD4 memory\nresting",
                      TRUE ~ celltype
                    ))
  } else {
    stop("This celltype_type is not available.")
  }
  
  if ("celltype_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    celltype_index <- plot_data$celltype %in% match.call()[["celltype_list"]]
    plot_data <- plot_data[celltype_index, ]
  }
  
  plot_data <- plot_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!is.na(celltype))
  compare_data <- plot_data %>% 
    filter(cluster != "controls", cluster != "lifetime cases")
  p <- plot_data %>% 
    ggplot() +
    aes(x = cluster, y = percentage,) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(aes(group = cluster), width = 0.2) +
    stat_compare_means(
      data = compare_data,
      method = compare_means_test,
      vjust = compare_means_vjust,
      hjust = compare_means_hjust
    ) +
    theme_bw() +
    scale_fill_manual(values = cluster_colours) +
    facet_wrap(~celltype, scales = "free_y")
  
  add_dashes(p, plot_data)
}


plot_cytokines <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    violin = FALSE,
    pretty_names = FALSE,
    complete_names = FALSE,
    cytokine_version = c("corrected_1", "additional_unimputed"),
    ncol = NULL,
    legend.position = "right",
    sort_by_list = FALSE,
    correct_cytokine_names = FALSE,
    ...
) {
  if (cytokine_version == "corrected_1") {
    plot_data <- data$cytokine_corrected_1 %>% 
      pivot_longer(
        cols = bFGF_corrected_1:sIL_6R_corrected_1,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      ) %>% 
      dplyr::mutate(cytokine = str_extract(cytokine, ".+?(?=_corrected_1)"))
  } else if (cytokine_version == "additional_unimputed") {
    plot_data <- data$cytokine_uncorrected %>% 
      pivot_longer(
        cols = bFGF:sIL_6R,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      )
  } else {
    stop("This cytokine version does not exist.")
  }
  
  plot_data <- plot_data %>% 
    group_by(cytokine, ID) %>% 
    summarise(
      median_normalised_concentration = median(normalised_concentration,
                                               na.rm = TRUE)
    ) %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    arrange(cluster)
  
  if ("cytokine_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    cytokine_index <- plot_data$cytokine %in% match.call()[["cytokine_list"]]
    plot_data <- plot_data[cytokine_index, ]
    cytokine_order <- data.frame(cytokine = match.call()[["cytokine_list"]])
  }
  
  if (complete_names) {
    # only implemented for a subset so far
    plot_data <- plot_data %>% 
      dplyr::mutate(
        cytokine = str_replace(cytokine, "hsCRP", "high sensitivity C-reactive protein"),
        cytokine = str_replace(cytokine, "IL_1RA", "interleukin 1 receptor antagonist"),
        cytokine = str_replace(cytokine, "SAA", "serum amyloid A"),
        cytokine = str_replace(cytokine, "IL_7", "interleukin 7"),
        cytokine = str_replace(cytokine, "TARC", "chemokine (C-C motif) ligand 17\n(TARC)"),
        cytokine = str_replace(cytokine, "VEGF_A_HS", "vascular endothelial growth factor A\nhigh sensitivity"),
        cytokine = str_replace(cytokine, "VEGF_A_LS", "vascular endothelial growth factor A\nlow sensitivity"),
        cytokine = str_replace(cytokine, "VEGF_C", "vascular endothelial growth factor C"),
        cytokine = str_replace(cytokine, "IL_16", "interleukin 16"),
        cytokine = str_replace(cytokine, "IP_10", "C-X-C motif chemokine ligand 10\n(IP-10)"),
        cytokine = str_replace(cytokine, "MCP_1", "chemokine (C-C motif) ligand 2\n(MCP-1)"),
        cytokine = str_replace(cytokine, "MCP_4", "chemokine (C-C motif) ligand 13\n(MCP-4)"),
        cytokine = str_replace(cytokine, "MIP_1beta", "chemokine (C-C motif) ligand 4\n(MIP-1beta)"),
        cytokine = str_replace(cytokine, "MIP_1alpha", "chemokine (C-C motif) ligand 3\n(MIP-1beta)"),
        cytokine = str_replace(cytokine, "TNF_alpha", "tumor necrosis factor alpha"),
        cytokine = str_replace(cytokine, "MDC", "chemokine (C-C motif) ligand 22\n(MDC)"),
        cytokine = str_replace(cytokine, "IL_17A", "interleukin 17A"),
        cytokine = str_replace(cytokine, "IL_12_IL_23p40", "interleukin 12/23p40"),
        cytokine = str_replace(cytokine, "Eotaxin_3", "chemokine (C-C motif) ligand 26\n(Eotaxin 3)")
      )
    if (sort_by_list) {
      cytokine_order <- cytokine_order %>% 
        dplyr::mutate(
          cytokine = str_replace(cytokine, "hsCRP", "high sensitivity C-reactive protein"),
          cytokine = str_replace(cytokine, "IL_1RA", "interleukin 1 receptor antagonist"),
          cytokine = str_replace(cytokine, "SAA", "serum amyloid A"),
          cytokine = str_replace(cytokine, "IL_7", "interleukin 7"),
          cytokine = str_replace(cytokine, "TARC", "chemokine (C-C motif) ligand 17\n(TARC)"),
          cytokine = str_replace(cytokine, "VEGF_A_HS", "vascular endothelial growth factor A\nhigh sensitivity"),
          cytokine = str_replace(cytokine, "VEGF_A_LS", "vascular endothelial growth factor A\nlow sensitivity"),
          cytokine = str_replace(cytokine, "VEGF_C", "vascular endothelial growth factor C"),
          cytokine = str_replace(cytokine, "IL_16", "interleukin 16"),
          cytokine = str_replace(cytokine, "IP_10", "C-X-C motif chemokine ligand 10\n(IP-10)"),
          cytokine = str_replace(cytokine, "MCP_1", "chemokine (C-C motif) ligand 2\n(MCP-1)"),
          cytokine = str_replace(cytokine, "MCP_4", "chemokine (C-C motif) ligand 13\n(MCP-4)"),
          cytokine = str_replace(cytokine, "MIP_1beta", "chemokine (C-C motif) ligand 4\n(MIP-1beta)"),
          cytokine = str_replace(cytokine, "MIP_1alpha", "chemokine (C-C motif) ligand 3\n(MIP-1beta)"),
          cytokine = str_replace(cytokine, "TNF_alpha", "tumor necrosis factor alpha"),
          cytokine = str_replace(cytokine, "MDC", "chemokine (C-C motif) ligand 22\n(MDC)"),
          cytokine = str_replace(cytokine, "IL_17A", "interleukin 17A"),
          cytokine = str_replace(cytokine, "IL_12_IL_23p40", "interleukin 12/23p40"),
          cytokine = str_replace(cytokine, "Eotaxin_3", "chemokine (C-C motif) ligand 26\n(Eotaxin 3)")
        )
    }
    
  }
  
  if (correct_cytokine_names) {
    plot_data <- correct_cytokine_names(plot_data)
    
    if (sort_by_list) {
      cytokine_order <- correct_cytokine_names(cytokine_order)
    }
  }
  
  if (pretty_names) {
    plot_data <- pretty_cytokine_names(plot_data)
    
    if (sort_by_list) {
      cytokine_order <- pretty_cytokine_names(cytokine_order)
    }
  }
  
  if (sort_by_list) {
    cytokine_order_vector <- cytokine_order %>% 
      pull(cytokine)
    plot_data <- plot_data %>% 
      dplyr::mutate(cytokine = factor(cytokine, levels = cytokine_order_vector))
  }
  
  if (cytokine_version == "corrected_1") {
    ylab_string <- "normalised conc. corrected for age/sex/BMI"
  } else {
    ylab_string <- "median normalised conc."
  }
  
  if (violin) {
    p <- plot_data %>% 
      ggplot(aes(x = cluster, y = median_normalised_concentration)) +
      geom_violin(aes(fill = cluster)) +
      geom_boxplot(aes(group = cluster), width = 0.2) +
      theme_bw() +
      theme(legend.position = legend.position) +
      scale_fill_manual(values = cluster_colours) +
      facet_wrap(~cytokine, ncol = ncol)
  } else {
    p <- plot_data %>% 
      ggplot(aes(x = cluster, y = median_normalised_concentration, fill = cluster)) +
      geom_boxplot() +
      theme_bw() +
      theme(legend.position = legend.position) +
      scale_fill_manual(values = cluster_colours) +
      facet_wrap(~cytokine, ncol = ncol)
  }
  
  p <- add_dashes(p, plot_data)
  
  p +
    ylab(ylab_string)
  
}

plot_BDI_sleeping_new <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    bdi_colours = c("-3" = "#3182bd", "-2" = "#9ecae1", "-1" = "#deebf7",
                    "0" = "#f0f0f0", "1" = "#fee0d2", "2" = "#fc9272", "3" = "#de2d26"),
    ...
) {
  # sleeping
  data$bdi_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    dplyr::mutate(BDI_16 = as.factor(BDI_16)) %>%
    ggplot(aes(x = cluster, fill = BDI_16)) +
    geom_bar(position = "fill") +
    theme_bw() +
    ggtitle("sleeping") +
    scale_fill_manual(values = bdi_colours)
}

plot_BDI_eating_new <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    bdi_colours = c("-3" = "#3182bd", "-2" = "#9ecae1", "-1" = "#deebf7",
                    "0" = "#f0f0f0", "1" = "#fee0d2", "2" = "#fc9272", "3" = "#de2d26"),
    ...
) {
  # eating
  data$bdi_data %>%
    right_join(cluster_assignments_df, by = "ID") %>%
    dplyr::mutate(BDI_18 = as.factor(BDI_18)) %>%
    ggplot(aes(x = cluster, fill = BDI_18)) +
    geom_bar(position = "fill") +
    theme_bw() +
    ggtitle("eating") +
    scale_fill_manual(values = bdi_colours)
}

plot_variable_importance <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    omic_type = c("cytokines_RNA"),
    cytokine_version = c("corrected_1", "additional_unimputed"),
    rna_version = c("standard", "corrected_1"),
    celltype_version = c("standard", "updated"),
    correct_cytokine_names = FALSE,
    ...
) {
  
  if (omic_type == "cytokines_RNA") {
    plot_variable_importance_cytokines_rna(
      data = data,
      cluster_assignments_df = cluster_assignments_df,
      cluster_colours = cluster_colours,
      cytokine_version = cytokine_version,
      rna_version = rna_version,
      correct_cytokine_names = correct_cytokine_names,
      ...
    )
  } else {
    stop("This omic type is not implemented yet")
  }
}

calculate_variable_importance_cytokines_rna <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    cytokine_version = c("corrected_1", "additional_unimputed"),
    rna_version = c("standard", "corrected_1"),
    show_gene_names = FALSE,
    include_age_bmi = FALSE,
    include_celltypes = FALSE,
    correct_cytokine_names = FALSE
) {
  rename_correction <- function(name) {
    stringr::str_extract(name, ".+?(?=_corrected)")
  }
  rename_correction_2 <- function(name) {
    stringr::str_extract(name, ".+?(?=_median)")
  }
  if (cytokine_version == "corrected_1") {
    variable_importance_data <- data$cytokine_corrected_1 %>% 
      rename_with(rename_correction, bFGF_corrected_1:sIL_6R_corrected_1) %>% 
      group_by(ID) %>% 
      summarise(across(bFGF:sIL_6R, median, na.rm = TRUE, .names = "{.col}_median")) %>% 
      right_join(cluster_assignments_df, by = "ID") %>% 
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      dplyr::select(bFGF_median:sIL_6R_median, cluster, ID) %>% 
      dplyr::rename_with(rename_correction_2, bFGF_median:sIL_6R_median)
  } else if (cytokine_version == "additional_unimputed") {
    variable_importance_data <- data$cytokine_uncorrected %>% 
      group_by(ID) %>% 
      summarise(across(bFGF:sIL_6R, median, na.rm = TRUE, .names = "{.col}_median")) %>% 
      right_join(cluster_assignments_df, by = "ID") %>% 
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      dplyr::select(bFGF_median:sIL_6R_median, cluster, ID) %>% 
      dplyr::rename_with(rename_correction_2, bFGF_median:sIL_6R_median)
  } else {
    stop("Not implemented yet.")
  }
  
  if (rna_version == "standard") {
    variable_importance_data_2 <- data$rna %>% 
      pivot_wider(
        id_cols = ID,
        names_from = "gene",
        values_from = "expression"
      ) %>% 
      right_join(cluster_assignments_df, by = "ID") %>%
      filter(!cluster %in% c("controls", "lifetime cases"))
  } else if ((rna_version == "corrected_1")) {
    variable_importance_data_2 <- data$rna_corrected_1 %>% 
      pivot_wider(
        id_cols = ID,
        names_from = "gene",
        values_from = "expression"
      ) %>% 
      rename_with(rename_correction, starts_with("ENSG")) %>% 
      right_join(cluster_assignments_df, by = "ID") %>%
      filter(!cluster %in% c("controls", "lifetime cases"))
  } else {
    stop("Not implemented yet.")
  }
  
  if (include_age_bmi) {
    variable_importance_data_3 <- data$imputed_data %>% 
      group_by(ID) %>%
      # for the latest data, the BMI is complete and the summary step is
      # not necessarily needed
      summarise(BMI = median(BMI)) %>% 
      right_join(cluster_assignments_df, by = "ID") %>% 
      left_join(data$age, by = "ID") %>% 
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      dplyr::select(-cluster)
    
    variable_importance_data <- variable_importance_data %>% 
      full_join(variable_importance_data_2 %>% dplyr::select(-cluster), by = "ID") %>% 
      full_join(variable_importance_data_3, by = "ID")
  } else {
    variable_importance_data <- variable_importance_data %>% 
      full_join(variable_importance_data_2 %>% dplyr::select(-cluster), by = "ID")
  }
  
  if (include_celltypes) {
    variable_importance_data_4 <- data$celltype_data_rna_updated_14 %>% 
      right_join(cluster_assignments_df, by = "ID") %>%
      filter(!cluster %in% c("controls", "lifetime cases"))
    
    variable_importance_data <- variable_importance_data %>% 
      full_join(variable_importance_data_4 %>% dplyr::select(-cluster), by = "ID")
  }
  
  variable_importance_data <- variable_importance_data %>% 
    dplyr::select(-ID)
  
  variable_names <- colnames(variable_importance_data)
  variable_names <- variable_names[variable_names != "cluster"]
  
  variable_importance <- set_names(variable_names) %>% 
    map_dfr(function(x) {
      lm_res <- lm(variable_importance_data[[x]] ~ variable_importance_data[["cluster"]])
      anova_res <- anova(lm_res)
      data.frame(F_value = anova_res[["F value"]][1])
    }, .id = "variable")
  
  if (show_gene_names) {
    variable_importance <- variable_importance %>% 
      left_join(data$rna_info, by = c("variable" = "ensembl_gene_id")) %>% 
      dplyr::mutate(
        # replace one gene name that is now available
        external_gene_name = if_else(
          external_gene_name == "AC005747.1",
          "GAS7",
          external_gene_name
        ),
        variable = if_else(str_detect(variable, "^ENSG"), paste0(external_gene_name, " ", variable),
                                       variable))
  }
  
  if (correct_cytokine_names) {
    variable_importance <- correct_cytokine_names(
      variable_importance,
      var_name = variable
    )
  }
  
  variable_importance
}

plot_variable_importance_cytokines_rna <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    cytokine_version = c("corrected_1", "additional_unimputed"),
    rna_version = c("standard", "corrected_1"),
    number_variables_shown = 40,
    show_gene_names = FALSE,
    force_cytokines_in_plot = FALSE,
    include_age_bmi = FALSE,
    include_celltypes = FALSE,
    show_annotation_string = TRUE,
    show_gene_id = TRUE,
    # order: clinical, immune marker, RNAseq
    variable_type_colours = c("#F78154", "#9ec91d", "#4196A3"),
    legend_position = "none",
    correct_cytokine_names,
    
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
    include_celltypes = include_celltypes,
    correct_cytokine_names = correct_cytokine_names
  )
  
  importance_order <- variable_importance %>% 
    arrange(desc(F_value)) %>% 
    pull(variable)
  
  if (cytokine_version == "corrected_1") {
    annotation_string <- "Based on normalised values corrected for age/sex/BMI."
  } else {
    annotation_string <- "Based on normalised values."
  }
  
  if (!show_annotation_string) {
    annotation_string <- NULL
  }
  
  plot_data_all <- variable_importance %>% 
    arrange(desc(F_value)) %>% 
    dplyr::mutate(variable_sort = variable)
  
  if (force_cytokines_in_plot) {
    # take the 3/4th of the most important variables (not necessarily RNA)
    plot_data_rna <- plot_data_all %>% 
      slice(seq_len(floor(number_variables_shown * 3/4)))
    plot_data_cytokines <- plot_data_all %>% 
      filter(str_detect(variable_sort, "ENSG", negate = TRUE),
             !variable_sort %in% plot_data_rna$variable_sort) %>% 
      slice(seq_len(floor(number_variables_shown * 1/4)))
    plot_data <- bind_rows(plot_data_rna, plot_data_cytokines)
  } else {
    plot_data <- plot_data_all %>% 
      # only take the first 40 most important variables as default
      slice(seq_len(number_variables_shown))
  }
  
  if (!show_gene_id) {
    plot_data <- plot_data %>% 
      dplyr::mutate(variable = str_replace(variable, " ENSG.*$", ""))
  }
  
  # I need to extract the order and use it as the levels, because otherwise it
  # will be ordered alphabetically
  importance_order <- plot_data %>% 
    pull(variable)

  p <- plot_data %>% 
    dplyr::mutate(`variable type` = case_when(
      # the -1 removes the ID colname
      variable_sort %in% c("age", "BMI",
                           colnames(data$celltype_data_rna_updated_14)[-1]) ~ "biological",
      str_detect(variable_sort, "ENSG") ~ "gene",
      TRUE ~ "immune marker"
      ),
      `variable type` = factor(`variable type`, levels = c("biological",
                                                           "immune marker",
                                                           "gene"))
    ) %>% 
    # this assumes that all the gene names shown are unique (shouldn't be a
    # problem for 40 genes, can be a problem for the complete data set)
    dplyr::mutate(variable = factor(variable, levels = importance_order)) %>% 
    ggplot() +
    aes(x = variable, y = F_value, fill = `variable type`) +
    geom_col() +
    scale_fill_manual(values = variable_type_colours) +
    theme_bw() +
    theme(
      legend.position = legend_position,
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    ylab("variable importance (F value)") +
    labs(
      caption = annotation_string
    )
  
  if (force_cytokines_in_plot) {
    p +
      geom_vline(aes(xintercept = floor(number_variables_shown * 3/4) + 0.5),
                 linetype = "longdash",
                 colour = "#555555")
  } else {
    p
  }
}

plot_RNA <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    rna_type = c("standard", "corrected_1"),
    compare_means_vjust = 1,
    compare_means_hjust = 0,
    compare_means_test = set_defaults()$compare_means_test,
    sort_by_list = FALSE,
    ncol = NULL,
    facet_scale = "free_y",
    # if the gene IDs are not used, the gene names have to be unique
    show_gene_id = TRUE,
    ...
) {
  if (rna_type == "standard") {
    plot_data <- data$rna
    lab_string <- "normalised expression"
  } else if (rna_type == "corrected_1") {
    plot_data <- data$rna_corrected_1 %>% 
      dplyr::mutate(gene = str_extract(gene, ".+?(?=_corrected)"))
    lab_string <- "normalised expression corrected for age/sex/BMI"
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
  
  # add the gene name info
  plot_data <- plot_data %>% 
    left_join(data$rna_info, by = c("gene" = "ensembl_gene_id")) %>% 
    dplyr::mutate(gene = paste0(gene, "\n", external_gene_name))
  
  # always set up the gene order and only decide later if it is used
  gene_order <- data.frame(gene = match.call()[["gene_list"]]) %>% 
      left_join(data$rna_info, by = c("gene" = "ensembl_gene_id")) %>% 
      dplyr::mutate(gene = paste0(gene, "\n", external_gene_name))
  
  if (!show_gene_id) {
    plot_data <- plot_data %>% 
      dplyr::mutate(gene = str_replace(gene, "^.*\\\n", ""))
    gene_order <- gene_order %>% 
      dplyr::mutate(gene = str_replace(gene, "^.*\\\n", ""))
  }
  
  if (sort_by_list) {
    gene_order_vector <- gene_order %>% 
      pull(gene)
    plot_data <- plot_data %>% 
      dplyr::mutate(gene = factor(gene, levels = gene_order_vector))
  }
  
  plot_data <- plot_data %>% 
    left_join(cluster_assignments_df, by = "ID")
  compare_data <- plot_data %>% 
    filter(cluster != "controls", cluster != "lifetime cases")
  p <- plot_data %>% 
    ggplot(aes(x = cluster, y = expression)) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(aes(group = cluster), width = 0.2) +
    stat_compare_means(
      data = compare_data,
      method = compare_means_test,
      vjust = compare_means_vjust,
      hjust = compare_means_hjust
    ) +
    theme_bw() +
    theme_bw() +
    scale_fill_manual(values = cluster_colours) +
    facet_wrap(~gene, scale = facet_scale, ncol = ncol) +
    ylab(lab_string)
  
  add_dashes(p, plot_data)
}

plot_cidi_single <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    separate_by = c("cluster", "diagnosis"),
    cidi_type = c("all", "current", "lifetime"),
    remove_type_string = FALSE,
    ...
) {
  separate_by <- match.arg(separate_by)
  cidi_type <- match.arg(cidi_type)
  
  plot_data <- data$cidi_single %>% 
    right_join(cluster_assignments_df, by = "ID")
  
  # default: cidi_type == "all"
  if (cidi_type == "current") {
    plot_data <- plot_data %>% 
      dplyr::select(ID, cluster, contains("curr"))
  }
  
  if (cidi_type == "lifetime") {
    plot_data <- plot_data %>% 
      dplyr::select(ID, cluster, contains("lt"))
  }
  
  cluster_sizes <- cluster_assignments_df %>% 
    group_by(cluster) %>% 
    summarise(cluster_size = n())
  
  plot_data <- plot_data %>% 
    pivot_longer(
      cols = -c(ID, cluster),
      names_to = "diagnosis",
      values_to = "has_diagnosis"
    ) %>% 
    filter(has_diagnosis != 0) %>% 
    group_by(cluster, diagnosis) %>% 
    summarise(count = n()) %>% 
    left_join(cluster_sizes, by = "cluster") %>% 
    dplyr::mutate(fraction = count / cluster_size,
                  # correct typos
                  diagnosis = str_replace(diagnosis, "Agorophobia", "Agoraphobia"),
                  diagnosis = str_replace(diagnosis, "Socialphobia", "Social phobia")
                  )
  
  if (remove_type_string) {
    plot_data <- plot_data %>% 
      mutate(
        diagnosis = str_replace(diagnosis, "(_curr|_lt)$", ""),
        diagnosis = str_replace(diagnosis, "_", " ")
      )
  }
  
  if (separate_by == "cluster") {
    plot_data %>% 
      ggplot() +
      aes(x = diagnosis, y = fraction, fill = diagnosis, label = count) +
      geom_col() +
      geom_text(vjust = 0.5) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      facet_wrap(~cluster)
  } else {
    plot_data %>% 
      ggplot() +
      aes(x = cluster, y = fraction, fill = cluster, label = count) +
      geom_col() +
      geom_text(vjust = 0.5) +
      theme_bw() +
      scale_fill_manual(values = cluster_colours) +
      facet_wrap(~diagnosis)
  }
}


calculate_RNA_correlation <- function(
    data,
    cluster_assignments_df,
    rna_type = c("standard", "corrected_1"),
    sample_type = c("cases"),
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
  
  # now that I put this code into it's own calculate_ function so that I can
  # call it outside of the plot function, match.call()[["gene_list"]]
  # does not work anymore if I call calculate_ outside of the plot function
  # therefore, I've added eval
  if ("gene_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    rna_index <- plot_data$gene %in% eval(match.call()[["gene_list"]])
    plot_data <- plot_data[rna_index, ]
  }
  
  # add the gene name info
  plot_data <- plot_data %>% 
    left_join(data$rna_info, by = c("gene" = "ensembl_gene_id")) %>% 
    dplyr::mutate(gene = paste0(gene, "\n", external_gene_name))
  
  plot_data <- plot_data %>% 
    left_join(cluster_assignments_df, by = "ID")
  
  if (sample_type == "cases") {
    plot_data <- plot_data %>% 
      filter(cluster != "controls" & cluster != "lifetime cases" & !is.na(cluster))
  } else {
    stop("This sample_type is not implemented.")
  }
  
  plot_data <- plot_data %>% 
    dplyr::select(-c(gene, description, is_available, cluster)) %>% 
    pivot_wider(
      id_cols = ID,
      names_from = external_gene_name,
      values_from = expression
    ) %>% 
    dplyr::select(-ID)
  
  # there are 8 observations with missing RNAseq data (completely misssing),
  # remove them
  corr_data <- cor(plot_data, use = "complete.obs")
  
  corr_data
}

# currently, this function only uses the gene names which can lead to problems
# if there are duplications (then I would need to include the ENSEMBL ID again)
plot_RNA_correlation <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    rna_type = c("standard", "corrected_1"),
    sample_type = c("cases"),
    tl.pos = "td",
    tl.cex = 12,
    ...
) {
  corr_data <- calculate_RNA_correlation(
    data = data,
    cluster_assignments_df = cluster_assignments_df,
    rna_type = rna_type,
    sample_type = sample_type,
    ...
  )
  
  ggcorrplot(corr_data,
             type = "upper",
             outline.color = "white",
             method = "square",
             tl.cex = tl.cex)
}

calculate_cytokine_correlation <- function(
    data,
    cluster_assignments_df,
    cytokine_version = c("corrected_1", "additional_unimputed"),
    sample_type = c("cases"),
    pretty_names = FALSE,
    include_age_bmi = FALSE,
    ...
) {
  if (cytokine_version == "corrected_1") {
    plot_data <- data$cytokine_corrected_1 %>% 
      pivot_longer(
        cols = bFGF_corrected_1:sIL_6R_corrected_1,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      ) %>% 
      dplyr::mutate(cytokine = str_extract(cytokine, ".+?(?=_corrected_1)"))
  } else if (cytokine_version == "additional_unimputed") {
    plot_data <- data$cytokine_uncorrected %>% 
      pivot_longer(
        cols = bFGF:sIL_6R,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      )
  } else {
    stop("This cytokine version does not exist.")
  }
  
  if ("cytokine_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    # now that I put this code into it's own calculate_ function so that I can
    # call it outside of the plot function, match.call()[["cytokine_list"]]
    # does not work anymore if I call calculate_ outside of the plot function
    # therefore, I've added eval
    cytokine_index <- plot_data$cytokine %in% eval(match.call()[["cytokine_list"]])
    plot_data <- plot_data[cytokine_index, ]
  }
  
  plot_data <- plot_data %>% 
    left_join(cluster_assignments_df, by = "ID")
  
  if (sample_type == "cases") {
    plot_data <- plot_data %>% 
      filter(cluster != "controls" & cluster != "lifetime cases" & !is.na(cluster))
  } else {
    stop("This sample_type is not implemented.")
  }
  
  if (pretty_names) {
    plot_data <- pretty_cytokine_names(plot_data)
  }
  
  if (include_age_bmi) {
    id_cols <- c("ID", "age", "BMI")
  } else {
    id_cols <- c("ID")
  }
  
  plot_data <- plot_data %>% 
    dplyr::select(-c(cluster)) %>% 
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = cytokine,
      values_from = normalised_concentration
    ) %>% 
    dplyr::select(-ID)
  
  # for every combination of cytokines, use the subset of complete observations
  # for this combination
  corr_data <- cor(plot_data, use = "pairwise.complete.obs")
  
  corr_data
}

plot_cytokine_correlation <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    cytokine_version = c("corrected_1", "additional_unimputed"),
    sample_type = c("cases"),
    tl.pos = "td",
    tl.cex = 12,
    pretty_names = FALSE,
    include_age_bmi = FALSE,
    ...
) {
  corr_data <- calculate_cytokine_correlation(
    data = data,
    cluster_assignments_df = cluster_assignments_df,
    cytokine_version = cytokine_version,
    sample_type = sample_type,
    pretty_names = pretty_names,
    include_age_bmi = include_age_bmi,
    ...
  )
  
  ggcorrplot(corr_data,
             type = "upper",
             outline.color = "white",
             method = "square",
             tl.cex = tl.cex)
}

calculate_cytokine_RNA_correlation <- function(
    data,
    cluster_assignments_df,
    cytokine_version = c("corrected_1", "additional_unimputed"),
    rna_type = c("standard", "corrected_1"),
    sample_type = c("cases"),
    pretty_names = FALSE,
    correct_cytokine_names = FALSE,
    include_age_bmi = FALSE,
    ...
) {
  if (cytokine_version == "corrected_1") {
    plot_data <- data$cytokine_corrected_1 %>% 
      pivot_longer(
        cols = bFGF_corrected_1:sIL_6R_corrected_1,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      ) %>% 
      dplyr::mutate(cytokine = str_extract(cytokine, ".+?(?=_corrected_1)"))
  } else if (cytokine_version == "additional_unimputed") {
    plot_data <- data$cytokine_uncorrected %>% 
      pivot_longer(
        cols = bFGF:sIL_6R,
        names_to = "cytokine",
        values_to = "normalised_concentration"
      )
  } else {
    stop("This cytokine version does not exist.")
  }
  
  if ("cytokine_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    # now that I put this code into it's own calculate_ function so that I can
    # call it outside of the plot function, match.call()[["cytokine_list"]]
    # does not work anymore if I call calculate_ outside of the plot function
    # therefore, I've added eval
    cytokine_index <- plot_data$cytokine %in% eval(match.call()[["cytokine_list"]])
    plot_data <- plot_data[cytokine_index, ]
  }
  
  plot_data <- plot_data %>% 
    left_join(cluster_assignments_df, by = "ID")
  
  if (rna_type == "standard") {
    plot_data_rna <- data$rna
  } else if (rna_type == "corrected_1") {
    plot_data_rna <- data$rna_corrected_1 %>% 
      dplyr::mutate(gene = str_extract(gene, ".+?(?=_corrected)"))
  } else {
    stop("This rna_type is not supported.")
  }
  
  # now that I put this code into it's own calculate_ function so that I can
  # call it outside of the plot function, match.call()[["gene_list"]]
  # does not work anymore if I call calculate_ outside of the plot function
  # therefore, I've added eval
  if ("gene_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    rna_index <- plot_data_rna$gene %in% eval(match.call()[["gene_list"]])
    plot_data_rna <- plot_data_rna[rna_index, ]
  }
  
  # add the gene name info
  plot_data_rna <- plot_data_rna %>% 
    left_join(data$rna_info, by = c("gene" = "ensembl_gene_id")) %>% 
    dplyr::mutate(gene = paste0(gene, "\n", external_gene_name))
  
  plot_data_rna <- plot_data_rna %>% 
    left_join(cluster_assignments_df, by = "ID")
  
  
  if (sample_type == "cases") {
    plot_data <- plot_data %>% 
      filter(cluster != "controls" & cluster != "lifetime cases" & !is.na(cluster))
    plot_data_rna <- plot_data_rna %>% 
      filter(cluster != "controls" & cluster != "lifetime cases" & !is.na(cluster))
  } else {
    stop("This sample_type is not implemented.")
  }
  
  if (pretty_names) {
    plot_data <- pretty_cytokine_names(plot_data)
  }
  if (correct_cytokine_names) {
    plot_data <- correct_cytokine_names(plot_data)
  }
  
  plot_data <- plot_data %>% 
    dplyr::select(-c(cluster)) %>% 
    pivot_wider(
      id_cols = ID,
      names_from = cytokine,
      values_from = normalised_concentration
    )
  
  plot_data_rna <- plot_data_rna %>% 
    dplyr::select(-c(gene, description, is_available, cluster)) %>% 
    pivot_wider(
      id_cols = ID,
      names_from = external_gene_name,
      values_from = expression
    )
  
  if (include_age_bmi) {
    plot_data_pheno <- data$cytokine_uncorrected %>%
      dplyr::select(ID, age, BMI)
    
    plot_data <- plot_data %>% 
      left_join(plot_data_pheno, by = "ID")
  }
  
  plot_data_all <- plot_data %>% 
    inner_join(plot_data_rna, by = "ID") %>% 
    dplyr::select(-ID)
  
  # for every combination of cytokines and genes, use the subset of complete observations
  # for this combination
  corr_data <- cor(plot_data_all, use = "pairwise.complete.obs")
  
  corr_data
}

plot_cytokine_RNA_correlation <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    cytokine_version = "additional_unimputed",
    rna_type = "standard",
    sample_type = c("cases"),
    tl.pos = "td",
    tl.cex = 12,
    pretty_names = FALSE,
    correct_cytokine_names = FALSE,
    ...
) {
  corr_data <- calculate_cytokine_RNA_correlation(
    data = data,
    cluster_assignments_df = cluster_assignments_df,
    cytokine_version = cytokine_version,
    rna_type = rna_type,
    sample_type = sample_type,
    pretty_names = pretty_names,
    correct_cytokine_names = correct_cytokine_names,
    ...
  )
  
  ggcorrplot(corr_data,
             type = "upper",
             outline.color = "white",
             method = "square",
             tl.cex = tl.cex)
}


plot_heatmap_genesets <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    previous_cluster_colours,
    id_order,
    previous_assignments = NULL,
    annotation_legend = TRUE,
    legend = TRUE,
    cluster_cols = FALSE,
    show_cases_controls = FALSE,
    show_current_depression = TRUE,
    show_sex = TRUE,
    show_age = TRUE,
    show_psychopharmaca = TRUE,
    show_medication = TRUE,
    scale_expression = FALSE,
    colour_palette = colorRampPalette(c("white", "red"))(100),
    geneset_data,
    cluster_rows = TRUE,
    geneset_list,
    ...
) {
  
  if (!is.null(previous_assignments)) {
    # have the possibility to use a data.frame for previous_assignments
    # (which then can also have a different order) than using cluster and
    # id_order vectors separately
    if (inherits(previous_assignments, "data.frame")) {
      previous_assignments_clean <- get_cluster_assignments(
        cluster_assignments = previous_assignments,
        add_healthy_controls = FALSE
      )
      
      previous_assignments_df <- previous_assignments_clean[["cluster_assignments_df"]]
      previous_assignments_df <- previous_assignments_df %>% 
        dplyr::rename(previous_cluster = cluster)
    } else {
      # have the ability to pass different ID order for previous assignments
      # because the data with standard cytokines was differently ordered than
      # the data with additional cytokines
      if (is.null(previous_assignments_order)) {
        old_id_order <- id_order
      } else {
        old_id_order <- previous_assignments_order
      }
      
      previous_assignments_df <- data.frame(
        ID = old_id_order,
        previous_cluster = previous_assignments
      )
    }
    
    cluster_assignments_df <- cluster_assignments_df %>% 
      # only use ids for clustering, no healthy controls or lifetime cases
      filter(cluster != "controls", cluster != "lifetime cases") %>% 
      left_join(previous_assignments_df, by = "ID")
  } else {
    # exclude healthy controls from cluster assignments
    # (needs to be done additionally to removing them from the cytokine values,
    # because otherwise these IDs are added by the right_join with the
    # cluster assignments later)
    cluster_assignments_df <- cluster_assignments_df %>% 
      # only use ids for clustering, no healthy controls
      filter(cluster != "controls", cluster != "lifetime cases")
  }
  plot_data <- geneset_data
  lab_string <- "normalised conc. corrected for age/sex/BMI"
  
  if (!"geneset_list" %in% names(match.call())) {
    stop("Please provide a `geneset_list` for which genes the expression should be shown.")
  }
  
  # somehow dplyr has problems filtering
  # -> should be resolved now because I made geneset_list a mandatory argument,
  # however I leave the code as is
  geneset_index <- colnames(plot_data) %in% c(match.call()[["geneset_list"]], "ID")
  plot_data <- plot_data[, geneset_index]
  plot_data <- plot_data %>% 
    relocate(all_of(geneset_list))
  
  plot_data <- plot_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    arrange(cluster, ID)
  
  if (scale_expression) {
    plot_data <- plot_data %>% 
      dplyr::mutate(across(all_of(match.call()[["geneset_list"]]), scale))
  }
  
  
  if ("previous_cluster" %in% colnames(plot_data)) {
    heatmap_data <- t(as.matrix(plot_data %>%
                                  dplyr::select(-c(ID, cluster, previous_cluster))))
  } else {
    heatmap_data <- t(as.matrix(plot_data %>% dplyr::select(-c(ID, cluster))))
  }
  
  colnames(heatmap_data) <- plot_data$ID
  
  # remove the "controls" and "lifetime" colour from the cluster colours
  if ("controls" %in% names(cluster_colours)) {
    cluster_colours <- cluster_colours[names(cluster_colours) != "controls"]
  }
  if ("lifetime cases" %in% names(cluster_colours)) {
    cluster_colours <- cluster_colours[names(cluster_colours) != "lifetime cases"]
  }
  
  # if there are NAs in the clustering assignment, set the last cluster colour
  # for NA to grey
  if (any(is.na(cluster_assignments_df$cluster))) {
    names(cluster_colours)[length(cluster_colours)] <- NA
    cluster_colours[length(cluster_colours)] <- "grey"
  }
  
  # prepare additional data for margin visualisation
  col_add_data <- data$imputed_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(.imp == 1) %>% 
    left_join(data$age, by = "ID") %>% 
    dplyr::mutate(age = cut(age, breaks = c(18, 25, 40, 55, Inf))) %>% 
    dplyr::select(ID, age, sex)
  
  # don't check for NAs in cluster assignment when it is a factor,
  # because it can't handle the factor object well
  if (!is.factor(cluster_assignments_df$cluster)) {
    col_info <- plot_data %>% 
      dplyr::mutate(# check for NAs in cluster assignment
        cluster = if_else(cluster == "cluster NA", NA_character_, cluster))
  } else {
    col_info <- plot_data
  }
  col_info <- col_info %>% 
    dplyr::mutate(
      medication = if_else(grepl("^F", ID), "mostly not medicated", "mostly medicated")) %>%
    left_join(col_add_data, by = "ID") %>% 
    left_join(data$cidi_data, by = "ID") %>% 
    left_join(data$proband_status, by = "ID") %>% 
    dplyr::rename(`diagnosis status` = status) %>% 
    left_join(data$cidi_depression, by = "ID") %>% 
    dplyr::mutate(`current dep.` = if_else(depression_current, "yes", "no")) %>% 
    left_join(data$medication_data, by = "ID") %>% 
    dplyr::mutate(psychopharmaca = if_else(t0_any_psych_medication, "yes", "no"))
  
  diagnosis_status_colours <- c(`current case` = "#1B9E77")
  if ("control" %in% col_info$`diagnosis status`) {
    diagnosis_status_colours <- c(diagnosis_status_colours,
                                  `control` = "#D95F02")
  }
  if ("lifetime case" %in% col_info$`diagnosis status`) {
    diagnosis_status_colours <- c(diagnosis_status_colours,
                                  `lifetime case` = "#E6AB02")
  }
  
  heatmap_ann_colors <- list(
    `diagnosis status` = diagnosis_status_colours,
    `current dep.` = c(no = "#FEAF16", yes = "#B00068"),
    sex = c(male = "#7570B3", female = "#E7298A"),
    medication = c(`mostly medicated` = "#66A61E", `mostly not medicated` = "#E6AB02"),
    psychopharmaca = c(yes = "#E41A1C", no = "#377EB8"),
    cluster = cluster_colours,
    age = c(`(18,25]` = "#b3d1ff", `(25,40]` = "#66a3ff", `(40,55]` ="#0066ff",
            `(55,Inf]` = "#0047b3")
  )
  
  if ("previous_cluster" %in% colnames(plot_data)) {
    col_info <- col_info %>% 
      dplyr::mutate(`prev. cluster` = previous_cluster,
                    # check for NAs in cluster assignment
                    `prev. cluster` = if_else(`prev. cluster` == "cluster NA", NA_character_, `prev. cluster`)) %>% 
      dplyr::select(cluster, `prev. cluster`, medication, psychopharmaca,
                    age, sex, `current dep.`, `diagnosis status`)
    
    temp_previous_cluster_colours <- previous_cluster_colours[seq_len(length(unique(previous_assignments_df$previous_cluster)))]
    
    if (any(is.na(previous_assignments_df$previous_cluster))) {
      names(temp_previous_cluster_colours)[length(temp_previous_cluster_colours)] <- NA
      temp_previous_cluster_colours[length(temp_previous_cluster_colours)] <- "grey"
    }
    
    heatmap_ann_colors[["prev. cluster"]] <- temp_previous_cluster_colours
    
  } else {
    col_info <- col_info %>% 
      dplyr::select(cluster, medication, psychopharmaca,
                    age, sex, `current dep.`, `diagnosis status`)
  }
  
  if (!show_cases_controls) {
    col_info <- col_info %>% 
      dplyr::select(-`diagnosis status`)
  }
  
  if (!show_current_depression) {
    col_info <- col_info %>% 
      dplyr::select(-`current dep.`)
  }
  
  if (!show_sex) {
    col_info <- col_info %>% 
      dplyr::select(-sex)
  }
  
  if (!show_age) {
    col_info <- col_info %>% 
      dplyr::select(-age)
  }
  
  if (!show_psychopharmaca) {
    col_info <- col_info %>% 
      dplyr::select(-psychopharmaca)
  }
  
  if (!show_medication) {
    col_info <- col_info %>% 
      dplyr::select(-medication)
  }
  
  
  col_info <- as.data.frame(col_info)
  
  rownames(col_info) <- plot_data$ID
  
  pheatmap(heatmap_data,
           color = colour_palette,
           cluster_rows = cluster_rows,
           cluster_cols = cluster_cols,
           show_colnames = FALSE,
           annotation_col = col_info,
           annotation_colors = heatmap_ann_colors,
           annotation_legend = annotation_legend,
           legend = legend,
           ...)
  
}

plot_imaging_residuals <- function(
    data,
    cluster_assignments_df,
    cluster_colours,
    compare_means_test = set_defaults()$compare_means_test,
    compare_means_vjust = 1,
    compare_means_hjust = 0,
    remove_controls = FALSE,
    show_n = TRUE,
    ...
) {
  plot_data <- data$imaging_residuals_ventricles %>% 
    right_join(cluster_assignments_df, by = "ID")
  
  size_data_detailed <- plot_data %>% 
    dplyr::select(ID, Added_to_4, cluster) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    # filter out cases for which we don't have information
    filter(!is.na(Added_to_4)) %>% 
    group_by(cluster) %>% 
    summarise(n = n())
  
  size_data_all <- plot_data %>% 
    dplyr::select(ID, cluster) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    group_by(cluster) %>% 
    summarise(n_all = n())
  
  size_data_detailed <- size_data_detailed %>% 
    left_join(size_data_all, by = "cluster")
  
  n_detailed <- purrr::pmap(size_data_detailed,
                            ~paste0(..1, " = ", ..2, " (", ..3, ")"))
  
  if (show_n) {
    n_string <- paste0("n (all): ",
                       paste0(n_detailed, collapse = "; "))
  } else {
    n_string <- NULL
  }
  
  if (remove_controls) {
    cluster_colours <- cluster_colours[names(cluster_colours) != "controls"]
    plot_data <- plot_data %>% 
      dplyr::filter(cluster != "controls") %>% 
      dplyr::mutate(cluster = forcats::fct_drop(cluster))
  }
  
  
  y_label <- "lateral ventricle volume without the\neffect of age, sex, ICV"
  
  compare_data <- plot_data %>% 
    filter(cluster != "controls", cluster != "lifetime cases")
  p <- plot_data %>% 
    ggplot(aes(x = cluster, y = Added_to_4)) +
    geom_violin(aes(fill = cluster)) +
    geom_boxplot(aes(group = cluster), width = 0.2) +
    stat_compare_means(
      data = compare_data,
      method = compare_means_test,
      vjust = compare_means_vjust,
      hjust = compare_means_hjust) + 
    theme_bw() +
    scale_fill_manual(values = cluster_colours) +
    labs(caption = n_string,
         y = y_label)
  
  if (remove_controls) {
    p
  } else {
    add_dashes(p, plot_data)
  }
}

plot_chromosomes <- function(
    data,
    cluster_assignments_df,
    drop_exotic_chromosomes = TRUE,
    ...
) {
  
  plot_data <- data$rna_info %>% 
    rename(chromosome = chromosome_name)
  
  if ("gene_list" %in% names(match.call())) {
    # somehow dplyr has problems filtering
    rna_index <- plot_data$ensembl_gene_id %in% match.call()[["gene_list"]]
    plot_data <- plot_data[rna_index, ]
  } else {
    stop("Please provide a `gene_list` for which genes the expression should be shown.")
  }
  
  if (drop_exotic_chromosomes) {
    # except 1-22 and X, drop all unused factor levels
    needed_extra_levels <- setdiff(
      c(unique(plot_data$chromosome), as.character(1:22), "X"),
      c(as.character(1:22), "X")
    )
    
    levels_for_plot <- c(as.character(1:22), "X", needed_extra_levels)
    plot_data <- plot_data %>% 
      mutate(
        chromosome = factor(chromosome, levels = levels_for_plot)
      )
  } else {
    extra_levels <- setdiff(
      unique(data$rna_info$chromosome_name),
      c(as.character(1:22), "X")
    )
    levels_for_plot <- c(as.character(1:22), "X", extra_levels)
    plot_data <- plot_data %>% 
      mutate(
        chromosome = factor(chromosome, levels = levels_for_plot)
      )
  }
  
  plot_data %>% 
    ggplot(aes(x = chromosome, fill = chromosome)) +
    geom_bar() +
    theme_bw() +
    scale_y_continuous(
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
      ) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(
      values = rep_len(palette.colors(palette = "Paired")[1:2],
                       length(levels(plot_data$chromosome))),
      drop = FALSE,
      guide = "none"
    )
}
