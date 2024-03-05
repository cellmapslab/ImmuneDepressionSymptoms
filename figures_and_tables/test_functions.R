library(dplyr)
library(broom)
library(purrr)
library(tidyr)
library(conover.test)

generate_test <- function(
    type = c("BDI_single"),
    cluster_assignments,
    general_data,
    ...
) {
  
  dots <- list(...)
  
  clean_cluster_assignments <-
    get_cluster_assignments(
      cluster_assignments = cluster_assignments,
      general_data = general_data
    )
  
  do.call(paste0("test_", type), c(list(data = general_data,
                                        cluster_assignments_df = clean_cluster_assignments$cluster_assignments_df),
                                   dots))
}

test_BDI_single <- function(
    data,
    cluster_assignments_df,
    ...
) {
  bdi_data_raw <- data$bdi_data %>% 
    left_join(cluster_assignments_df, by = "ID")
  var_name_bdi <- colnames(bdi_data_raw)
  var_name_bdi_index_exclude <- var_name_bdi %in% c("ID", "cluster")
  var_name_bdi <- var_name_bdi[!var_name_bdi_index_exclude]
  
  bdi_fisher_results <- purrr::set_names(var_name_bdi) %>% 
    purrr::map_dfr(function(curr_var_name) {
      data <- table(bdi_data_raw[[curr_var_name]], bdi_data_raw[["cluster"]])
      if (ncol(data) == 2) {
        res <- fisher.test(data, simulate.p.value = TRUE)
        ret <- data.frame(p_value = res[["p.value"]],
                          estimate = ifelse(is.null(res[["estimate"]]),
                                            NA, res[["estimate"]]))
        rownames(ret) <- NULL
        ret
      } else {
        data.frame(p_value = NA, estimate = NA)
      }
      
    }, .id = "disease")
  
  bdi_fisher_results %>% 
    filter(p_value < 0.05) %>% 
    arrange(p_value)
}

test_BDI <- function(
    data,
    cluster_assignments_df,
    method = c("conover", "Tukey"),
    ...
) {
  test_data <- data$bdi_sumscores %>% 
    group_by(ID) %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    mutate(cluster = paste0("cl. ", cluster))
  
  if (method == "conover") {
    conover.test(
      x = test_data$BDI_total,
      g = test_data$cluster,
      method = "bh",
      list = TRUE
    )
  } else if (method == "Tukey") {
    model <- aov(BDI_total ~ cluster, data = test_data)
    print(summary(model))
    print(TukeyHSD(model))
  } else {
    stop("This method is not available.")
  }
}

test_MADRS <- function(
    data,
    cluster_assignments_df,
    method = c("conover", "Tukey"),
    ...
) {
  test_data <- data$madrs %>% 
    group_by(ID) %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    mutate(cluster = paste0("cl. ", cluster))
  
  if (method == "conover") {
    conover.test(
      x = test_data$MADRS_total,
      g = test_data$cluster,
      method = "bh",
      list = TRUE
    )
  } else if (method == "Tukey") {
    model <- aov(MADRS_total ~ cluster, data = test_data)
    print(summary(model))
    print(TukeyHSD(model))
  } else {
    stop("This method is not available.")
  }
}

test_age <- function(
    data,
    cluster_assignments_df,
    method = c("conover", "Tukey"),
    ...
) {
  test_data <- data$age %>% 
    group_by(ID) %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    mutate(cluster = paste0("cl. ", cluster))
  
  if (method == "conover") {
    conover.test(
      x = test_data$age,
      g = test_data$cluster,
      method = "bh",
      list = TRUE
    )
  } else if (method == "Tukey") {
    model <- aov(age ~ cluster, data = test_data)
    print(summary(model))
    print(TukeyHSD(model))
  } else {
    stop("This method is not implemented.")
  }
}

test_BMI <- function(
    data,
    cluster_assignments_df,
    method = c("conover", "Tukey"),
    bmi_source = c("imputed", "with_missing"),
    ...
) {
  if (bmi_source == "imputed") {
    test_data <- data$imputed_data %>% 
      group_by(ID) %>% 
      summarise(median_BMI = median(BMI)) %>% 
      right_join(cluster_assignments_df, by = "ID") %>% 
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      mutate(cluster = paste0("cl. ", cluster))
  } else if (bmi_source == "with_missing") {
    test_data <- data$bmi %>% 
      # do the renaming so that it is consistent with the above test dataA
      rename(median_BMI = BMI) %>% 
      right_join(cluster_assignments_df, by = "ID") %>% 
      filter(!cluster %in% c("controls", "lifetime cases")) %>% 
      mutate(cluster = paste0("cl. ", cluster))
  } else {
    stop("This source is not available.")
  }
  
  if (method == "conover") {
    conover.test(
      x = test_data$median_BMI,
      g = test_data$cluster,
      method = "bh",
      list = TRUE
    )
  } else if (method == "Tukey") {
    model <- aov(median_BMI ~ cluster, data = test_data)
    print(summary(model))
    print(TukeyHSD(model))
  } else {
    stop("This method is not implemented.")
  }
}

test_cidi_diagnosis <- function(
    data,
    cluster_assignments_df,
    ...
) {
  test_data_raw <- data$cidi_single %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    mutate(cluster = paste0("cluster ", cluster)) %>% 
    pivot_longer(
      cols = -c(ID, cluster),
      names_to = "diagnosis",
      values_to = "has_diagnosis"
    ) %>% 
    filter(str_detect(diagnosis, "curr$"))
  
  diagnosis_values <- unique(test_data_raw$diagnosis)
  
  set_names(diagnosis_values) %>% 
    map_dfr(function(x) {
      data_single <- test_data_raw %>% 
        filter(diagnosis == x)
      
      number_cases <- data_single %>% 
        filter(has_diagnosis == 1) %>% 
        nrow()
      res <- data.frame(p_value = NA)
      if (number_cases > 0) {
        test_data <- with(data_single, table(has_diagnosis, cluster))
        res$p_value <- fisher.test(test_data, simulate.p.value = TRUE)$p.value
      }
      res
    }, .id = "diagnosis")
  
}

test_cidi_depression_dysthymia <- function(
    data,
    cluster_assignments_df,
    ...
) {
  test_data_raw <- data$cidi_single %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    mutate(cluster = paste0("cluster ", cluster),
           has_depression_dysthymia = Depression_full_curr == 1 |
             Dysthymia_full_curr == 1 | Depression_subthr_curr == 1) %>% 
    dplyr::select(cluster, has_depression_dysthymia)
  
  test_data <- with(test_data_raw, table(has_depression_dysthymia, cluster))
  
  print(fisher.test(test_data, simulate.p.value = TRUE))
}

test_ECG <- function(
    data,
    cluster_assignments_df,
    method = c("Tukey"),
    QC_filtering = c("with_exclusions", "without_exclusions"),
    ...
) {
  if (QC_filtering == "with_exclusions") {
    test_data <- data$ecg
  } else if (QC_filtering == "without_exclusions") {
    test_data <- data$ecg_wo_excl
  } else {
    stop("This QC_filtering type is not supported.")
  }
  
  test_data <- test_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    mutate(cluster = paste0("cl. ", cluster)) %>% 
    pivot_longer(
      cols = NNMin_mean:lfhf,
      names_to = "variable"
    )
  
  if (method == "conover") {
    stop("This method is not available.")
  } else if (method == "Tukey") {
    unique(test_data$variable) %>% 
      walk(function(current_var) {
        test_data_subset <- test_data %>% 
          filter(variable == current_var)
        model <- aov(value ~ cluster, data = test_data_subset)
        if (summary(model)[[1]]$`Pr(>F)`[1] < 0.05) {
          print(current_var)
          print(summary(model))
          print(TukeyHSD(model))
        }
      })
  } else {
    stop("This method is not available.")
  }
}

test_imaging_ancova <- function(
    data,
    cluster_assignments_df,
    variable_type = c("thickness_revised", "volumina"),
    model_type = c("without_age2"),
    ...
) {
  if (variable_type == "thickness_revised") {
    test_data <- data$imaging_thickness_revised %>% 
      dplyr::rename(value = thickness,
                    variable = region)
  } else if (variable_type == "volumina") {
    test_data <- data$imaging_volumina %>% 
      rename(value = volumina,
             variable = region)
  } else {
    stop("This variable_type is not supported.")
  }
  
  if (model_type == "without_age2") {
    model_formula <- as.formula("value ~ cluster + age + ICV + sex")
  } else {
    stop("This formula_type is not supported.")
  }
  
  sex_data <- data$imputed_data %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    dplyr::select(ID, sex)
  
  test_data <- test_data %>% 
    right_join(cluster_assignments_df, by = "ID") %>% 
    filter(!cluster %in% c("controls", "lifetime cases")) %>% 
    left_join(data$age, by = "ID") %>% 
    left_join(sex_data, by = "ID") %>% 
    filter(!is.na(variable)) %>% 
    mutate(cluster = as.factor(cluster))
  
  unique(test_data$variable) %>% 
    set_names() %>% 
    map_dfr(function(current_var) {
      test_data_subset <- test_data %>% 
        filter(variable == current_var)
      model <- aov(formula = model_formula,
                   data = test_data_subset,
                   contrasts = list(cluster = contr.sum, sex = contr.sum))
      res <- Anova(model, type = "III")
      broom::tidy(res) %>% 
        mutate(imaging_type = variable_type)
    }, .id = "variable")
}
