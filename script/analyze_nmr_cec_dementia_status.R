options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(logistf)
  library(pROC)
  library(purrr)
  library(readr)
  library(readxl)
  library(scales)
  library(stringr)
  library(tidyr)
})

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
if (is.na(script_path)) {
  script_dir <- normalizePath(file.path(getwd(), "script"), winslash = "/", mustWork = FALSE)
} else {
  script_dir <- dirname(script_path)
}
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results", "nmr_cec_dementia_status")
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

analysis_date <- "2026-09-03"
fdr_alpha <- 0.05
random_seed <- 20260903L
set.seed(random_seed)

input_paths <- list(
  nmr = file.path(project_root, "clean_data", "ASPREE_merged_lipoprofile_clean.csv"),
  nmr_codebook = file.path(project_root, "clean_data", "ASPREE_merged_lipoprofile_codebook.xlsx"),
  cec = file.path(project_root, "outputs", "cec_analysis", "cec_sample_summary.csv"),
  metadata = file.path(project_root, "raw_data", "subject_metadata.xlsx"),
  clinical = file.path(project_root, "raw_data", "clinical.csv"),
  medication = file.path(project_root, "raw_data", "med.csv"),
  allocation = file.path(project_root, "raw_data", "allocation.csv"),
  baseline_diagnoses = file.path(project_root, "raw_data", "BL_Dx.csv")
)
stopifnot(all(file.exists(unlist(input_paths))))

nmr_raw <- read_csv(input_paths$nmr, show_col_types = FALSE)
codebook <- read_excel(input_paths$nmr_codebook)
cec_raw <- read_csv(input_paths$cec, show_col_types = FALSE)
metadata_raw <- read_excel(input_paths$metadata)
clinical_raw <- read_csv(input_paths$clinical, show_col_types = FALSE)
med_raw <- read_csv(input_paths$medication, show_col_types = FALSE)
allocation_raw <- read_csv(input_paths$allocation, show_col_types = FALSE)
baseline_raw <- read_csv(input_paths$baseline_diagnoses, show_col_types = FALSE)

assert_unique <- function(df, key, label) {
  duplicate_keys <- df %>% count(.data[[key]], name = "n") %>% filter(n > 1)
  if (nrow(duplicate_keys) > 0) stop(label, " contains duplicate ", key, " values")
}

nmr <- nmr_raw %>%
  rename(
    subject_id = `Subject ID`,
    labcorp_accession = `Labcorp Accession Number`,
    collection_date = `Collection Date`,
    comment = Comment
  )

nmr_ids <- nmr$subject_id
cec <- cec_raw %>% filter(sample %in% nmr_ids) %>% rename(subject_id = sample)
metadata <- metadata_raw %>% filter(subject_id %in% nmr_ids)
clinical <- clinical_raw %>% filter(subject_id %in% nmr_ids)
med <- med_raw %>% filter(subject_id %in% nmr_ids)
allocation <- allocation_raw %>% filter(subject_id %in% nmr_ids)
baseline <- baseline_raw %>% filter(subject_id %in% nmr_ids)

assert_unique(nmr, "subject_id", "NMR data")
assert_unique(cec, "subject_id", "CEC summary")
assert_unique(metadata, "subject_id", "Subject metadata")
assert_unique(clinical, "subject_id", "Clinical data")
assert_unique(med, "subject_id", "Medication data")
assert_unique(allocation, "subject_id", "Allocation data")
assert_unique(baseline, "subject_id", "Baseline diagnosis data")

if (length(unique(nmr_ids)) != 174L) {
  stop("Expected 174 unique NMR subjects; found ", length(unique(nmr_ids)))
}

identifier_columns <- c("labcorp_accession", "subject_id", "collection_date", "comment")
all_nmr_features <- setdiff(names(nmr), identifier_columns)
constant_features <- c("LPX", "LPZ")
duplicate_aliases <- c("MVX_VAL", "MVX_LEU", "MVX_ILEU", "MVX_GLYCA")
ketone_features <- c("KetBod", "B-HB", "AcAc", "Acetone")
excluded_features <- c(constant_features, duplicate_aliases, ketone_features)
nmr_features <- setdiff(all_nmr_features, excluded_features)
cec_features <- c("cec_index_plate_qc", "cec_index_global_qc")

if (length(nmr_features) != 45L) {
  stop("Expected 45 included NMR features after exclusions; found ", length(nmr_features))
}

feature_lookup <- codebook %>%
  transmute(
    feature = abbreviated_variable_name,
    feature_label = variable_name,
    feature_group = variable_group_description,
    unit = measurement_units
  ) %>%
  distinct(feature, .keep_all = TRUE)

cec_lookup <- tibble(
  feature = cec_features,
  feature_label = c("CEC index normalized to plate QC", "CEC index normalized to global QC"),
  feature_group = "Cholesterol efflux capacity",
  unit = "Index"
)

feature_manifest <- tibble(feature = all_nmr_features) %>%
  left_join(feature_lookup, by = "feature") %>%
  mutate(
    modality = "NMR",
    analysis_status = case_when(
      feature %in% ketone_features ~ "Excluded: ketone marker affected by alcohol contamination",
      feature %in% constant_features ~ "Excluded: constant feature",
      feature %in% duplicate_aliases ~ "Excluded: exact duplicate alias",
      TRUE ~ "Included"
    ),
    duplicate_of = case_when(
      feature == "MVX_VAL" ~ "Val",
      feature == "MVX_LEU" ~ "Leu",
      feature == "MVX_ILEU" ~ "Ileu",
      feature == "MVX_GLYCA" ~ "GLYCA",
      TRUE ~ NA_character_
    )
  ) %>%
  bind_rows(
    cec_lookup %>%
      mutate(modality = "CEC", analysis_status = "Included", duplicate_of = NA_character_)
  ) %>%
  select(modality, feature, feature_label, feature_group, unit, analysis_status, duplicate_of)

analysis_df <- nmr %>%
  select(subject_id, all_of(nmr_features)) %>%
  inner_join(metadata, by = "subject_id") %>%
  inner_join(clinical %>% select(subject_id, CVD_protocol, Stroke), by = "subject_id") %>%
  inner_join(med %>% select(subject_id, AVx_HDL_Source, AVx_MedPres, AVx_Name_Med), by = "subject_id") %>%
  inner_join(allocation %>% select(subject_id, Treatment), by = "subject_id") %>%
  inner_join(baseline %>% select(subject_id, HTN_deriv, Diab_deriv, Dyslipidemia), by = "subject_id") %>%
  inner_join(cec %>% select(subject_id, all_of(cec_features)), by = "subject_id") %>%
  mutate(
    dementia_diagnosis = factor(dementia_diagnosis, levels = c("Normal", "Dementia")),
    dementia_bin = case_when(
      dementia_diagnosis == "Dementia" ~ 1L,
      dementia_diagnosis == "Normal" ~ 0L,
      TRUE ~ NA_integer_
    ),
    gender = relevel(factor(gender), ref = "Woman"),
    apoe_group = relevel(factor(apoe_group), ref = "apoe3e3"),
    Treatment = relevel(factor(Treatment), ref = "Placebo"),
    source_group = case_when(
      AVx_HDL_Source == "Yr 2" ~ "Yr 2",
      AVx_HDL_Source == "Yr 3" ~ "Yr 3",
      AVx_HDL_Source %in% c("Yr 5", "Yr 6", "Yr 7") ~ "Yr 5+",
      TRUE ~ NA_character_
    ),
    source_group = relevel(factor(source_group), ref = "Yr 3"),
    medication = if_else(!is.na(AVx_MedPres) & as.character(AVx_MedPres) == "1", 1L, 0L),
    CVD = as.integer(CVD_protocol),
    stroke = as.integer(Stroke),
    diabetes = as.integer(Diab_deriv),
    hypertension = as.integer(HTN_deriv),
    dyslipidemia = as.integer(Dyslipidemia)
  )

for (feature in c(nmr_features, cec_features)) {
  analysis_df[[feature]] <- suppressWarnings(as.numeric(analysis_df[[feature]]))
}

base_covariates <- c("gender", "apoe_group", "Treatment", "source_group")
expanded_covariates <- c(
  base_covariates, "medication", "CVD", "stroke", "diabetes", "hypertension", "dyslipidemia"
)
required_fields <- c("dementia_bin", base_covariates, expanded_covariates)
required_fields <- unique(required_fields)

if (nrow(analysis_df) != 174L) stop("Expected 174 linked subjects; found ", nrow(analysis_df))
if (anyNA(analysis_df[required_fields])) {
  missing_counts <- colSums(is.na(analysis_df[required_fields]))
  stop("Missing required model fields: ", paste(names(missing_counts)[missing_counts > 0], collapse = ", "))
}
if (!identical(sort(unique(analysis_df$dementia_bin)), c(0L, 1L))) {
  stop("Dementia outcome must contain exactly 0 and 1")
}
if (anyNA(analysis_df[c(nmr_features, cec_features)])) {
  stop("Included biomarker features unexpectedly contain missing values")
}
if (any(vapply(analysis_df[c(nmr_features, cec_features)], function(x) sd(x) <= 0, logical(1)))) {
  stop("At least one included biomarker has zero variance")
}

model_specifications <- list(base = base_covariates, expanded = expanded_covariates)

cohort_summary <- tibble(
  metric = c(
    "Linked analysis cohort", "Dementia cases", "Normal controls",
    "Included NMR features", "Excluded ketone features", "Included CEC features"
  ),
  value = c(
    as.character(nrow(analysis_df)),
    as.character(sum(analysis_df$dementia_bin == 1)),
    as.character(sum(analysis_df$dementia_bin == 0)),
    as.character(length(nmr_features)),
    paste(ketone_features, collapse = ", "),
    as.character(length(cec_features))
  )
)

categorical_variables <- c(
  gender = "Sex", apoe_group = "APOE genotype", Treatment = "Allocation",
  source_group = "AVx HDL source", medication = "Medication prescribed",
  CVD = "CVD", stroke = "Stroke", diabetes = "Diabetes",
  hypertension = "Hypertension", dyslipidemia = "Dyslipidemia"
)

covariate_by_dementia <- imap_dfr(categorical_variables, function(label, variable_name) {
  analysis_df %>%
    transmute(
      variable = label,
      level = as.character(.data[[variable_name]]),
      dementia_diagnosis = as.character(dementia_diagnosis)
    ) %>%
    count(variable, level, dementia_diagnosis, name = "n") %>%
    group_by(variable, dementia_diagnosis) %>%
    mutate(percent_within_status = 100 * n / sum(n)) %>%
    ungroup()
})

source_by_dementia <- analysis_df %>%
  transmute(
    source_group = factor(source_group, levels = c("Yr 2", "Yr 3", "Yr 5+")),
    dementia_diagnosis = factor(dementia_diagnosis, levels = c("Normal", "Dementia"))
  ) %>%
  count(source_group, dementia_diagnosis, .drop = FALSE, name = "n") %>%
  group_by(source_group) %>%
  mutate(percent_within_source = 100 * n / sum(n)) %>%
  ungroup()

separation_summary <- covariate_by_dementia %>%
  select(variable, level, dementia_diagnosis, n) %>%
  complete(
    nesting(variable, level),
    dementia_diagnosis = c("Normal", "Dementia"),
    fill = list(n = 0L)
  ) %>%
  pivot_wider(names_from = dementia_diagnosis, values_from = n, values_fill = 0) %>%
  mutate(zero_cell = Normal == 0 | Dementia == 0) %>%
  arrange(desc(zero_cell), variable, level)

medication_summary <- analysis_df %>%
  mutate(
    medication_name = if_else(
      is.na(AVx_Name_Med), "No statin-related medication recorded", AVx_Name_Med
    )
  ) %>%
  count(medication, medication_name, dementia_diagnosis, name = "n") %>%
  group_by(dementia_diagnosis) %>%
  mutate(percent_within_status = 100 * n / sum(n)) %>%
  ungroup()

feature_metadata <- feature_manifest %>%
  filter(analysis_status == "Included") %>%
  select(modality, feature, feature_label, feature_group, unit)

rank_biserial_results <- map_dfr(c(nmr_features, cec_features), function(feature) {
  x <- analysis_df[[feature]]
  y <- analysis_df$dementia_bin
  n_case <- sum(y == 1)
  n_control <- sum(y == 0)
  ranked <- rank(x, ties.method = "average")
  u_case <- sum(ranked[y == 1]) - n_case * (n_case + 1) / 2
  auc_case_higher <- u_case / (n_case * n_control)
  test <- suppressWarnings(wilcox.test(x[y == 1], x[y == 0], exact = FALSE))
  tibble(
    modality = if_else(feature %in% nmr_features, "NMR", "CEC"),
    feature = feature,
    n = length(x),
    cases = n_case,
    controls = n_control,
    median_dementia = median(x[y == 1]),
    median_normal = median(x[y == 0]),
    rank_biserial = 2 * auc_case_higher - 1,
    p_value = test$p.value
  )
}) %>%
  left_join(feature_metadata, by = c("modality", "feature")) %>%
  group_by(modality) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    direction = case_when(
      rank_biserial > 0 ~ "Higher in Dementia",
      rank_biserial < 0 ~ "Higher in Normal",
      TRUE ~ "No direction"
    )
  ) %>%
  ungroup() %>%
  arrange(modality, fdr, p_value)

empty_firth_row <- function(feature, modality, model_name, n, cases, controls, status) {
  tibble(
    modality = modality, feature = feature, model = model_name,
    n = n, cases = cases, controls = controls,
    odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    log_odds = NA_real_, p_value = NA_real_,
    full_iterations = NA_integer_, profile_iterations = NA_integer_,
    max_fit_convergence_metric = NA_real_, max_profile_convergence_metric = NA_real_,
    status = status
  )
}

fit_firth_feature <- function(feature, modality, model_name, covariates) {
  model_df <- analysis_df %>%
    select(dementia_bin, all_of(covariates), all_of(feature)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  model_df$z_feature <- as.numeric(scale(model_df[[feature]]))
  formula_text <- paste(
    "dementia_bin ~ z_feature +",
    paste(covariates, collapse = " + ")
  )
  fit_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      logistf(
        as.formula(formula_text), data = model_df,
        pl = TRUE, plconf = 2, firth = TRUE
      ),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(empty_firth_row(
      feature, modality, model_name, nrow(model_df), sum(model_df$dementia_bin == 1),
      sum(model_df$dementia_bin == 0), paste("ERROR:", conditionMessage(fit))
    ))
  }
  beta <- unname(fit$coefficients["z_feature"])
  lo <- unname(fit$ci.lower["z_feature"])
  hi <- unname(fit$ci.upper["z_feature"])
  p <- unname(fit$prob["z_feature"])
  finite_result <- all(is.finite(c(beta, lo, hi, p)))
  status <- if (!finite_result) {
    "ERROR: non-finite feature estimate"
  } else if (length(fit_warnings) > 0) {
    paste("WARNING:", paste(unique(fit_warnings), collapse = " | "))
  } else {
    "OK"
  }
  tibble(
    modality = modality,
    feature = feature,
    model = model_name,
    n = nrow(model_df),
    cases = sum(model_df$dementia_bin == 1),
    controls = sum(model_df$dementia_bin == 0),
    odds_ratio = exp(beta),
    conf_low = exp(lo),
    conf_high = exp(hi),
    log_odds = beta,
    p_value = p,
    full_iterations = unname(fit$iter["full"]),
    profile_iterations = sum(fit$pl.iter[2, ], na.rm = TRUE),
    max_fit_convergence_metric = max(abs(fit$conv), na.rm = TRUE),
    max_profile_convergence_metric = max(abs(fit$pl.conv), na.rm = TRUE),
    status = status
  )
}

firth_results <- map_dfr(names(model_specifications), function(model_name) {
  covariates <- model_specifications[[model_name]]
  map_dfr(c(nmr_features, cec_features), function(feature) {
    fit_firth_feature(
      feature = feature,
      modality = if_else(feature %in% nmr_features, "NMR", "CEC"),
      model_name = model_name,
      covariates = covariates
    )
  })
}) %>%
  left_join(feature_metadata, by = c("modality", "feature")) %>%
  group_by(modality, model) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    direction = case_when(
      odds_ratio > 1 ~ "Higher dementia odds",
      odds_ratio < 1 ~ "Lower dementia odds",
      TRUE ~ "No direction"
    ),
    significant_fdr = !is.na(fdr) & fdr < fdr_alpha
  ) %>%
  ungroup() %>%
  arrange(modality, model, fdr, p_value)

if (any(str_detect(firth_results$status, "^ERROR"))) {
  stop("At least one prespecified Firth model failed; inspect results before proceeding")
}

metabolic_vulnerability_map <- tibble(
  analysis_feature = c(
    "MVX_MVX", "MVX_IVX", "MVX_MMX", "GLYCA", "MVX_MSCHDLP",
    "MVX_CTR", "Val", "Leu", "Ileu"
  ),
  mvx_variable = c(
    "MVX_MVX", "MVX_IVX", "MVX_MMX", "MVX_GLYCA", "MVX_MSCHDLP",
    "MVX_CTR", "MVX_VAL", "MVX_LEU", "MVX_ILEU"
  ),
  display_label = c(
    "Metabolic vulnerability index", "Inflammation vulnerability index",
    "Metabolic malnutrition index", "GlycA", "Medium/small HDL particles",
    "Citrate", "Valine", "Leucine", "Isoleucine"
  ),
  display_order = seq_len(9)
)

metabolic_vulnerability_results <- firth_results %>%
  filter(modality == "NMR", feature %in% metabolic_vulnerability_map$analysis_feature) %>%
  left_join(metabolic_vulnerability_map, by = c("feature" = "analysis_feature")) %>%
  arrange(display_order, model) %>%
  select(
    mvx_variable, analysis_feature = feature, display_label, display_order,
    model, n, cases, controls, odds_ratio, conf_low, conf_high, log_odds,
    p_value, fdr, direction, significant_fdr, status
  )
stopifnot(
  nrow(metabolic_vulnerability_results) == 18L,
  n_distinct(metabolic_vulnerability_results$mvx_variable) == 9L,
  !any(metabolic_vulnerability_results$analysis_feature %in% ketone_features)
)

fit_hdl_profile <- function(model_name, covariates) {
  hdl_features <- paste0("H", 1:7, "P")
  model_df <- analysis_df %>%
    select(dementia_bin, all_of(covariates), all_of(hdl_features)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  for (feature in hdl_features) {
    model_df[[paste0("z_", feature)]] <- as.numeric(scale(model_df[[feature]]))
  }
  z_features <- paste0("z_", hdl_features)
  full_formula <- as.formula(paste(
    "dementia_bin ~",
    paste(c(z_features, covariates), collapse = " + ")
  ))
  reduced_formula <- as.formula(paste("dementia_bin ~", paste(covariates, collapse = " + ")))
  fit_warnings <- character()
  full_fit <- withCallingHandlers(
    logistf(full_formula, data = model_df, pl = TRUE, plconf = 2:8, firth = TRUE),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  reduced_fit <- logistf(reduced_formula, data = model_df, pl = FALSE, firth = TRUE)
  comparison <- anova(full_fit, reduced_fit, method = "nested")
  global <- tibble(
    model = model_name,
    n = nrow(model_df),
    cases = sum(model_df$dementia_bin == 1),
    controls = sum(model_df$dementia_bin == 0),
    degrees_freedom = comparison$df,
    penalized_likelihood_ratio_chisq = unname(comparison$chisq),
    p_value = unname(comparison$pval),
    max_fit_convergence_metric = max(abs(full_fit$conv), na.rm = TRUE),
    max_profile_convergence_metric = max(abs(full_fit$pl.conv), na.rm = TRUE),
    status = if_else(
      length(fit_warnings) == 0, "OK",
      paste("WARNING:", paste(unique(fit_warnings), collapse = " | "))
    )
  )
  coefficients <- map_dfr(seq_along(hdl_features), function(i) {
    feature <- hdl_features[[i]]
    term <- z_features[[i]]
    tibble(
      model = model_name,
      feature = feature,
      odds_ratio = exp(unname(full_fit$coefficients[term])),
      conf_low = exp(unname(full_fit$ci.lower[term])),
      conf_high = exp(unname(full_fit$ci.upper[term])),
      log_odds = unname(full_fit$coefficients[term]),
      p_value = unname(full_fit$prob[term])
    )
  })
  list(global = global, coefficients = coefficients)
}

hdl_profile_fits <- imap(model_specifications, ~ fit_hdl_profile(.y, .x))
hdl_profile_global <- bind_rows(map(hdl_profile_fits, "global"))
hdl_profile_coefficients <- bind_rows(map(hdl_profile_fits, "coefficients")) %>%
  left_join(feature_metadata %>% select(feature, feature_label, unit), by = "feature")
hdl_profile_correlation <- cor(analysis_df[paste0("H", 1:7, "P")], method = "spearman") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature")

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold_id <- integer(length(y))
  for (class_value in sort(unique(y))) {
    idx <- sample(which(y == class_value))
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold_id
}

safe_biomarker_names <- make.names(c(nmr_features, cec_features), unique = TRUE)
biomarker_name_map <- tibble(
  feature = c(nmr_features, cec_features),
  matrix_feature = paste0("bio__", safe_biomarker_names),
  modality = c(rep("NMR", length(nmr_features)), rep("CEC", length(cec_features)))
)
biomarker_matrix <- as.matrix(analysis_df[c(nmr_features, cec_features)])
storage.mode(biomarker_matrix) <- "double"
colnames(biomarker_matrix) <- biomarker_name_map$matrix_feature
covariate_matrix <- model.matrix(
  ~ gender + apoe_group + Treatment + source_group + medication + CVD + stroke +
    diabetes + hypertension + dyslipidemia - 1,
  data = analysis_df
)
full_prediction_matrix <- cbind(covariate_matrix, biomarker_matrix)
prediction_y <- analysis_df$dementia_bin

fit_cv_glmnet <- function(x_train, y_train, inner_fold_id) {
  keep_columns <- apply(x_train, 2, function(x) is.finite(sd(x)) && sd(x) > 0)
  x_train_kept <- x_train[, keep_columns, drop = FALSE]
  fit <- cv.glmnet(
    x = x_train_kept,
    y = y_train,
    family = "binomial",
    alpha = 0.5,
    foldid = inner_fold_id,
    type.measure = "auc",
    standardize = TRUE,
    intercept = TRUE,
    keep = FALSE
  )
  list(fit = fit, keep_columns = keep_columns)
}

predict_cv_glmnet <- function(fit_object, newx, lambda_rule) {
  kept_newx <- newx[, fit_object$keep_columns, drop = FALSE]
  as.numeric(predict(fit_object$fit, newx = kept_newx, s = lambda_rule, type = "response"))
}

extract_selected_biomarkers <- function(fit_object, repeat_id, outer_fold, lambda_rule) {
  lambda_value <- if (lambda_rule == "lambda.min") {
    fit_object$fit$lambda.min
  } else {
    fit_object$fit$lambda.1se
  }
  coefficients <- as.matrix(coef(fit_object$fit, s = lambda_rule))[, 1]
  tibble(matrix_feature = names(coefficients), coefficient = as.numeric(coefficients)) %>%
    filter(matrix_feature %in% biomarker_name_map$matrix_feature, abs(coefficient) > 1e-12) %>%
    mutate(repeat_id = repeat_id, outer_fold = outer_fold, rule = lambda_rule, lambda = lambda_value)
}

outer_repeats <- 10L
outer_folds <- 5L
inner_folds <- 5L
prediction_rows <- list()
selection_rows <- list()
prediction_counter <- 1L
selection_counter <- 1L

for (repeat_id in seq_len(outer_repeats)) {
  outer_fold_id <- make_stratified_folds(prediction_y, outer_folds, random_seed + repeat_id)
  for (outer_fold in seq_len(outer_folds)) {
    train_index <- outer_fold_id != outer_fold
    test_index <- !train_index
    inner_fold_id <- make_stratified_folds(
      prediction_y[train_index], inner_folds,
      random_seed + repeat_id * 100L + outer_fold
    )
    covariate_fit <- fit_cv_glmnet(
      covariate_matrix[train_index, , drop = FALSE], prediction_y[train_index], inner_fold_id
    )
    full_fit <- fit_cv_glmnet(
      full_prediction_matrix[train_index, , drop = FALSE], prediction_y[train_index], inner_fold_id
    )
    for (lambda_rule in c("lambda.min", "lambda.1se")) {
      prediction_rows[[prediction_counter]] <- tibble(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        subject_id = analysis_df$subject_id[test_index],
        dementia_bin = prediction_y[test_index],
        model = "Covariates only",
        rule = lambda_rule,
        predicted_probability = predict_cv_glmnet(
          covariate_fit, covariate_matrix[test_index, , drop = FALSE], lambda_rule
        )
      )
      prediction_counter <- prediction_counter + 1L
      prediction_rows[[prediction_counter]] <- tibble(
        repeat_id = repeat_id,
        outer_fold = outer_fold,
        subject_id = analysis_df$subject_id[test_index],
        dementia_bin = prediction_y[test_index],
        model = "Covariates + biomarkers",
        rule = lambda_rule,
        predicted_probability = predict_cv_glmnet(
          full_fit, full_prediction_matrix[test_index, , drop = FALSE], lambda_rule
        )
      )
      prediction_counter <- prediction_counter + 1L
      selected <- extract_selected_biomarkers(full_fit, repeat_id, outer_fold, lambda_rule)
      if (nrow(selected) > 0) {
        selection_rows[[selection_counter]] <- selected
        selection_counter <- selection_counter + 1L
      }
    }
  }
}

elastic_net_outer_predictions <- bind_rows(prediction_rows)
elastic_net_repeat_performance <- elastic_net_outer_predictions %>%
  group_by(repeat_id, model, rule) %>%
  summarise(
    n = n(),
    cases = sum(dementia_bin == 1),
    controls = sum(dementia_bin == 0),
    auc = as.numeric(pROC::auc(
      response = dementia_bin, predictor = predicted_probability,
      levels = c(0, 1), direction = "<", quiet = TRUE
    )),
    .groups = "drop"
  )

elastic_net_summary <- elastic_net_repeat_performance %>%
  group_by(model, rule) %>%
  summarise(
    repeats = n(),
    mean_auc = mean(auc),
    sd_auc = sd(auc),
    median_auc = median(auc),
    min_auc = min(auc),
    max_auc = max(auc),
    .groups = "drop"
  )

elastic_net_delta_auc <- elastic_net_repeat_performance %>%
  select(repeat_id, rule, model, auc) %>%
  pivot_wider(names_from = model, values_from = auc) %>%
  mutate(delta_auc = `Covariates + biomarkers` - `Covariates only`) %>%
  arrange(rule, repeat_id)

selected_raw <- if (length(selection_rows) == 0) {
  tibble(
    matrix_feature = character(), coefficient = double(), repeat_id = integer(),
    outer_fold = integer(), rule = character(), lambda = double()
  )
} else {
  bind_rows(selection_rows)
}
selected_counts <- selected_raw %>%
  group_by(rule, matrix_feature) %>%
  summarise(
    selection_count = n(),
    mean_coefficient_when_selected = mean(coefficient),
    .groups = "drop"
  )

elastic_net_selection_stability <- expand_grid(
  rule = c("lambda.min", "lambda.1se"),
  matrix_feature = biomarker_name_map$matrix_feature
) %>%
  left_join(selected_counts, by = c("rule", "matrix_feature")) %>%
  mutate(
    selection_count = replace_na(selection_count, 0L),
    selection_frequency = selection_count / (outer_repeats * outer_folds)
  ) %>%
  left_join(biomarker_name_map, by = "matrix_feature") %>%
  left_join(feature_metadata, by = c("modality", "feature")) %>%
  arrange(rule, desc(selection_frequency), feature)

prediction_integrity <- elastic_net_outer_predictions %>%
  count(repeat_id, model, rule, subject_id, name = "prediction_count")
stopifnot(
  nrow(firth_results) == 2L * (length(nmr_features) + length(cec_features)),
  all(firth_results$n == 174L),
  all(is.finite(firth_results$odds_ratio)),
  all(is.finite(firth_results$conf_low)),
  all(is.finite(firth_results$conf_high)),
  all(is.finite(firth_results$p_value)),
  !any(firth_results$feature %in% excluded_features),
  nrow(hdl_profile_global) == 2L,
  nrow(hdl_profile_coefficients) == 14L,
  nrow(elastic_net_outer_predictions) == outer_repeats * nrow(analysis_df) * 4L,
  nrow(prediction_integrity) == nrow(elastic_net_outer_predictions),
  all(prediction_integrity$prediction_count == 1L),
  all(is.finite(elastic_net_outer_predictions$predicted_probability)),
  all(elastic_net_outer_predictions$predicted_probability >= 0),
  all(elastic_net_outer_predictions$predicted_probability <= 1),
  nrow(elastic_net_repeat_performance) == outer_repeats * 4L,
  nrow(elastic_net_selection_stability) == 2L * (length(nmr_features) + length(cec_features))
)

figure_palette <- c(
  ink = "#252525", grey_dark = "#5B5B5B", grey_mid = "#A6A6A6",
  grey_light = "#D9D9D9", blue = "#2F75B5", blue_light = "#C9DDEC",
  red = "#B33A3A", normal = "#9C9C9C", dementia = "#2F75B5"
)
width_mm = 183

theme_nature <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      line = element_line(linewidth = 0.35, colour = figure_palette[["ink"]]),
      axis.line = element_line(linewidth = 0.35, colour = figure_palette[["ink"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = figure_palette[["ink"]]),
      axis.ticks.length = grid::unit(1.8, "pt"),
      axis.title = element_text(size = base_size, colour = figure_palette[["ink"]]),
      axis.text = element_text(size = base_size - 0.5, colour = figure_palette[["ink"]]),
      legend.position = "top",
      legend.justification = "left",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.6),
      legend.key.width = grid::unit(10, "pt"),
      legend.key.height = grid::unit(7, "pt"),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_text(size = base_size - 0.7, colour = figure_palette[["grey_dark"]], hjust = 0),
      plot.margin = margin(5.5, 6.5, 5.5, 6.5, unit = "pt"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

save_plot <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svg_path <- file.path(plots_dir, paste0(stem, ".svg"))
  pdf_path <- file.path(plots_dir, paste0(stem, ".pdf"))
  tiff_path <- file.path(plots_dir, paste0(stem, ".tiff"))
  png_path <- file.path(plots_dir, paste0(stem, ".png"))

  svglite::svglite(svg_path, width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()

  pdf_device <- tolower(Sys.getenv("NATURE_FIGURE_PDF_DEVICE", "postscript"))
  if (identical(pdf_device, "cairo")) {
    grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Helvetica", bg = "white")
    print(plot)
    grDevices::dev.off()
  } else {
    gs_bin <- Sys.which("gs")
    if (!nzchar(gs_bin)) stop("Ghostscript is required for editable PDF export when Cairo is unavailable")
    eps_path <- tempfile(fileext = ".eps")
    grDevices::postscript(
      eps_path, onefile = FALSE, horizontal = FALSE, paper = "special",
      width = width_in, height = height_in, family = "Helvetica", bg = "white"
    )
    print(plot)
    grDevices::dev.off()
    gs_status <- suppressWarnings(system2(
      gs_bin,
      args = c(
        "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dQUIET", "-dEPSCrop",
        "-dAutoRotatePages=/None", "-sDEVICE=pdfwrite",
        "-dCompatibilityLevel=1.4", paste0("-sOutputFile=", shQuote(pdf_path)), shQuote(eps_path)
      ),
      stdout = FALSE, stderr = FALSE
    ))
    unlink(eps_path)
    if (!identical(gs_status, 0L) || !file.exists(pdf_path)) stop("PDF export failed for ", stem)
  }

  ragg::agg_tiff(tiff_path, width = width_in, height = height_in, units = "in", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in", res = 300, background = "white")
  print(plot)
  grDevices::dev.off()
}

log_limits <- function(low, high, margin = 0.08) {
  bounds <- range(c(low, high), na.rm = TRUE)
  stopifnot(all(is.finite(bounds)), all(bounds > 0))
  log_bounds <- log(bounds)
  span <- diff(log_bounds)
  if (span == 0) span <- 0.2
  exp(c(log_bounds[1] - margin * span, log_bounds[2] + margin * span))
}

outcome_counts <- analysis_df %>%
  count(dementia_diagnosis, name = "n") %>%
  mutate(dementia_diagnosis = factor(dementia_diagnosis, levels = c("Normal", "Dementia")))
outcome_plot <- ggplot(outcome_counts, aes(x = dementia_diagnosis, y = n, fill = dementia_diagnosis)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.45, size = 2.3, family = "Helvetica", colour = figure_palette[["ink"]]) +
  scale_fill_manual(values = c("Normal" = figure_palette[["normal"]], "Dementia" = figure_palette[["dementia"]]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), breaks = pretty_breaks(4)) +
  labs(x = NULL, y = "Participants") +
  theme_nature()
save_plot(outcome_plot, "dementia_status_counts", 89, 60)

source_plot <- ggplot(
  source_by_dementia,
  aes(x = source_group, y = n, fill = dementia_diagnosis)
) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, colour = "white", linewidth = 0.2) +
  geom_text(
    aes(label = n), position = position_dodge(width = 0.72),
    vjust = -0.35, size = 2.1, family = "Helvetica"
  ) +
  scale_fill_manual(values = c("Normal" = figure_palette[["normal"]], "Dementia" = figure_palette[["dementia"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), breaks = pretty_breaks(5)) +
  labs(x = "Sample-source group", y = "Participants") +
  theme_nature()
save_plot(source_plot, "sample_source_by_dementia", 89, 67)

rank_plot_data <- rank_biserial_results %>%
  filter(modality == "NMR") %>%
  slice_min(order_by = fdr, n = 15, with_ties = FALSE) %>%
  arrange(rank_biserial) %>%
  mutate(
    feature = factor(feature, levels = feature),
    evidence = if_else(fdr < fdr_alpha, "FDR < 0.05", "Other")
  )
rank_limit <- max(abs(rank_plot_data$rank_biserial), na.rm = TRUE) * 1.18
rank_plot <- ggplot(rank_plot_data, aes(x = rank_biserial, y = feature)) +
  geom_vline(xintercept = 0, colour = figure_palette[["grey_mid"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rank_biserial, yend = feature), colour = figure_palette[["grey_mid"]], linewidth = 0.4) +
  geom_point(aes(fill = evidence), shape = 21, colour = "white", stroke = 0.25, size = 2.0) +
  scale_fill_manual(values = c("FDR < 0.05" = figure_palette[["red"]], "Other" = figure_palette[["grey_dark"]]), guide = "none") +
  scale_x_continuous(limits = c(-rank_limit, rank_limit), breaks = pretty_breaks(5)) +
  labs(x = "Rank-biserial effect (positive = higher in Dementia)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.4))
save_plot(rank_plot, "unadjusted_nmr_rank_biserial", 89, 95)

nmr_feature_rank <- firth_results %>%
  filter(modality == "NMR") %>%
  group_by(feature) %>%
  summarise(rank_fdr = min(fdr), .groups = "drop") %>%
  slice_min(order_by = rank_fdr, n = 12, with_ties = FALSE) %>%
  arrange(rank_fdr)
nmr_forest_data <- firth_results %>%
  filter(modality == "NMR", feature %in% nmr_feature_rank$feature) %>%
  mutate(
    model = recode(model, base = "Base", expanded = "Expanded"),
    feature = factor(feature, levels = rev(nmr_feature_rank$feature))
  )
nmr_xlim <- log_limits(nmr_forest_data$conf_low, nmr_forest_data$conf_high, margin = 0.10)
nmr_forest_plot <- ggplot(nmr_forest_data, aes(x = odds_ratio, y = feature, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.52), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.52), size = 1.9, stroke = 0.5) +
  geom_point(
    data = nmr_forest_data %>% filter(significant_fdr),
    position = position_dodge(width = 0.52), shape = 1, size = 2.9,
    stroke = 0.65, colour = figure_palette[["red"]], show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(limits = nmr_xlim, breaks = log_breaks(n = 5), labels = label_number(accuracy = 0.1)) +
  labs(x = "Firth odds ratio per 1-SD higher biomarker (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.3), legend.position = "top")
save_plot(nmr_forest_plot, "nmr_firth_top_forest", 183, 112)

nmr_stability_data <- firth_results %>%
  filter(modality == "NMR") %>%
  select(feature, model, log_odds, p_value, fdr) %>%
  pivot_wider(names_from = model, values_from = c(log_odds, p_value, fdr)) %>%
  mutate(
    base_log2_or = log_odds_base / log(2),
    expanded_log2_or = log_odds_expanded / log(2),
    evidence = case_when(
      pmin(fdr_base, fdr_expanded, na.rm = TRUE) < fdr_alpha ~ "q < 0.05",
      pmin(p_value_base, p_value_expanded, na.rm = TRUE) < 0.05 ~ "P < 0.05 only",
      TRUE ~ "Other"
    )
  )
stability_limit <- max(abs(c(nmr_stability_data$base_log2_or, nmr_stability_data$expanded_log2_or))) * 1.18
nmr_stability_plot <- ggplot(nmr_stability_data, aes(x = base_log2_or, y = expanded_log2_or)) +
  geom_abline(slope = 1, intercept = 0, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_hline(yintercept = 0, colour = figure_palette[["grey_light"]], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = figure_palette[["grey_light"]], linewidth = 0.3) +
  geom_point(aes(fill = evidence), shape = 21, colour = "white", stroke = 0.25, size = 2.05) +
  scale_fill_manual(values = c(
    "q < 0.05" = figure_palette[["red"]],
    "P < 0.05 only" = figure_palette[["grey_dark"]],
    "Other" = figure_palette[["grey_light"]]
  )) +
  scale_x_continuous(limits = c(-stability_limit, stability_limit), breaks = pretty_breaks(4)) +
  scale_y_continuous(limits = c(-stability_limit, stability_limit), breaks = pretty_breaks(4)) +
  coord_fixed() +
  labs(x = "Base-model log2(OR per SD)", y = "Expanded-model log2(OR per SD)") +
  theme_nature() +
  theme(legend.position = "top")
save_plot(nmr_stability_plot, "nmr_firth_adjustment_stability", 89, 85)

mvx_display_order <- metabolic_vulnerability_map %>% arrange(display_order) %>% pull(display_label)
mvx_forest_data <- metabolic_vulnerability_results %>%
  mutate(
    model = recode(model, base = "Base", expanded = "Expanded"),
    display_label = factor(display_label, levels = rev(mvx_display_order))
  )
mvx_xlim <- log_limits(mvx_forest_data$conf_low, mvx_forest_data$conf_high, margin = 0.10)
mvx_forest_plot <- ggplot(mvx_forest_data, aes(x = odds_ratio, y = display_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.50), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.50), size = 1.9, stroke = 0.5) +
  geom_point(
    data = mvx_forest_data %>% filter(significant_fdr),
    position = position_dodge(width = 0.50), shape = 1, size = 2.9,
    stroke = 0.65, colour = figure_palette[["red"]], show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(limits = mvx_xlim, breaks = log_breaks(n = 5), labels = label_number(accuracy = 0.1)) +
  labs(x = "Firth odds ratio per 1-SD higher biomarker (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.3), legend.position = "top")
save_plot(mvx_forest_plot, "metabolic_vulnerability_firth_forest", 183, 100)

cec_forest_data <- firth_results %>%
  filter(modality == "CEC") %>%
  mutate(
    model = recode(model, base = "Base", expanded = "Expanded"),
    feature_label = recode(
      feature,
      cec_index_plate_qc = "Plate-QC index",
      cec_index_global_qc = "Global-QC index"
    ),
    feature_label = factor(feature_label, levels = c("Plate-QC index", "Global-QC index"))
  )
cec_xlim <- log_limits(cec_forest_data$conf_low, cec_forest_data$conf_high, margin = 0.10)
cec_forest_plot <- ggplot(cec_forest_data, aes(x = odds_ratio, y = feature_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.42), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.42), size = 1.9, stroke = 0.5) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(limits = cec_xlim, breaks = log_breaks(n = 4), labels = label_number(accuracy = 0.1)) +
  labs(x = "Firth odds ratio per 1-SD higher CEC (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.4), legend.position = "top")
save_plot(cec_forest_plot, "cec_firth_forest", 89, 60)

hdl_forest_data <- hdl_profile_coefficients %>%
  mutate(
    model = recode(model, base = "Base", expanded = "Expanded"),
    feature = factor(feature, levels = rev(paste0("H", 1:7, "P")))
  )
hdl_xlim <- log_limits(hdl_forest_data$conf_low, hdl_forest_data$conf_high, margin = 0.10)
hdl_forest_plot <- ggplot(hdl_forest_data, aes(x = odds_ratio, y = feature, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.46), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.46), size = 1.9, stroke = 0.5) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(limits = hdl_xlim, breaks = log_breaks(n = 5), labels = label_number(accuracy = 0.1)) +
  labs(x = "Joint-profile Firth odds ratio per 1-SD higher subclass (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.4), legend.position = "top")
save_plot(hdl_forest_plot, "hdl_profile_joint_firth_forest", 89, 82)

elastic_auc_plot_data <- elastic_net_repeat_performance %>%
  filter(rule == "lambda.1se") %>%
  mutate(model = factor(model, levels = c("Covariates only", "Covariates + biomarkers")))
elastic_auc_lower <- max(0, min(elastic_auc_plot_data$auc, na.rm = TRUE) - 0.03)
elastic_auc_upper <- min(1, max(elastic_auc_plot_data$auc, na.rm = TRUE) + 0.03)
elastic_auc_plot <- ggplot(elastic_auc_plot_data, aes(x = model, y = auc, group = repeat_id)) +
  geom_line(colour = figure_palette[["grey_mid"]], linewidth = 0.35) +
  geom_point(aes(fill = model), shape = 21, colour = "white", stroke = 0.2, size = 1.55) +
  stat_summary(
    aes(group = model), fun = mean, geom = "point", shape = 23,
    size = 2.7, stroke = 0.6, colour = figure_palette[["ink"]], fill = "white"
  ) +
  scale_fill_manual(values = c(
    "Covariates only" = figure_palette[["grey_dark"]],
    "Covariates + biomarkers" = figure_palette[["blue"]]
  ), guide = "none") +
  scale_y_continuous(
    limits = c(elastic_auc_lower, elastic_auc_upper),
    breaks = pretty_breaks(5), labels = label_number(accuracy = 0.01)
  ) +
  labs(x = NULL, y = "Nested-CV AUC") +
  theme_nature() +
  theme(axis.text.x = element_text(size = 6.2))
save_plot(elastic_auc_plot, "elastic_net_nested_auc", 89, 70)

selection_plot_data <- elastic_net_selection_stability %>%
  filter(rule == "lambda.1se") %>%
  slice_max(order_by = selection_frequency, n = 15, with_ties = FALSE) %>%
  arrange(selection_frequency) %>%
  mutate(feature = factor(feature, levels = feature))
if (max(selection_plot_data$selection_frequency) > 0) {
  selection_plot <- ggplot(selection_plot_data, aes(x = selection_frequency, y = feature)) +
    geom_segment(aes(x = 0, xend = selection_frequency, yend = feature), colour = figure_palette[["grey_mid"]], linewidth = 0.4) +
    geom_point(fill = figure_palette[["blue"]], colour = "white", shape = 21, stroke = 0.25, size = 2.0) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = label_percent(accuracy = 1)) +
    labs(x = "Selection frequency across outer training sets", y = NULL) +
    theme_nature() +
    theme(axis.text.y = element_text(size = 6.4))
} else {
  selection_plot <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "No biomarker was selected at lambda.1se",
      size = 2.4, family = "Helvetica", colour = figure_palette[["ink"]]
    ) +
    xlim(0, 1) + ylim(0, 1) + theme_void(base_family = "Helvetica")
}
save_plot(selection_plot, "elastic_net_selection_stability", 89, 90)

write_csv(analysis_df, file.path(results_dir, "analysis_dataset.csv"))
write_csv(feature_manifest, file.path(results_dir, "feature_manifest.csv"))
write_csv(cohort_summary, file.path(results_dir, "cohort_summary.csv"))
write_csv(covariate_by_dementia, file.path(results_dir, "covariate_by_dementia_summary.csv"))
write_csv(source_by_dementia, file.path(results_dir, "source_by_dementia_summary.csv"))
write_csv(separation_summary, file.path(results_dir, "separation_summary.csv"))
write_csv(medication_summary, file.path(results_dir, "medication_distribution.csv"))
write_csv(rank_biserial_results, file.path(results_dir, "unadjusted_rank_results.csv"))
write_csv(firth_results, file.path(results_dir, "firth_single_feature_results.csv"))
write_csv(metabolic_vulnerability_results, file.path(results_dir, "metabolic_vulnerability_results.csv"))
write_csv(hdl_profile_global, file.path(results_dir, "hdl_profile_global_tests.csv"))
write_csv(hdl_profile_coefficients, file.path(results_dir, "hdl_profile_joint_coefficients.csv"))
write_csv(hdl_profile_correlation, file.path(results_dir, "hdl_profile_spearman_correlation.csv"))
write_csv(elastic_net_outer_predictions, file.path(results_dir, "elastic_net_outer_predictions.csv"))
write_csv(elastic_net_repeat_performance, file.path(results_dir, "elastic_net_repeat_performance.csv"))
write_csv(elastic_net_summary, file.path(results_dir, "elastic_net_summary.csv"))
write_csv(elastic_net_delta_auc, file.path(results_dir, "elastic_net_delta_auc.csv"))
write_csv(elastic_net_selection_stability, file.path(results_dir, "elastic_net_selection_stability.csv"))

analysis_assumptions <- tibble(
  item = c(
    "Dementia endpoint", "Time-to-event data", "Estimator", "Age",
    "Sampling-source adjustment", "Medication", "Ketone panel", "Multiplicity",
    "Prediction interpretation", "Causal interpretation"
  ),
  specification = c(
    "Dementia is coded 1 and Normal is coded 0 for all 174 participants.",
    "Dementia diagnosis time and censoring time are unavailable; no survival model is fitted.",
    "Firth penalized-likelihood logistic regression is used because prespecified covariates contain complete or quasi-complete separation.",
    "Participant age and sample age are not included.",
    "AVx_HDL_Source is included in both models as Yr 2, Yr 3 or Yr 5+.",
    "AVx_MedPres equal to 1 is medication=1; missing is medication=0. Drug names are descriptive only.",
    "KetBod, B-HB, AcAc and Acetone are excluded because of alcohol contamination.",
    "BH FDR is controlled separately for NMR and CEC within each base/expanded model.",
    "Nested-CV AUC is exploratory and is not an absolute-risk or population-calibration estimate.",
    "All estimates are associations with recorded dementia status, not causal effects."
  )
)
write_csv(analysis_assumptions, file.path(results_dir, "analysis_assumptions.csv"))

table_files <- c(
  "analysis_dataset.csv", "analysis_assumptions.csv", "feature_manifest.csv", "cohort_summary.csv",
  "covariate_by_dementia_summary.csv", "source_by_dementia_summary.csv", "separation_summary.csv",
  "medication_distribution.csv", "unadjusted_rank_results.csv", "firth_single_feature_results.csv",
  "metabolic_vulnerability_results.csv", "hdl_profile_global_tests.csv",
  "hdl_profile_joint_coefficients.csv", "hdl_profile_spearman_correlation.csv",
  "elastic_net_outer_predictions.csv", "elastic_net_repeat_performance.csv",
  "elastic_net_summary.csv", "elastic_net_delta_auc.csv", "elastic_net_selection_stability.csv"
)
all_plot_files <- list.files(plots_dir, full.names = FALSE)
plot_files <- all_plot_files[tolower(tools::file_ext(all_plot_files)) %in% c("png", "pdf", "svg", "tiff")]
output_manifest <- bind_rows(
  tibble(file = table_files, type = "table", description = "Authoritative numerical analysis output"),
  tibble(file = file.path("plots", plot_files), type = "figure", description = "Report figure"),
  tibble(
    file = c(
      "nmr_cec_dementia_status_report.html", "nmr_cec_dementia_status_report.Rmd",
      "session_info.txt", "figure_contract.md", "figure_qa_notes.md"
    ),
    type = c("report", "source", "reproducibility", "figure protocol", "figure QA"),
    description = c(
      "Self-contained English HTML report", "Reproducible report source", "R session information",
      "Publication figure claim and export contract", "Publication figure validation record"
    )
  )
) %>%
  mutate(generated_on = analysis_date)
write_csv(output_manifest, file.path(results_dir, "output_manifest.csv"))

capture.output(sessionInfo(), file = file.path(results_dir, "session_info.txt"))

report_input <- file.path(results_dir, "nmr_cec_dementia_status_report.Rmd")
if (!file.exists(report_input)) stop("Report template not found: ", report_input)
rmarkdown::render(
  input = report_input,
  output_file = "nmr_cec_dementia_status_report.html",
  output_dir = results_dir,
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

cat("Analysis complete. Report:", file.path(results_dir, "nmr_cec_dementia_status_report.html"), "\n")
