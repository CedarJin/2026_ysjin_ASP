options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(purrr)
  library(readr)
  library(readxl)
  library(scales)
  library(stringr)
  library(survival)
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
results_dir <- file.path(project_root, "results", "nmr_cec_death_dsr")
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

analysis_date <- "2026-09-02"
fdr_alpha <- 0.05
set.seed(20260902)

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
  if (nrow(duplicate_keys) > 0) {
    stop(label, " contains duplicate ", key, " values")
  }
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
      mutate(
        modality = "CEC",
        analysis_status = "Included",
        duplicate_of = NA_character_
      )
  ) %>%
  select(modality, feature, feature_label, feature_group, unit, analysis_status, duplicate_of)

analysis_df <- nmr %>%
  select(subject_id, all_of(nmr_features)) %>%
  inner_join(metadata, by = "subject_id") %>%
  inner_join(clinical %>% select(subject_id, Death_DSR, CVD_protocol, Stroke), by = "subject_id") %>%
  inner_join(med %>% select(subject_id, AVx_HDL_Source, AVx_MedPres, AVx_Name_Med), by = "subject_id") %>%
  inner_join(allocation %>% select(subject_id, Treatment), by = "subject_id") %>%
  inner_join(baseline %>% select(subject_id, HTN_deriv, Diab_deriv, Dyslipidemia), by = "subject_id") %>%
  inner_join(cec %>% select(subject_id, all_of(cec_features)), by = "subject_id") %>%
  mutate(
    death_event = 1L,
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

required_fields <- c(
  "Death_DSR", "death_event", "gender", "apoe_group", "Treatment", "source_group",
  "medication", "CVD", "stroke", "diabetes", "hypertension", "dyslipidemia"
)
if (anyNA(analysis_df[required_fields])) {
  missing_counts <- colSums(is.na(analysis_df[required_fields]))
  stop("Missing required model fields: ", paste(names(missing_counts)[missing_counts > 0], collapse = ", "))
}
if (nrow(analysis_df) != 174L) {
  stop("Expected 174 complete linked subjects; found ", nrow(analysis_df))
}
if (any(analysis_df$Death_DSR <= 0)) {
  stop("Death_DSR must be positive")
}
if (any(vapply(analysis_df[c(nmr_features, cec_features)], function(x) sd(x, na.rm = TRUE) == 0, logical(1)))) {
  stop("At least one included feature has zero variance")
}

source_summary <- analysis_df %>%
  group_by(source_group) %>%
  summarise(
    n = n(),
    median_death_dsr = median(Death_DSR),
    q1_death_dsr = quantile(Death_DSR, 0.25),
    q3_death_dsr = quantile(Death_DSR, 0.75),
    min_death_dsr = min(Death_DSR),
    max_death_dsr = max(Death_DSR),
    .groups = "drop"
  )

cohort_summary <- tibble(
  metric = c(
    "Linked analysis cohort", "Assumed deaths", "Median Death_DSR", "Death_DSR IQR",
    "Death_DSR range", "Included NMR features", "Excluded ketone features", "Included CEC features"
  ),
  value = c(
    as.character(nrow(analysis_df)),
    as.character(sum(analysis_df$death_event)),
    sprintf("%.1f days", median(analysis_df$Death_DSR)),
    sprintf("%.1f–%.1f days", quantile(analysis_df$Death_DSR, 0.25), quantile(analysis_df$Death_DSR, 0.75)),
    sprintf("%d–%d days", min(analysis_df$Death_DSR), max(analysis_df$Death_DSR)),
    as.character(length(nmr_features)),
    paste(ketone_features, collapse = ", "),
    as.character(length(cec_features))
  )
)

categorical_summary <- bind_rows(
  analysis_df %>% count(variable = "Sex", level = as.character(gender), name = "n"),
  analysis_df %>% count(variable = "APOE genotype", level = as.character(apoe_group), name = "n"),
  analysis_df %>% count(variable = "Allocation", level = as.character(Treatment), name = "n"),
  analysis_df %>% count(variable = "AVx HDL source", level = as.character(source_group), name = "n"),
  analysis_df %>% count(variable = "Medication prescribed", level = as.character(medication), name = "n"),
  analysis_df %>% count(variable = "CVD", level = as.character(CVD), name = "n"),
  analysis_df %>% count(variable = "Stroke", level = as.character(stroke), name = "n"),
  analysis_df %>% count(variable = "Diabetes", level = as.character(diabetes), name = "n"),
  analysis_df %>% count(variable = "Hypertension", level = as.character(hypertension), name = "n"),
  analysis_df %>% count(variable = "Dyslipidemia", level = as.character(dyslipidemia), name = "n")
) %>%
  group_by(variable) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

medication_summary <- analysis_df %>%
  mutate(medication_name = if_else(is.na(AVx_Name_Med), "No statin-related medication recorded", AVx_Name_Med)) %>%
  count(medication, medication_name, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

feature_metadata <- feature_manifest %>%
  filter(analysis_status == "Included") %>%
  select(modality, feature, feature_label, feature_group, unit)

spearman_results <- map_dfr(c(nmr_features, cec_features), function(feature) {
  modality <- if_else(feature %in% nmr_features, "NMR", "CEC")
  x <- analysis_df[[feature]]
  keep <- is.finite(x) & is.finite(analysis_df$Death_DSR)
  test <- suppressWarnings(cor.test(x[keep], analysis_df$Death_DSR[keep], method = "spearman", exact = FALSE))
  tibble(
    modality = modality,
    feature = feature,
    n = sum(keep),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  left_join(feature_metadata, by = c("modality", "feature")) %>%
  group_by(modality) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    direction = case_when(rho > 0 ~ "Longer time to death", rho < 0 ~ "Shorter time to death", TRUE ~ "No direction")
  ) %>%
  ungroup() %>%
  arrange(modality, fdr, p_value)

base_covariates <- c("gender", "apoe_group", "Treatment", "source_group")
expanded_covariates <- c(
  base_covariates, "medication", "CVD", "stroke", "diabetes", "hypertension", "dyslipidemia"
)

fit_cox_feature <- function(feature, modality, model_name, covariates) {
  model_df <- analysis_df %>%
    select(Death_DSR, death_event, all_of(covariates), all_of(feature)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  model_df$z_feature <- as.numeric(scale(model_df[[feature]]))
  formula_text <- paste(
    "Surv(Death_DSR, death_event) ~ z_feature +",
    paste(covariates, collapse = " + ")
  )
  fit_warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      coxph(as.formula(formula_text), data = model_df, ties = "efron", x = TRUE, model = TRUE),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(tibble(
      modality = modality, feature = feature, model = model_name, n = nrow(model_df), events = sum(model_df$death_event),
      hr = NA_real_, conf_low = NA_real_, conf_high = NA_real_, log_hr = NA_real_, std_error = NA_real_,
      p_value = NA_real_, concordance = NA_real_, aic = NA_real_, ph_p_feature = NA_real_, ph_p_global = NA_real_,
      status = paste("ERROR:", conditionMessage(fit))
    ))
  }
  coef_row <- summary(fit)$coefficients["z_feature", , drop = FALSE]
  beta <- unname(coef_row[1, "coef"])
  se <- unname(coef_row[1, "se(coef)"])
  ph <- tryCatch(cox.zph(fit, transform = "km", terms = TRUE, singledf = TRUE), error = function(e) NULL)
  ph_feature <- if (is.null(ph) || !"z_feature" %in% rownames(ph$table)) NA_real_ else ph$table["z_feature", "p"]
  ph_global <- if (is.null(ph) || !"GLOBAL" %in% rownames(ph$table)) NA_real_ else ph$table["GLOBAL", "p"]
  concordance <- summary(fit)$concordance
  concordance_value <- if (length(concordance) >= 1) unname(concordance[1]) else NA_real_
  tibble(
    modality = modality,
    feature = feature,
    model = model_name,
    n = fit$n,
    events = fit$nevent,
    hr = exp(beta),
    conf_low = exp(beta - 1.96 * se),
    conf_high = exp(beta + 1.96 * se),
    log_hr = beta,
    std_error = se,
    p_value = unname(coef_row[1, "Pr(>|z|)"]),
    concordance = concordance_value,
    aic = extractAIC(fit)[2],
    ph_p_feature = ph_feature,
    ph_p_global = ph_global,
    status = if_else(length(fit_warnings) == 0, "OK", paste(unique(fit_warnings), collapse = " | "))
  )
}

model_specifications <- list(
  base = base_covariates,
  expanded = expanded_covariates
)

cox_results <- map_dfr(names(model_specifications), function(model_name) {
  covariates <- model_specifications[[model_name]]
  map_dfr(c(nmr_features, cec_features), function(feature) {
    fit_cox_feature(
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
    ph_fdr_feature = p.adjust(ph_p_feature, method = "BH"),
    direction = case_when(hr > 1 ~ "Earlier death / higher hazard", hr < 1 ~ "Later death / lower hazard", TRUE ~ "No direction"),
    significant_fdr = !is.na(fdr) & fdr < fdr_alpha
  ) %>%
  ungroup() %>%
  arrange(modality, model, fdr, p_value)

# The vendor codebook defines nine metabolic-vulnerability (MVX) measures.
# Four MVX-prefixed fields are exact aliases of canonical NMR variables already
# tested above; map them back here so the biological panel is complete without
# duplicating observations or altering the original 45-feature FDR family.
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

metabolic_vulnerability_results <- cox_results %>%
  filter(modality == "NMR", feature %in% metabolic_vulnerability_map$analysis_feature) %>%
  left_join(metabolic_vulnerability_map, by = c("feature" = "analysis_feature")) %>%
  arrange(display_order, model) %>%
  select(
    mvx_variable, analysis_feature = feature, display_label, display_order,
    model, n, events, hr, conf_low, conf_high, log_hr, std_error,
    p_value, fdr, ph_p_feature, ph_fdr_feature, direction, significant_fdr
  )

stopifnot(
  nrow(metabolic_vulnerability_results) == 18L,
  dplyr::n_distinct(metabolic_vulnerability_results$mvx_variable) == 9L,
  !any(metabolic_vulnerability_results$analysis_feature %in% ketone_features)
)

fit_hdl_profile <- function(model_name, covariates) {
  hdl_features <- paste0("H", 1:7, "P")
  model_df <- analysis_df %>%
    select(Death_DSR, death_event, all_of(covariates), all_of(hdl_features)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  for (feature in hdl_features) {
    model_df[[paste0("z_", feature)]] <- as.numeric(scale(model_df[[feature]]))
  }
  reduced_formula <- as.formula(paste("Surv(Death_DSR, death_event) ~", paste(covariates, collapse = " + ")))
  full_formula <- as.formula(paste(
    "Surv(Death_DSR, death_event) ~",
    paste(c(covariates, paste0("z_", hdl_features)), collapse = " + ")
  ))
  reduced_fit <- coxph(reduced_formula, data = model_df, ties = "efron", x = TRUE)
  full_fit <- coxph(full_formula, data = model_df, ties = "efron", x = TRUE)
  lrt <- 2 * (as.numeric(logLik(full_fit)) - as.numeric(logLik(reduced_fit)))
  df_diff <- attr(logLik(full_fit), "df") - attr(logLik(reduced_fit), "df")
  global <- tibble(
    model = model_name,
    n = full_fit$n,
    events = full_fit$nevent,
    degrees_freedom = df_diff,
    likelihood_ratio_chisq = lrt,
    p_value = pchisq(lrt, df = df_diff, lower.tail = FALSE),
    reduced_aic = extractAIC(reduced_fit)[2],
    full_aic = extractAIC(full_fit)[2]
  )
  coefficients <- map_dfr(hdl_features, function(feature) {
    term <- paste0("z_", feature)
    row <- summary(full_fit)$coefficients[term, , drop = FALSE]
    beta <- unname(row[1, "coef"])
    se <- unname(row[1, "se(coef)"])
    tibble(
      model = model_name,
      feature = feature,
      hr = exp(beta),
      conf_low = exp(beta - 1.96 * se),
      conf_high = exp(beta + 1.96 * se),
      p_value = unname(row[1, "Pr(>|z|)"])
    )
  })
  list(global = global, coefficients = coefficients)
}

hdl_profile_fits <- imap(model_specifications, ~ fit_hdl_profile(.y, .x))
hdl_profile_global <- bind_rows(map(hdl_profile_fits, "global")) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))
hdl_profile_coefficients <- bind_rows(map(hdl_profile_fits, "coefficients")) %>%
  left_join(feature_metadata %>% select(feature, feature_label, unit), by = "feature")

safe_biomarker_names <- make.names(c(nmr_features, cec_features), unique = TRUE)
biomarker_name_map <- tibble(
  feature = c(nmr_features, cec_features),
  matrix_feature = paste0("bio__", safe_biomarker_names),
  modality = c(rep("NMR", length(nmr_features)), rep("CEC", length(cec_features)))
)
biomarker_matrix <- analysis_df %>%
  select(all_of(c(nmr_features, cec_features))) %>%
  as.matrix()
storage.mode(biomarker_matrix) <- "double"
biomarker_matrix <- scale(biomarker_matrix)
colnames(biomarker_matrix) <- biomarker_name_map$matrix_feature

covariate_matrix <- model.matrix(
  ~ gender + apoe_group + Treatment + source_group + medication + CVD + stroke + diabetes + hypertension + dyslipidemia - 1,
  data = analysis_df
)
penalized_x <- cbind(covariate_matrix, biomarker_matrix)
penalty_factor <- c(rep(0, ncol(covariate_matrix)), rep(1, ncol(biomarker_matrix)))
penalized_y <- Surv(analysis_df$Death_DSR, analysis_df$death_event)

cv_fit <- cv.glmnet(
  x = penalized_x,
  y = penalized_y,
  family = "cox",
  alpha = 0.5,
  nfolds = 10,
  type.measure = "C",
  standardize = FALSE,
  penalty.factor = penalty_factor,
  cox.ties = "efron",
  keep = TRUE
)

extract_penalized <- function(lambda_name, lambda_value) {
  coefficients <- as.matrix(coef(cv_fit$glmnet.fit, s = lambda_value))[, 1]
  biomarker_coefficients <- coefficients[biomarker_name_map$matrix_feature]
  tibble(
    lambda_rule = lambda_name,
    matrix_feature = names(biomarker_coefficients),
    coefficient = as.numeric(biomarker_coefficients)
  ) %>%
    left_join(biomarker_name_map, by = "matrix_feature") %>%
    left_join(feature_metadata, by = c("modality", "feature")) %>%
    filter(abs(coefficient) > 1e-12) %>%
    mutate(
      hazard_ratio_per_sd = exp(coefficient),
      direction = if_else(coefficient > 0, "Earlier death / higher hazard", "Later death / lower hazard")
    ) %>%
    arrange(desc(abs(coefficient)))
}

elastic_net_selected <- bind_rows(
  extract_penalized("lambda.min", cv_fit$lambda.min),
  extract_penalized("lambda.1se", cv_fit$lambda.1se)
)

index_min <- which.min(abs(cv_fit$lambda - cv_fit$lambda.min))
index_1se <- which.min(abs(cv_fit$lambda - cv_fit$lambda.1se))
elastic_net_summary <- tibble(
  lambda_rule = c("lambda.min", "lambda.1se"),
  lambda = c(cv_fit$lambda.min, cv_fit$lambda.1se),
  cross_validated_c_index = c(cv_fit$cvm[index_min], cv_fit$cvm[index_1se]),
  cross_validated_c_index_se = c(cv_fit$cvsd[index_min], cv_fit$cvsd[index_1se]),
  selected_biomarkers = c(
    sum(elastic_net_selected$lambda_rule == "lambda.min"),
    sum(elastic_net_selected$lambda_rule == "lambda.1se")
  ),
  alpha = 0.5,
  folds = 10L,
  seed = 20260902L
)

figure_palette <- c(
  ink = "#222222",
  grey_dark = "#666666",
  grey_mid = "#A6A6A6",
  grey_light = "#D9D9D9",
  grey_pale = "#F2F2F2",
  blue = "#2C6EAA",
  blue_light = "#A9C9E2",
  red = "#C4473A",
  red_pale = "#F8E6E3"
)
figure_width_mm = 183

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
    if (!nzchar(gs_bin)) {
      stop("Ghostscript is required for the editable PDF fallback when Cairo is unavailable")
    }
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
    if (!identical(gs_status, 0L) || !file.exists(pdf_path)) {
      stop("Ghostscript PDF export failed for ", stem)
    }
  }

  ragg::agg_tiff(tiff_path, width = width_in, height = height_in, units = "in", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in", res = 300, background = "white")
  print(plot)
  grDevices::dev.off()
}

death_median <- median(analysis_df$Death_DSR)
death_distribution_plot <- ggplot(analysis_df, aes(x = Death_DSR)) +
  geom_histogram(
    binwidth = 150, boundary = 0,
    fill = figure_palette[["grey_light"]], colour = "white", linewidth = 0.25
  ) +
  geom_vline(xintercept = death_median, colour = figure_palette[["blue"]], linewidth = 0.7) +
  annotate(
    "text", x = death_median + 45, y = Inf,
    label = paste0("Median = ", format(round(death_median), big.mark = ","), " d"),
    vjust = 1.15, hjust = 0, colour = figure_palette[["blue"]], size = 2.25, family = "Helvetica"
  ) +
  scale_x_continuous(breaks = seq(1000, 2500, 500), labels = label_number(big.mark = ","), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Days since randomization", y = "Participants") +
  theme_nature()
save_plot(death_distribution_plot, "death_dsr_distribution", 89, 62)

source_counts <- analysis_df %>%
  mutate(source_group = factor(source_group, levels = c("Yr 2", "Yr 3", "Yr 5+"))) %>%
  count(source_group, name = "n") %>%
  arrange(source_group)
source_label_levels <- paste0(source_counts$source_group, "\nn = ", source_counts$n)
source_plot_data <- analysis_df %>%
  mutate(source_group = factor(source_group, levels = c("Yr 2", "Yr 3", "Yr 5+"))) %>%
  left_join(source_counts, by = "source_group") %>%
  mutate(source_label = factor(
    paste0(source_group, "\nn = ", n),
    levels = source_label_levels
  ))
source_colours <- setNames(c("#9ABFDA", "#6E9CC0", "#3E6C8F"), source_label_levels)
source_fills <- setNames(c("#DCEAF4", "#C5D8E7", "#B5C8D7"), source_label_levels)
source_plot <- ggplot(source_plot_data, aes(x = source_label, y = Death_DSR)) +
  geom_boxplot(
    aes(fill = source_label), width = 0.42, outlier.shape = NA,
    colour = figure_palette[["ink"]], linewidth = 0.4
  ) +
  geom_point(
    aes(colour = source_label),
    position = position_jitter(width = 0.10, height = 0, seed = 20260902),
    size = 0.72, stroke = 0
  ) +
  scale_fill_manual(values = source_fills, guide = "none") +
  scale_colour_manual(values = source_colours, guide = "none") +
  scale_y_continuous(breaks = seq(1000, 2500, 500), labels = label_number(big.mark = ","), expand = expansion(mult = c(0.03, 0.04))) +
  labs(x = "Sample-source group", y = "Time to death (days)") +
  theme_nature()
save_plot(source_plot, "death_dsr_by_source", 89, 67)

spearman_plot_data <- spearman_results %>%
  filter(modality == "NMR") %>%
  slice_min(order_by = fdr, n = 15, with_ties = FALSE) %>%
  arrange(rho) %>%
  mutate(
    feature = factor(feature, levels = feature),
    highlight = feature == "H4P"
  )
spearman_plot <- ggplot(spearman_plot_data, aes(x = rho, y = feature)) +
  geom_vline(xintercept = 0, colour = figure_palette[["grey_mid"]], linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = rho, yend = feature), colour = figure_palette[["grey_mid"]], linewidth = 0.4) +
  geom_point(aes(fill = highlight), shape = 21, colour = "white", stroke = 0.25, size = 2.0) +
  scale_fill_manual(values = c(`TRUE` = figure_palette[["red"]], `FALSE` = figure_palette[["grey_dark"]]), guide = "none") +
  scale_x_continuous(limits = c(-0.28, 0.22), breaks = seq(-0.2, 0.2, 0.1)) +
  labs(x = "Spearman rho with time to death", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.4))
save_plot(spearman_plot, "spearman_nmr_top", 89, 95)

nmr_feature_rank <- cox_results %>%
  filter(modality == "NMR") %>%
  group_by(feature) %>%
  summarise(rank_fdr = min(fdr), .groups = "drop") %>%
  slice_min(order_by = rank_fdr, n = 12, with_ties = FALSE) %>%
  arrange(rank_fdr)
nmr_forest_data <- cox_results %>%
  filter(modality == "NMR", feature %in% nmr_feature_rank$feature) %>%
  mutate(
    model = recode(model, base = "Base", expanded = "Expanded"),
    feature = factor(feature, levels = rev(nmr_feature_rank$feature))
  )
nmr_forest_plot <- ggplot(nmr_forest_data, aes(x = hr, y = feature, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.52), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.52), size = 1.9, stroke = 0.5) +
  geom_text(
    data = nmr_forest_data %>% filter(feature == "H4P", model == "Base"),
    aes(x = conf_high * 1.025, y = feature, label = "q = 0.027"),
    hjust = 0, vjust = 0.72, size = 2.15, family = "Helvetica",
    colour = figure_palette[["red"]], inherit.aes = FALSE
  ) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(
    limits = c(0.68, 1.82), breaks = c(0.7, 1.0, 1.4, 1.8),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(x = "Hazard ratio per 1-SD higher biomarker (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.3), legend.position = "top")
save_plot(nmr_forest_plot, "nmr_cox_top_forest", 183, 112)

mvx_display_order <- metabolic_vulnerability_map %>%
  arrange(display_order) %>%
  pull(display_label)
mvx_forest_data <- metabolic_vulnerability_results %>%
  mutate(
    model = recode(model, base = "Base", expanded = "Expanded"),
    display_label = factor(display_label, levels = rev(mvx_display_order))
  )
mvx_forest_plot <- ggplot(mvx_forest_data, aes(x = hr, y = display_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.50), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.50), size = 1.9, stroke = 0.5) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(
    limits = c(0.72, 1.48), breaks = c(0.8, 1.0, 1.2, 1.4),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(x = "Hazard ratio per 1-SD higher biomarker (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.3), legend.position = "top")
save_plot(mvx_forest_plot, "metabolic_vulnerability_cox_forest", 183, 100)

nmr_stability_data <- cox_results %>%
  filter(modality == "NMR") %>%
  select(feature, model, log_hr, p_value, fdr) %>%
  pivot_wider(names_from = model, values_from = c(log_hr, p_value, fdr)) %>%
  mutate(
    base_log2_hr = log_hr_base / log(2),
    expanded_log2_hr = log_hr_expanded / log(2),
    evidence = case_when(
      fdr_base < fdr_alpha ~ "Base q < 0.05",
      p_value_base < 0.05 ~ "Base P < 0.05",
      TRUE ~ "Other"
    )
  )
# Cox hazard ratios and glmnet penalties must be strictly positive before log-scale display.
stopifnot(all(cox_results$hr > 0), all(cv_fit$lambda > 0))
stability_limit <- max(abs(c(nmr_stability_data$base_log2_hr, nmr_stability_data$expanded_log2_hr))) * 1.30
nmr_volcano_plot <- ggplot(nmr_stability_data, aes(x = base_log2_hr, y = expanded_log2_hr)) +
  geom_abline(slope = 1, intercept = 0, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_hline(yintercept = 0, colour = figure_palette[["grey_light"]], linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = figure_palette[["grey_light"]], linewidth = 0.3) +
  geom_point(aes(fill = evidence), shape = 21, colour = "white", stroke = 0.25, size = 2.05) +
  geom_text(
    data = nmr_stability_data %>% filter(feature == "H4P"),
    aes(label = feature), colour = figure_palette[["red"]],
    size = 2.25, family = "Helvetica", nudge_x = -0.015, nudge_y = 0.075,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "Base q < 0.05" = figure_palette[["red"]],
    "Base P < 0.05" = figure_palette[["grey_dark"]],
    "Other" = figure_palette[["grey_light"]]
  )) +
  scale_x_continuous(limits = c(-stability_limit, stability_limit), breaks = pretty_breaks(4)) +
  scale_y_continuous(limits = c(-stability_limit, stability_limit), breaks = pretty_breaks(4)) +
  coord_fixed() +
  labs(x = "Base-model log2(HR per SD)", y = "Expanded-model log2(HR per SD)") +
  theme_nature() +
  theme(legend.position = "top")
save_plot(nmr_volcano_plot, "nmr_cox_volcano", 89, 85)

cec_forest_data <- cox_results %>%
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
cec_forest_plot <- ggplot(cec_forest_data, aes(x = hr, y = feature_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = figure_palette[["grey_mid"]], linetype = 2, linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0,
    position = position_dodge(width = 0.42), linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.42), size = 1.9, stroke = 0.5) +
  scale_colour_manual(values = c("Base" = figure_palette[["ink"]], "Expanded" = figure_palette[["blue"]])) +
  scale_shape_manual(values = c("Base" = 21, "Expanded" = 19)) +
  scale_x_log10(limits = c(0.72, 1.32), breaks = c(0.8, 1.0, 1.2), labels = label_number(accuracy = 0.1)) +
  labs(x = "Hazard ratio per 1-SD higher CEC (95% CI)", y = NULL) +
  theme_nature() +
  theme(axis.text.y = element_text(size = 6.4), legend.position = "top")
save_plot(cec_forest_plot, "cec_cox_forest", 89, 60)

coef_path <- as.matrix(cv_fit$glmnet.fit$beta)
biomarker_path_rows <- intersect(rownames(coef_path), biomarker_name_map$matrix_feature)
selected_path_counts <- colSums(abs(coef_path[biomarker_path_rows, , drop = FALSE]) > 1e-12)
fit_lambda <- cv_fit$glmnet.fit$lambda
cv_plot_data <- tibble(
  lambda = cv_fit$lambda,
  log10_lambda = log10(cv_fit$lambda),
  c_index = cv_fit$cvm,
  c_index_se = cv_fit$cvsd,
  selected_biomarkers = vapply(
    cv_fit$lambda,
    function(value) selected_path_counts[which.min(abs(fit_lambda - value))],
    numeric(1)
  )
)
elastic_plot <- ggplot(cv_plot_data, aes(x = log10_lambda, y = c_index)) +
  geom_ribbon(
    aes(ymin = c_index - c_index_se, ymax = c_index + c_index_se),
    fill = "#E6E6E6"
  ) +
  geom_line(colour = figure_palette[["ink"]], linewidth = 0.5) +
  geom_point(colour = figure_palette[["ink"]], size = 0.65) +
  geom_vline(xintercept = log10(cv_fit$lambda.min), colour = figure_palette[["red"]], linewidth = 0.55) +
  geom_vline(xintercept = log10(cv_fit$lambda.1se), colour = figure_palette[["blue"]], linetype = 2, linewidth = 0.55) +
  annotate(
    "text", x = log10(cv_fit$lambda.min), y = min(cv_plot_data$c_index - cv_plot_data$c_index_se),
    label = "lambda.min: H4P", colour = figure_palette[["red"]], angle = 90,
    hjust = -0.04, vjust = -0.35, size = 2.15, family = "Helvetica"
  ) +
  annotate(
    "text", x = log10(cv_fit$lambda.1se), y = min(cv_plot_data$c_index - cv_plot_data$c_index_se),
    label = "lambda.1se: none", colour = figure_palette[["blue"]], angle = 90,
    hjust = -0.04, vjust = 1.15, size = 2.15, family = "Helvetica"
  ) +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), expand = expansion(mult = c(0.04, 0.07))) +
  labs(x = "log10(lambda)", y = "Cross-validated C-index") +
  theme_nature()
save_plot(elastic_plot, "elastic_net_selected_features", 89, 67)

write_csv(analysis_df, file.path(results_dir, "analysis_dataset.csv"))
write_csv(feature_manifest, file.path(results_dir, "feature_manifest.csv"))
write_csv(cohort_summary, file.path(results_dir, "cohort_summary.csv"))
write_csv(categorical_summary, file.path(results_dir, "covariate_distribution.csv"))
write_csv(medication_summary, file.path(results_dir, "medication_distribution.csv"))
write_csv(source_summary, file.path(results_dir, "death_dsr_by_source_summary.csv"))
write_csv(spearman_results, file.path(results_dir, "spearman_results.csv"))
write_csv(cox_results, file.path(results_dir, "cox_single_feature_results.csv"))
write_csv(metabolic_vulnerability_results, file.path(results_dir, "metabolic_vulnerability_results.csv"))
write_csv(hdl_profile_global, file.path(results_dir, "hdl_profile_global_tests.csv"))
write_csv(hdl_profile_coefficients, file.path(results_dir, "hdl_profile_joint_coefficients.csv"))
write_csv(elastic_net_summary, file.path(results_dir, "elastic_net_summary.csv"))
write_csv(elastic_net_selected, file.path(results_dir, "elastic_net_selected_features.csv"))

analysis_assumptions <- tibble(
  item = c(
    "Death endpoint", "Time origin", "Sampling time", "Age", "Sampling-source adjustment",
    "Medication", "Ketone panel", "Multiplicity", "Interpretation"
  ),
  specification = c(
    "All 174 participants are assigned death_event = 1; Death_DSR is treated as time from randomization to death.",
    "Randomization (day 0).",
    "Exact specimen-collection DSR was unavailable; delayed entry could not be implemented.",
    "Participant age and sample age were not included.",
    "AVx_HDL_Source was included in both models as Yr 2, Yr 3, or Yr 5+.",
    "AVx_MedPres equal to 1 was coded as medication = 1; missing was coded as 0. Drug names were descriptive only.",
    "KetBod, B-HB, AcAc, and Acetone were excluded because of alcohol contamination.",
    "Benjamini-Hochberg FDR was controlled separately for NMR and CEC within each model.",
    "Hazard ratios describe earlier versus later death within this presumed-decedent cohort and are not population absolute-risk estimates."
  )
)
write_csv(analysis_assumptions, file.path(results_dir, "analysis_assumptions.csv"))

table_files <- c(
  "analysis_dataset.csv", "feature_manifest.csv", "cohort_summary.csv", "covariate_distribution.csv",
  "medication_distribution.csv", "death_dsr_by_source_summary.csv", "spearman_results.csv",
  "cox_single_feature_results.csv", "metabolic_vulnerability_results.csv",
  "hdl_profile_global_tests.csv", "hdl_profile_joint_coefficients.csv",
  "elastic_net_summary.csv", "elastic_net_selected_features.csv", "analysis_assumptions.csv"
)
all_plot_files <- list.files(plots_dir, full.names = FALSE)
plot_files <- all_plot_files[tolower(tools::file_ext(all_plot_files)) %in% c("png", "pdf", "svg", "tiff")]
output_manifest <- bind_rows(
  tibble(file = table_files, type = "table", description = "Authoritative numerical analysis output"),
  tibble(file = file.path("plots", plot_files), type = "figure", description = "Report figure"),
  tibble(
    file = c(
      "nmr_cec_death_dsr_report.html", "nmr_cec_death_dsr_report.Rmd", "session_info.txt",
      "figure_contract.md", "figure_qa_notes.md"
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

report_input <- file.path(results_dir, "nmr_cec_death_dsr_report.Rmd")
if (!file.exists(report_input)) {
  stop("Report template not found: ", report_input)
}
rmarkdown::render(
  input = report_input,
  output_file = "nmr_cec_death_dsr_report.html",
  output_dir = results_dir,
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)

cat("Analysis complete. Report:", file.path(results_dir, "nmr_cec_death_dsr_report.html"), "\n")
