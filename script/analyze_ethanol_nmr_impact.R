options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(purrr)
  library(readr)
  library(readxl)
  library(tidyr)
})

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (is.na(script_path)) {
  script_dir <- normalizePath(file.path(getwd(), "script"), winslash = "/", mustWork = FALSE)
} else {
  script_dir <- dirname(script_path)
}
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
clean_dir <- file.path(project_root, "clean_data")
results_root <- file.path(project_root, "results")
results_dir <- file.path(results_root, "ethanol_nmr_impact_analysis")
plots_dir <- file.path(results_dir, "plots")
if (dir.exists(results_dir)) {
  unlink(results_dir, recursive = TRUE, force = TRUE)
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

threshold_200 <- 200
threshold_400 <- 400
fdr_alpha <- 0.05
x_var <- "plasma_spectrum__ch2_ethanol"
x_unit <- "mg/dL"

ethanol_df <- read_csv(file.path(clean_dir, "ASPREE_ethanol_contamination_clean.csv"), show_col_types = FALSE)
lipid_df <- read_csv(file.path(clean_dir, "ASPREE_merged_lipoprofile_clean.csv"), show_col_types = FALSE)
metadata_df <- read_excel(file.path(project_root, "raw_data", "subject_metadata.xlsx"))
codebook_df <- read_excel(file.path(clean_dir, "ASPREE_merged_lipoprofile_codebook.xlsx"))

lipid_df <- lipid_df %>%
  rename(
    subject_id = `Subject ID`,
    labcorp_accession = `Labcorp Accession Number`,
    collection_date = `Collection Date`,
    comment = Comment
  )

metadata_df <- metadata_df %>%
  rename(
    subject_id = subject_id,
    gender = gender,
    dementia_diagnosis = dementia_diagnosis,
    apoe_group = apoe_group,
    avx_hdl_mgdl = AVx_hdl_mgDL,
    hdl_category = hdl_category
  )

analysis_df <- lipid_df %>%
  inner_join(
    ethanol_df %>% select(subject_id, all_of(x_var)),
    by = "subject_id"
  ) %>%
  left_join(metadata_df, by = "subject_id") %>%
  mutate(
    contamination_status_200 = if_else(.data[[x_var]] >= threshold_200, "Contaminated", "Clean"),
    contamination_status_400 = if_else(.data[[x_var]] >= threshold_400, "Contaminated", "Clean")
  )

sample_set_display <- c(
  all_samples = "All samples",
  clean_threshold_200 = "Clean samples (< 200 mg/dL)",
  clean_threshold_400 = "Clean samples (< 400 mg/dL)"
)

codebook_lookup <- codebook_df %>%
  transmute(
    variable = abbreviated_variable_name,
    variable_name = variable_name,
    variable_group = variable_group_description,
    unit = measurement_units
  )

excluded_vars <- c(
  "subject_id",
  "labcorp_accession",
  "collection_date",
  "comment"
)

numeric_lipid_vars <- names(analysis_df)[vapply(analysis_df, is.numeric, logical(1))]
numeric_lipid_vars <- setdiff(numeric_lipid_vars, c(x_var, "avx_hdl_mgdl", excluded_vars))

describe_var <- function(var_name) {
  hit <- codebook_lookup %>% filter(variable == var_name)
  if (nrow(hit) == 0) {
    tibble(variable = var_name, variable_name = var_name, variable_group = NA_character_, unit = NA_character_)
  } else {
    hit[1, ]
  }
}

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(x) > 5000 || dplyr::n_distinct(x) < 3) {
    return(list(p_value = NA_real_, is_normal = NA))
  }
  p <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
  list(p_value = p, is_normal = ifelse(is.na(p), NA, p > 0.05))
}

analyze_one_round <- function(df, round_label) {
  x_norm <- safe_shapiro(df[[x_var]])
  x_normality <- tibble(
    sample_set = round_label,
    variable = x_var,
    variable_name = x_var,
    variable_group = "Ethanol contamination",
    unit = x_unit,
    n = sum(is.finite(df[[x_var]])),
    shapiro_p = x_norm$p_value,
    is_normal = x_norm$is_normal
  )

  results <- map_dfr(numeric_lipid_vars, function(var_name) {
    sub <- df %>%
      select(all_of(c(x_var, var_name))) %>%
      filter(if_all(everything(), ~ is.finite(.x)))

    var_meta <- describe_var(var_name)
    y_norm <- safe_shapiro(sub[[var_name]])
    method <- if (!is.na(x_norm$is_normal) && !is.na(y_norm$is_normal) && x_norm$is_normal && y_norm$is_normal) "pearson" else "spearman"

    cor_res <- if (nrow(sub) >= 3 && dplyr::n_distinct(sub[[x_var]]) >= 3 && dplyr::n_distinct(sub[[var_name]]) >= 3) {
      tryCatch(
        suppressWarnings(cor.test(sub[[x_var]], sub[[var_name]], method = method, exact = FALSE)),
        error = function(e) NULL
      )
    } else {
      NULL
    }

    estimate_value <- if (is.null(cor_res)) NA_real_ else unname(cor_res$estimate[[1]])
    p_value <- if (is.null(cor_res)) NA_real_ else cor_res$p.value

    tibble(
      sample_set = round_label,
      variable = var_name,
      variable_name = var_meta$variable_name,
      variable_group = var_meta$variable_group,
      unit = var_meta$unit,
      n = nrow(sub),
      x_shapiro_p = x_norm$p_value,
      x_is_normal = x_norm$is_normal,
      y_shapiro_p = y_norm$p_value,
      y_is_normal = y_norm$is_normal,
      method = method,
      estimate = estimate_value,
      p_value = p_value
    )
  }) %>%
    mutate(
      fdr = p.adjust(p_value, method = "fdr"),
      significant_fdr = !is.na(fdr) & fdr < fdr_alpha
    )

  y_normality <- results %>%
    transmute(
      sample_set = sample_set,
      variable = variable,
      variable_name = variable_name,
      variable_group = variable_group,
      unit = unit,
      n = n,
      shapiro_p = y_shapiro_p,
      is_normal = y_is_normal
    )

  list(
    correlations = results,
    normality = bind_rows(x_normality, y_normality)
  )
}

get_sample_df <- function(sample_set_label) {
  if (sample_set_label == "all_samples") {
    analysis_df
  } else if (sample_set_label == "clean_threshold_200") {
    filter(analysis_df, contamination_status_200 == "Clean")
  } else if (sample_set_label == "clean_threshold_400") {
    filter(analysis_df, contamination_status_400 == "Clean")
  } else {
    stop("Unknown sample set: ", sample_set_label)
  }
}

sample_sets <- c("all_samples", "clean_threshold_200", "clean_threshold_400")
round_results <- setNames(
  lapply(sample_sets, function(label) analyze_one_round(get_sample_df(label), label)),
  sample_sets
)

normality_results <- bind_rows(lapply(round_results, `[[`, "normality"))
correlation_results <- bind_rows(lapply(round_results, `[[`, "correlations"))

write_csv(normality_results, file.path(results_dir, "normality_results.csv"))
write_csv(filter(correlation_results, sample_set == "all_samples"), file.path(results_dir, "correlation_results_all_samples.csv"))
write_csv(filter(correlation_results, sample_set == "clean_threshold_200"), file.path(results_dir, "correlation_results_clean_threshold_200.csv"))
write_csv(filter(correlation_results, sample_set == "clean_threshold_400"), file.path(results_dir, "correlation_results_clean_threshold_400.csv"))
write_csv(analysis_df, file.path(results_dir, "analysis_dataset_with_metadata.csv"))

make_distribution <- function(df, sample_set_label) {
  vars <- c("dementia_diagnosis", "gender", "apoe_group", "hdl_category")
  map_dfr(vars, function(var_name) {
    df %>%
      count(value = .data[[var_name]], name = "n") %>%
      mutate(
        sample_set = sample_set_label,
        variable = var_name,
        percent = 100 * n / sum(n)
      ) %>%
      select(sample_set, variable, value, n, percent)
  })
}

distribution_summary <- bind_rows(
  make_distribution(get_sample_df("all_samples"), "all_samples"),
  make_distribution(get_sample_df("clean_threshold_200"), "clean_threshold_200"),
  make_distribution(get_sample_df("clean_threshold_400"), "clean_threshold_400")
)
write_csv(distribution_summary, file.path(results_dir, "metadata_distribution_summary.csv"))

sample_counts <- tibble(
  sample_set = sample_sets,
  n = c(
    nrow(get_sample_df("all_samples")),
    nrow(get_sample_df("clean_threshold_200")),
    nrow(get_sample_df("clean_threshold_400"))
  )
)
write_csv(sample_counts, file.path(results_dir, "sample_counts.csv"))

make_metadata_plot <- function(sample_set_label) {
  plot_df <- distribution_summary %>%
    filter(sample_set == sample_set_label) %>%
    mutate(
      variable = recode(
        variable,
        dementia_diagnosis = "Dementia vs Normal",
        gender = "Female vs Male",
        apoe_group = "APOE genotype distribution",
        hdl_category = "HDL-category"
      ),
      value = factor(value, levels = unique(value))
    )

  file_name <- paste0("metadata_distribution_", sample_set_label, ".png")
  file_path <- file.path(plots_dir, file_name)

  p <- ggplot(plot_df, aes(x = value, y = n, fill = value)) +
    geom_col(width = 0.75, show.legend = FALSE) +
    geom_text(aes(label = paste0(n, "\n", sprintf("%.1f%%", percent))), vjust = -0.2, size = 3.3) +
    facet_wrap(~ variable, scales = "free_x", ncol = 2) +
    labs(
      title = paste0(sample_set_display[[sample_set_label]], ": metadata distribution"),
      x = NULL,
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )

  ggsave(file_path, plot = p, width = 12, height = 8, dpi = 300)
  file_name
}

metadata_plot_manifest <- tibble(
  sample_set = sample_sets,
  plot_file = vapply(sample_sets, make_metadata_plot, character(1))
)
write_csv(metadata_plot_manifest, file.path(results_dir, "metadata_plot_manifest.csv"))

plot_sig_result <- function(result_row) {
  sample_set_label <- result_row$sample_set[[1]]
  var_name <- result_row$variable[[1]]
  cor_method <- result_row$method[[1]]
  estimate <- result_row$estimate[[1]]
  p_value <- result_row$p_value[[1]]
  fdr <- result_row$fdr[[1]]
  significant_fdr <- result_row$significant_fdr[[1]]
  y_label_name <- result_row$variable_name[[1]] %||% var_name
  y_unit <- result_row$unit[[1]]

  plot_df <- get_sample_df(sample_set_label) %>%
    transmute(
      ethanol = .data[[x_var]],
      y = .data[[var_name]]
    ) %>%
    filter(if_all(everything(), ~ is.finite(.x)))

  if (nrow(plot_df) < 3) {
    return(NA_character_)
  }

  corr_label <- if (cor_method == "pearson") "Pearson r" else "Spearman rho"
  annotation <- paste0(
    corr_label, " = ", sprintf("%.3f", estimate),
    "\nP = ", format.pval(p_value, digits = 3, eps = 1e-4),
    "\nFDR = ", format.pval(fdr, digits = 3, eps = 1e-4),
    "\nN = ", nrow(plot_df),
    ifelse(significant_fdr, "\nFDR-significant", "")
  )

  file_name <- paste0(sample_set_label, "__", gsub("[^A-Za-z0-9_]+", "_", var_name), ".png")
  file_path <- file.path(plots_dir, file_name)

  p <- ggplot(plot_df, aes(x = ethanol, y = y)) +
    geom_point(alpha = 0.75, color = "#1F4E78", size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#B22222", fill = "#F4A6A6") +
    annotate(
      "label",
      x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
      label = annotation,
      size = 3.5,
      label.size = 0.2
    ) +
    labs(
      title = paste0(
        sample_set_display[[sample_set_label]], ": ", y_label_name, " vs plasma spectrum CH2 ethanol",
        ifelse(significant_fdr, " [FDR-significant]", "")
      ),
      x = paste0("plasma_spectrum__ch2_ethanol (", x_unit, ")"),
      y = ifelse(is.na(y_unit) || y_unit == "", y_label_name, paste0(y_label_name, " (", y_unit, ")"))
    ) +
    theme_minimal(base_size = 12)

  ggsave(file_path, plot = p, width = 7.2, height = 5.2, dpi = 300)
  file_name
}

sig_results <- correlation_results %>%
  filter(!is.na(p_value), p_value < 0.05) %>%
  arrange(sample_set, fdr, p_value)

plot_manifest <- if (nrow(sig_results) > 0) {
  sig_results %>%
    mutate(plot_file = map_chr(split(., seq_len(n())), plot_sig_result))
} else {
  sig_results %>% mutate(plot_file = character())
}

if (nrow(plot_manifest) > 0) {
  write_csv(plot_manifest, file.path(results_dir, "significant_results_with_plots.csv"))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_pct <- function(x) {
  ifelse(is.na(x), "", paste0(formatC(x, format = "f", digits = 1), "%"))
}

result_summary_table <- function(sample_set_label) {
  df <- correlation_results %>% filter(sample_set == sample_set_label)
  nominal <- df %>%
    filter(!is.na(p_value), p_value < 0.05) %>%
    arrange(p_value, fdr) %>%
    select(variable, variable_name, variable_group, unit, method, estimate, p_value, fdr)
  fdr_pass <- df %>%
    filter(significant_fdr) %>%
    arrange(fdr, p_value) %>%
    select(variable, variable_name, variable_group, unit, method, estimate, p_value, fdr)
  list(nominal = nominal, fdr = fdr_pass)
}

correlation_report_table <- function(sample_set_label) {
  correlation_results %>%
    filter(sample_set == sample_set_label) %>%
    mutate(
      estimate = fmt_num(estimate, 3),
      p_value = ifelse(is.na(p_value), "", format.pval(p_value, digits = 3, eps = 1e-4)),
      fdr = ifelse(is.na(fdr), "", format.pval(fdr, digits = 3, eps = 1e-4))
    ) %>%
    select(variable, variable_name, variable_group, unit, n, method, estimate, p_value, fdr, significant_fdr)
}

normality_report_table <- normality_results %>%
  mutate(
    shapiro_p = ifelse(is.na(shapiro_p), "", format.pval(shapiro_p, digits = 3, eps = 1e-4))
  ) %>%
  select(sample_set, variable, variable_name, variable_group, unit, n, shapiro_p, is_normal)

distribution_report_table <- distribution_summary %>%
  mutate(percent = fmt_pct(percent)) %>%
  arrange(variable, sample_set, desc(n))

all_nominal <- result_summary_table("all_samples")
clean_200_nominal <- result_summary_table("clean_threshold_200")
clean_400_nominal <- result_summary_table("clean_threshold_400")

all_sig_n <- sum(filter(correlation_results, sample_set == "all_samples")$significant_fdr, na.rm = TRUE)
clean_200_sig_n <- sum(filter(correlation_results, sample_set == "clean_threshold_200")$significant_fdr, na.rm = TRUE)
clean_400_sig_n <- sum(filter(correlation_results, sample_set == "clean_threshold_400")$significant_fdr, na.rm = TRUE)

summary_lines <- c(
  paste0("- Correlation target: `", x_var, "`"),
  paste0("- Exploratory thresholds: `", threshold_200, "` and `", threshold_400, "` ", x_unit),
  paste0("- Threshold 200 definition: `>= ", threshold_200, "` contaminated, `< ", threshold_200, "` clean"),
  paste0("- Threshold 400 definition: `>= ", threshold_400, "` contaminated, `< ", threshold_400, "` clean"),
  "- Normality test: Shapiro-Wilk",
  "- Correlation rule: Pearson only when both ethanol and the lipoprofile variable passed normality in that round; otherwise Spearman",
  paste0("- FDR threshold for significance: ", fdr_alpha),
  paste0("- Sample counts: all = ", sample_counts$n[sample_counts$sample_set == "all_samples"],
         ", clean <200 = ", sample_counts$n[sample_counts$sample_set == "clean_threshold_200"],
         ", clean <400 = ", sample_counts$n[sample_counts$sample_set == "clean_threshold_400"]),
  paste0("- Nominal associations with p < 0.05: all samples = ", nrow(all_nominal$nominal),
         ", clean <200 = ", nrow(clean_200_nominal$nominal),
         ", clean <400 = ", nrow(clean_400_nominal$nominal)),
  paste0("- Significant associations after FDR: all samples = ", all_sig_n,
         ", clean <200 = ", clean_200_sig_n,
         ", clean <400 = ", clean_400_sig_n)
)

metadata_plot_md_lines <- c(
  "### All Samples",
  "",
  "![](plots/metadata_distribution_all_samples.png)",
  "",
  "### Clean Samples (< 200 mg/dL)",
  "",
  "![](plots/metadata_distribution_clean_threshold_200.png)",
  "",
  "### Clean Samples (< 400 mg/dL)",
  "",
  "![](plots/metadata_distribution_clean_threshold_400.png)",
  ""
)

plot_md_lines <- if (nrow(plot_manifest) == 0) {
  c("No nominal p < 0.05 correlations were detected, so no scatterplots were generated.")
} else {
  unlist(map(seq_len(nrow(plot_manifest)), function(i) {
    row <- plot_manifest[i, ]
    badge <- if (isTRUE(row$significant_fdr)) " [FDR-significant]" else ""
    c(
      paste0("### ", row$sample_set, " - ", row$variable, badge),
      "",
      paste0("![](plots/", row$plot_file, ")"),
      ""
    )
  }))
}

md_lines <- c(
  "# Ethanol Contamination Impact on NMR Lipoprofile",
  "",
  "## Summary",
  "",
  summary_lines,
  "",
  "## Inputs",
  "",
  "- `clean_data/ASPREE_ethanol_contamination_clean.csv`",
  "- `clean_data/ASPREE_merged_lipoprofile_clean.csv`",
  "- `raw_data/subject_metadata.xlsx`",
  "- `clean_data/ASPREE_merged_lipoprofile_codebook.xlsx`",
  "",
  "## Sample Counts",
  "",
  knitr::kable(sample_counts, format = "pipe"),
  "",
  "## Metadata Distribution Summary",
  "",
  metadata_plot_md_lines,
  "",
  knitr::kable(distribution_report_table, format = "pipe"),
  "",
  "## Normality Results",
  "",
  knitr::kable(normality_report_table, format = "pipe"),
  "",
  "## Correlation Summary: All Samples",
  "",
  paste0("Nominal associations with p < 0.05: ", nrow(all_nominal$nominal)),
  "",
  if (nrow(all_nominal$nominal) > 0) knitr::kable(all_nominal$nominal %>% mutate(estimate = fmt_num(estimate), p_value = format.pval(p_value, digits = 3, eps = 1e-4), fdr = format.pval(fdr, digits = 3, eps = 1e-4)), format = "pipe") else "No variables met nominal p < 0.05.",
  "",
  paste0("FDR-significant associations: ", nrow(all_nominal$fdr)),
  "",
  if (nrow(all_nominal$fdr) > 0) knitr::kable(all_nominal$fdr %>% mutate(estimate = fmt_num(estimate), p_value = format.pval(p_value, digits = 3, eps = 1e-4), fdr = format.pval(fdr, digits = 3, eps = 1e-4)), format = "pipe") else "No variables passed FDR correction.",
  "",
  "## Correlation Results: All Samples",
  "",
  knitr::kable(correlation_report_table("all_samples"), format = "pipe"),
  "",
  "## Correlation Summary: Clean Samples (< 200 mg/dL)",
  "",
  paste0("Nominal associations with p < 0.05: ", nrow(clean_200_nominal$nominal)),
  "",
  if (nrow(clean_200_nominal$nominal) > 0) knitr::kable(clean_200_nominal$nominal %>% mutate(estimate = fmt_num(estimate), p_value = format.pval(p_value, digits = 3, eps = 1e-4), fdr = format.pval(fdr, digits = 3, eps = 1e-4)), format = "pipe") else "No variables met nominal p < 0.05.",
  "",
  paste0("FDR-significant associations: ", nrow(clean_200_nominal$fdr)),
  "",
  if (nrow(clean_200_nominal$fdr) > 0) knitr::kable(clean_200_nominal$fdr %>% mutate(estimate = fmt_num(estimate), p_value = format.pval(p_value, digits = 3, eps = 1e-4), fdr = format.pval(fdr, digits = 3, eps = 1e-4)), format = "pipe") else "No variables passed FDR correction.",
  "",
  "## Correlation Results: Clean Samples (< 200 mg/dL)",
  "",
  knitr::kable(correlation_report_table("clean_threshold_200"), format = "pipe"),
  "",
  "## Correlation Summary: Clean Samples (< 400 mg/dL)",
  "",
  paste0("Nominal associations with p < 0.05: ", nrow(clean_400_nominal$nominal)),
  "",
  if (nrow(clean_400_nominal$nominal) > 0) knitr::kable(clean_400_nominal$nominal %>% mutate(estimate = fmt_num(estimate), p_value = format.pval(p_value, digits = 3, eps = 1e-4), fdr = format.pval(fdr, digits = 3, eps = 1e-4)), format = "pipe") else "No variables met nominal p < 0.05.",
  "",
  paste0("FDR-significant associations: ", nrow(clean_400_nominal$fdr)),
  "",
  if (nrow(clean_400_nominal$fdr) > 0) knitr::kable(clean_400_nominal$fdr %>% mutate(estimate = fmt_num(estimate), p_value = format.pval(p_value, digits = 3, eps = 1e-4), fdr = format.pval(fdr, digits = 3, eps = 1e-4)), format = "pipe") else "No variables passed FDR correction.",
  "",
  "## Correlation Results: Clean Samples (< 400 mg/dL)",
  "",
  knitr::kable(correlation_report_table("clean_threshold_400"), format = "pipe"),
  "",
  "## Significant Correlation Plots",
  "",
  plot_md_lines,
  "",
  "## Output Files",
  "",
  "- `normality_results.csv`",
  "- `correlation_results_all_samples.csv`",
  "- `correlation_results_clean_threshold_200.csv`",
  "- `correlation_results_clean_threshold_400.csv`",
  "- `metadata_distribution_summary.csv`",
  "- `metadata_plot_manifest.csv`",
  "- `sample_counts.csv`",
  "- `analysis_dataset_with_metadata.csv`",
  "- `significant_results_with_plots.csv` if any nominal p < 0.05 correlations are present",
  "- `plots/` directory with ggplot figures for all nominal p < 0.05 correlations; FDR-significant plots are specially labeled",
  ""
)

writeLines(md_lines, file.path(results_dir, "ethanol_nmr_impact_report.md"))

html_table <- function(df) {
  knitr::kable(df, format = "html", table.attr = "class='table table-striped table-sm'")
}

plot_html_blocks <- if (nrow(plot_manifest) == 0) {
  "<p>No nominal p &lt; 0.05 correlations were detected, so no scatterplots were generated.</p>"
} else {
  paste(
    map_chr(seq_len(nrow(plot_manifest)), function(i) {
      row <- plot_manifest[i, ]
      badge <- if (isTRUE(row$significant_fdr)) {
        "<div style='display:inline-block;background:#b22222;color:#fff;font-weight:bold;padding:4px 8px;border-radius:4px;margin-bottom:8px;'>FDR-significant in this sample set</div>"
      } else {
        "<div style='display:inline-block;background:#eaf2f8;color:#1f4e78;padding:4px 8px;border-radius:4px;margin-bottom:8px;'>Nominal p &lt; 0.05</div>"
      }
      paste0(
        "<h3>", row$sample_set, " - ", row$variable, ifelse(isTRUE(row$significant_fdr), " [FDR-significant]", ""), "</h3>",
        badge,
        "<img src='plots/", row$plot_file, "' alt='", row$variable, "' style='max-width:900px; width:100%; height:auto;'/>"
      )
    }),
    collapse = "\n"
  )
}

summary_html_table <- function(df) {
  if (nrow(df) == 0) {
    return("<p>None.</p>")
  }
  html_table(
    df %>%
      mutate(
        estimate = fmt_num(estimate),
        p_value = format.pval(p_value, digits = 3, eps = 1e-4),
        fdr = format.pval(fdr, digits = 3, eps = 1e-4)
      )
  )
}

html_body <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>Ethanol Contamination Impact on NMR Lipoprofile</title>",
  "<style>",
  "body{font-family:Arial,sans-serif;max-width:1400px;margin:40px auto;padding:0 24px;line-height:1.45;color:#222;}",
  "table{border-collapse:collapse;width:100%;margin:16px 0 28px 0;font-size:13px;}",
  "th,td{border:1px solid #ddd;padding:6px 8px;text-align:left;vertical-align:top;}",
  "th{background:#eaf2f8;position:sticky;top:0;}",
  "h1,h2,h3{color:#1f4e78;}",
  "code{background:#f4f4f4;padding:2px 4px;border-radius:3px;}",
  "ul{margin-top:0;}",
  "</style></head><body>",
  "<h1>Ethanol Contamination Impact on NMR Lipoprofile</h1>",
  "<h2>Summary</h2><ul><li>", paste(gsub("^- ", "", summary_lines), collapse = "</li><li>"), "</li></ul>",
  "<h2>Inputs</h2><ul>",
  "<li><code>clean_data/ASPREE_ethanol_contamination_clean.csv</code></li>",
  "<li><code>clean_data/ASPREE_merged_lipoprofile_clean.csv</code></li>",
  "<li><code>raw_data/subject_metadata.xlsx</code></li>",
  "<li><code>clean_data/ASPREE_merged_lipoprofile_codebook.xlsx</code></li>",
  "</ul>",
  "<h2>Sample Counts</h2>", html_table(sample_counts),
  "<h2>Metadata Distribution Summary</h2>",
  "<h3>All Samples</h3><img src='plots/metadata_distribution_all_samples.png' alt='metadata all samples' style='max-width:1100px; width:100%; height:auto;'/>",
  "<h3>Clean Samples (&lt; 200 mg/dL)</h3><img src='plots/metadata_distribution_clean_threshold_200.png' alt='metadata clean threshold 200' style='max-width:1100px; width:100%; height:auto;'/>",
  "<h3>Clean Samples (&lt; 400 mg/dL)</h3><img src='plots/metadata_distribution_clean_threshold_400.png' alt='metadata clean threshold 400' style='max-width:1100px; width:100%; height:auto;'/>",
  html_table(distribution_report_table),
  "<h2>Normality Results</h2>", html_table(normality_report_table),
  "<h2>Correlation Summary: All Samples</h2>",
  "<p>Nominal associations with p &lt; 0.05: ", nrow(all_nominal$nominal), "</p>",
  summary_html_table(all_nominal$nominal),
  "<p>FDR-significant associations: ", nrow(all_nominal$fdr), "</p>",
  summary_html_table(all_nominal$fdr),
  "<h2>Correlation Results: All Samples</h2>", html_table(correlation_report_table("all_samples")),
  "<h2>Correlation Summary: Clean Samples (&lt; 200 mg/dL)</h2>",
  "<p>Nominal associations with p &lt; 0.05: ", nrow(clean_200_nominal$nominal), "</p>",
  summary_html_table(clean_200_nominal$nominal),
  "<p>FDR-significant associations: ", nrow(clean_200_nominal$fdr), "</p>",
  summary_html_table(clean_200_nominal$fdr),
  "<h2>Correlation Results: Clean Samples (&lt; 200 mg/dL)</h2>", html_table(correlation_report_table("clean_threshold_200")),
  "<h2>Correlation Summary: Clean Samples (&lt; 400 mg/dL)</h2>",
  "<p>Nominal associations with p &lt; 0.05: ", nrow(clean_400_nominal$nominal), "</p>",
  summary_html_table(clean_400_nominal$nominal),
  "<p>FDR-significant associations: ", nrow(clean_400_nominal$fdr), "</p>",
  summary_html_table(clean_400_nominal$fdr),
  "<h2>Correlation Results: Clean Samples (&lt; 400 mg/dL)</h2>", html_table(correlation_report_table("clean_threshold_400")),
  "<h2>Significant Correlation Plots</h2>", plot_html_blocks,
  "</body></html>"
)

writeLines(html_body, file.path(results_dir, "ethanol_nmr_impact_report.html"))

message("Analysis complete.")
message("Results directory: ", results_dir)
