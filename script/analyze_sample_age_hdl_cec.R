options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(htmltools)
  library(knitr)
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
script_dir <- if (is.na(script_path)) {
  normalizePath(file.path(getwd(), "script"), winslash = "/", mustWork = FALSE)
} else {
  dirname(script_path)
}
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

clean_dir <- file.path(project_root, "clean_data")
outputs_dir <- file.path(project_root, "outputs", "cec_analysis")
raw_dir <- file.path(project_root, "raw_data")
results_dir <- file.path(project_root, "results", "sample_age_hdl_cec_correlation")
plots_dir <- file.path(results_dir, "plots")

if (dir.exists(results_dir)) {
  unlink(results_dir, recursive = TRUE, force = TRUE)
}
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "grey25", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "grey15"),
      axis.line = element_line(linewidth = 0.35, color = "grey15"),
      axis.ticks = element_line(linewidth = 0.3, color = "grey15"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", hjust = 0),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey90"),
      plot.margin = margin(7, 7, 7, 7)
    )
}

save_plot <- function(plot, filename, width, height) {
  png_path <- file.path(plots_dir, paste0(filename, ".png"))
  pdf_path <- file.path(plots_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 320, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  tibble(plot = filename, png = file.path("plots", paste0(filename, ".png")), pdf = file.path("plots", paste0(filename, ".pdf")))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 0.001 ~ "<0.001",
    TRUE ~ formatC(x, format = "f", digits = 3)
  )
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(na_if(as.character(x), "NA")))
}

run_correlation <- function(df, x_var, y_var, method, family, y_label) {
  sub <- df %>%
    select(all_of(c("subject_id", x_var, y_var))) %>%
    mutate(
      x = .data[[x_var]],
      y = .data[[y_var]]
    ) %>%
    filter(is.finite(x), is.finite(y))

  test <- if (nrow(sub) >= 3 && n_distinct(sub$x) >= 3 && n_distinct(sub$y) >= 3) {
    tryCatch(
      suppressWarnings(cor.test(sub$x, sub$y, method = method, exact = FALSE)),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  tibble(
    analysis_family = family,
    x_variable = x_var,
    y_variable = y_var,
    y_label = y_label,
    method = method,
    n = nrow(sub),
    estimate = if (is.null(test)) NA_real_ else unname(test$estimate[[1]]),
    p_value = if (is.null(test)) NA_real_ else test$p.value
  )
}

nmr_raw <- read_csv(file.path(clean_dir, "ASPREE_merged_lipoprofile_clean.csv"), show_col_types = FALSE)
cec_raw <- read_csv(file.path(outputs_dir, "cec_sample_summary.csv"), show_col_types = FALSE)
age_raw <- read_excel(file.path(raw_dir, "sample_age.xlsx"), sheet = "Extracted Data")
codebook_raw <- read_excel(file.path(clean_dir, "ASPREE_merged_lipoprofile_codebook.xlsx"), sheet = "Codebook")

nmr <- nmr_raw %>%
  rename(
    subject_id = `Subject ID`,
    labcorp_accession = `Labcorp Accession Number`,
    collection_date_nmr = `Collection Date`
  )

age <- age_raw %>%
  transmute(
    subject_id = as.character(subject_id),
    collection_year_age_file = safe_numeric(collection_year),
    sample_age = safe_numeric(sample_age)
  )

cec <- cec_raw %>%
  filter(str_detect(sample, "^ASP-[0-9]{3}$")) %>%
  transmute(
    subject_id = sample,
    cec_plate = plate,
    cec_n_wells = n_wells,
    cec_wells = wells,
    pct_cec_mean,
    pct_cec_sd,
    cec_index_plate_qc,
    cec_index_global_qc
  )

hdl_vars <- c(
  "cHDLP", "L-cHDLP", "M-cHDLP", "S-cHDLP",
  "H7P", "H6P", "H5P", "H4P", "H3P", "H2P", "H1P",
  "HDLZ", "MVX_MSCHDLP", "ApoA1"
)
requested_lipoprotein_vars <- c(
  "TRLP", "VL-TRLP", "L-TRLP", "M-TRLP", "S-TRLP", "VS-TRLP",
  "cLDLP", "L-cLDLP", "M-cLDLP", "S-cLDLP",
  "cHDLP", "L-cHDLP", "M-cHDLP", "S-cHDLP",
  "H7P", "H6P", "H5P", "H4P", "H3P", "H2P", "H1P",
  "TRLZ", "LDLZ", "HDLZ", "NTRLTG", "NTRLC"
)
cec_vars <- c("cec_index_plate_qc", "cec_index_global_qc")
context_cec_vars <- c("pct_cec_mean")

codebook <- codebook_raw %>%
  transmute(
    variable = abbreviated_variable_name,
    variable_name,
    variable_group = variable_group_description,
    unit = measurement_units
  )

all_nmr_vars <- codebook %>%
  filter(
    variable %in% names(nmr_raw),
    !variable %in% c("Labcorp Accession Number", "Subject ID", "Collection Date", "Comment")
  ) %>%
  pull(variable)

all_nmr_variable_manifest <- codebook %>%
  filter(variable %in% all_nmr_vars) %>%
  mutate(
    analysis_family = "All NMR variables",
    inclusion_note = "Included in the all-variable NMR screen because it is a non-identifier variable in the clean NMR/lipoprotein codebook."
  )

variable_manifest <- bind_rows(
  codebook %>%
    filter(variable %in% hdl_vars) %>%
    mutate(
      analysis_family = "NMR HDL characteristics",
      inclusion_note = case_when(
        variable == "ApoA1" ~ "Included as an HDL-related apolipoprotein reported in the NMR/lipoprotein file.",
        variable == "MVX_MSCHDLP" ~ "Included because the codebook defines it as Medium/Small HDL Particles.",
        TRUE ~ "Included because the codebook label directly refers to HDL particles, HDL subclasses, or HDL size."
      )
    ),
  codebook %>%
    filter(variable %in% requested_lipoprotein_vars) %>%
    mutate(
      analysis_family = "Lipoprotein characteristics",
      inclusion_note = "Included in the TRL/LDL/HDL lipoprotein characteristic panel."
    ),
  tibble(
    variable = c(cec_vars, context_cec_vars),
    variable_name = c("CEC index normalized to plate QC", "CEC index normalized to global QC", "NC-corrected percent CEC mean"),
    variable_group = c("CEC", "CEC", "CEC"),
    unit = c("index", "index", "%"),
    analysis_family = c("CEC index", "CEC index", "CEC context"),
    inclusion_note = c(
      "Primary CEC index requested by the user; plate-specific QC normalization.",
      "Primary CEC index requested by the user; global QC normalization.",
      "Supplemental CEC measure retained for interpretation but not the primary CEC-index family."
    )
  )
) %>%
  arrange(analysis_family, variable)

excluded_hdl_like <- codebook %>%
  filter(variable == "NHDLC") %>%
  mutate(
    analysis_family = "Excluded",
    inclusion_note = "Excluded from the HDL-characteristics family because the codebook defines it as non-HDL cholesterol."
  )

analysis_df <- nmr %>%
  mutate(across(all_of(all_nmr_vars), safe_numeric)) %>%
  left_join(age, by = "subject_id") %>%
  left_join(cec, by = "subject_id")

sample_summary <- tibble(
  metric = c(
    "NMR samples",
    "NMR samples with sample_age",
    "NMR samples with CEC index",
    "NMR samples with both sample_age and CEC index",
    "CEC ASP sample rows used",
    "CEC non-ASP/control rows excluded"
  ),
  n = c(
    nrow(nmr),
    sum(is.finite(analysis_df$sample_age)),
    sum(is.finite(analysis_df$cec_index_plate_qc) | is.finite(analysis_df$cec_index_global_qc)),
    sum(is.finite(analysis_df$sample_age) & (is.finite(analysis_df$cec_index_plate_qc) | is.finite(analysis_df$cec_index_global_qc))),
    nrow(cec),
    nrow(cec_raw) - nrow(cec)
  )
)

sample_age_distribution <- analysis_df %>%
  count(sample_age, name = "n") %>%
  arrange(sample_age)

correlations <- bind_rows(
  map_dfr(hdl_vars, function(v) {
    label <- variable_manifest$variable_name[match(v, variable_manifest$variable)]
    bind_rows(
      run_correlation(analysis_df, "sample_age", v, "spearman", "NMR HDL characteristics", label),
      run_correlation(analysis_df, "sample_age", v, "pearson", "NMR HDL characteristics", label)
    )
  }),
  map_dfr(cec_vars, function(v) {
    label <- variable_manifest$variable_name[match(v, variable_manifest$variable)]
    bind_rows(
      run_correlation(analysis_df, "sample_age", v, "spearman", "CEC index", label),
      run_correlation(analysis_df, "sample_age", v, "pearson", "CEC index", label)
    )
  }),
  map_dfr(context_cec_vars, function(v) {
    label <- variable_manifest$variable_name[match(v, variable_manifest$variable)]
    bind_rows(
      run_correlation(analysis_df, "sample_age", v, "spearman", "CEC context", label),
      run_correlation(analysis_df, "sample_age", v, "pearson", "CEC context", label)
    )
  })
) %>%
  group_by(analysis_family, method) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    estimate_label = fmt_num(estimate),
    p_label = fmt_p(p_value),
    fdr_label = fmt_p(fdr),
    significant_fdr_0_05 = !is.na(fdr) & fdr < 0.05
  )

hdl_correlations <- correlations %>% filter(analysis_family == "NMR HDL characteristics")
cec_correlations <- correlations %>% filter(analysis_family == "CEC index")
requested_lipoprotein_correlations <- bind_rows(
  map_dfr(requested_lipoprotein_vars, function(v) {
    label <- variable_manifest$variable_name[
      match(v, variable_manifest$variable)
    ]
    bind_rows(
      run_correlation(analysis_df, "sample_age", v, "spearman", "Lipoprotein characteristics", label),
      run_correlation(analysis_df, "sample_age", v, "pearson", "Lipoprotein characteristics", label)
    )
  })
) %>%
  group_by(analysis_family, method) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    estimate_label = fmt_num(estimate),
    p_label = fmt_p(p_value),
    fdr_label = fmt_p(fdr),
    significant_fdr_0_05 = !is.na(fdr) & fdr < 0.05
  )
correlations <- bind_rows(correlations, requested_lipoprotein_correlations)
hdl_correlations <- correlations %>% filter(analysis_family == "NMR HDL characteristics")
cec_correlations <- correlations %>% filter(analysis_family == "CEC index")
requested_lipoprotein_correlations <- correlations %>% filter(analysis_family == "Lipoprotein characteristics")

all_nmr_correlations <- bind_rows(
  map_dfr(all_nmr_vars, function(v) {
    label <- codebook$variable_name[match(v, codebook$variable)]
    bind_rows(
      run_correlation(analysis_df, "sample_age", v, "spearman", "All NMR variables", label),
      run_correlation(analysis_df, "sample_age", v, "pearson", "All NMR variables", label)
    )
  })
) %>%
  group_by(analysis_family, method) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    estimate_label = fmt_num(estimate),
    p_label = fmt_p(p_value),
    fdr_label = fmt_p(fdr),
    significant_fdr_0_05 = !is.na(fdr) & fdr < 0.05
  )

write_csv(analysis_df, file.path(results_dir, "analysis_dataset_sample_age_hdl_cec.csv"))
write_csv(sample_summary, file.path(results_dir, "sample_inclusion_summary.csv"))
write_csv(sample_age_distribution, file.path(results_dir, "sample_age_distribution.csv"))
write_csv(variable_manifest, file.path(results_dir, "variable_manifest.csv"))
write_csv(all_nmr_variable_manifest, file.path(results_dir, "all_nmr_variable_manifest.csv"))
write_csv(excluded_hdl_like, file.path(results_dir, "excluded_hdl_like_variables.csv"))
write_csv(all_nmr_correlations, file.path(results_dir, "correlation_results_all_nmr_variables.csv"))
write_csv(bind_rows(all_nmr_correlations, cec_correlations), file.path(results_dir, "correlation_results_all.csv"))
write_csv(hdl_correlations, file.path(results_dir, "correlation_results_sample_age_vs_hdl.csv"))
write_csv(cec_correlations, file.path(results_dir, "correlation_results_sample_age_vs_cec_index.csv"))
write_csv(requested_lipoprotein_correlations, file.path(results_dir, "correlation_results_sample_age_vs_lipoprotein_panel.csv"))

plot_manifest <- list()

p_age <- analysis_df %>%
  filter(is.finite(sample_age)) %>%
  ggplot(aes(x = factor(sample_age))) +
  geom_bar(fill = "#4C78A8", width = 0.72) +
  labs(
    title = "Sample age distribution",
    x = "Sample age (years)",
    y = "Number of NMR samples"
  ) +
  theme_pub()
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_age, "sample_age_distribution", 4.8, 3.4)

hdl_spearman <- hdl_correlations %>%
  filter(method == "spearman") %>%
  mutate(
    y_label = factor(y_label, levels = y_label[order(estimate)]),
    sig = if_else(significant_fdr_0_05, "FDR < 0.05", "Not FDR-significant")
  )

p_lollipop <- ggplot(hdl_spearman, aes(x = estimate, y = y_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_segment(aes(x = 0, xend = estimate, yend = y_label), linewidth = 0.45, color = "grey55") +
  geom_point(aes(fill = sig, size = -log10(p_value)), shape = 21, color = "grey10", stroke = 0.25) +
  scale_fill_manual(values = c("FDR < 0.05" = "#B24745", "Not FDR-significant" = "#4C78A8")) +
  scale_size_continuous(range = c(1.8, 5), name = expression(-log[10](italic(P)))) +
  labs(
    title = "Sample age versus NMR HDL characteristics",
    subtitle = "Spearman correlations; P values adjusted across HDL variables by Benjamini-Hochberg FDR",
    x = "Spearman rho",
    y = NULL,
    fill = NULL
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_pub(9)
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_lollipop, "sample_age_vs_hdl_spearman_lollipop", 7.0, 5.1)

p_heat <- hdl_spearman %>%
  mutate(label = glue("{fmt_num(estimate, 2)}\nFDR {fmt_p(fdr)}")) %>%
  ggplot(aes(x = "Sample age", y = y_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.6, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#B24745",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(
    title = "Correlation heatmap",
    x = NULL,
    y = NULL
  ) +
  theme_pub(9) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_heat, "sample_age_vs_hdl_spearman_heatmap", 5.7, 5.1)

hdl_long <- analysis_df %>%
  select(subject_id, sample_age, all_of(hdl_vars)) %>%
  pivot_longer(all_of(hdl_vars), names_to = "variable", values_to = "value") %>%
  left_join(codebook %>% select(variable, variable_name, unit), by = "variable") %>%
  filter(is.finite(sample_age), is.finite(value)) %>%
  mutate(variable_name = factor(variable_name, levels = hdl_spearman$y_label[order(as.character(hdl_spearman$y_label))]))

p_hdl_scatter <- ggplot(hdl_long, aes(x = sample_age, y = value)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), alpha = 0.58, size = 1.35, color = "#2F4B5B") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, color = "#B24745", fill = "#B24745", alpha = 0.15) +
  facet_wrap(~ variable_name, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = sort(unique(analysis_df$sample_age[is.finite(analysis_df$sample_age)]))) +
  labs(
    title = "Sample age and NMR HDL measurements",
    x = "Sample age (years)",
    y = "NMR measurement"
  ) +
  theme_pub(8) +
  theme(panel.grid.major.x = element_line(linewidth = 0.18, color = "grey92"))
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_hdl_scatter, "sample_age_vs_hdl_scatter_facets", 9.4, 7.0)

requested_spearman <- requested_lipoprotein_correlations %>%
  filter(method == "spearman") %>%
  mutate(
    y_label = factor(y_label, levels = y_label[order(estimate)]),
    sig = if_else(significant_fdr_0_05, "FDR < 0.05", "Not FDR-significant")
  )

p_requested_lollipop <- ggplot(requested_spearman, aes(x = estimate, y = y_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_segment(aes(x = 0, xend = estimate, yend = y_label), linewidth = 0.45, color = "grey55") +
  geom_point(aes(fill = sig, size = -log10(p_value)), shape = 21, color = "grey10", stroke = 0.25) +
  scale_fill_manual(values = c("FDR < 0.05" = "#B24745", "Not FDR-significant" = "#4C78A8")) +
  scale_size_continuous(range = c(1.8, 5), name = expression(-log[10](italic(P)))) +
  labs(
    title = "Sample age versus lipoprotein panel",
    subtitle = "Spearman correlations; P values adjusted across lipoprotein-panel variables by Benjamini-Hochberg FDR",
    x = "Spearman rho",
    y = NULL,
    fill = NULL
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_pub(9)
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_requested_lollipop, "sample_age_vs_lipoprotein_panel_spearman_lollipop", 7.4, 7.3)

p_requested_heat <- requested_spearman %>%
  mutate(label = glue("{fmt_num(estimate, 2)}\nFDR {fmt_p(fdr)}")) %>%
  ggplot(aes(x = "Sample age", y = y_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.4, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "white", high = "#B24745",
    midpoint = 0, limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(
    title = "Lipoprotein panel heatmap",
    x = NULL,
    y = NULL
  ) +
  theme_pub(8) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_requested_heat, "sample_age_vs_lipoprotein_panel_spearman_heatmap", 5.9, 7.3)

requested_long <- analysis_df %>%
  select(subject_id, sample_age, all_of(requested_lipoprotein_vars)) %>%
  pivot_longer(all_of(requested_lipoprotein_vars), names_to = "variable", values_to = "value") %>%
  left_join(codebook %>% select(variable, variable_name, unit), by = "variable") %>%
  filter(is.finite(sample_age), is.finite(value)) %>%
  mutate(variable_name = factor(variable_name, levels = requested_spearman$y_label[order(as.character(requested_spearman$y_label))]))

p_requested_scatter <- ggplot(requested_long, aes(x = sample_age, y = value)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), alpha = 0.58, size = 1.2, color = "#2F4B5B") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.48, color = "#B24745", fill = "#B24745", alpha = 0.15) +
  facet_wrap(~ variable_name, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = sort(unique(analysis_df$sample_age[is.finite(analysis_df$sample_age)]))) +
  labs(
    title = "Sample age and lipoprotein panel",
    x = "Sample age (years)",
    y = "NMR measurement"
  ) +
  theme_pub(7.5) +
  theme(panel.grid.major.x = element_line(linewidth = 0.18, color = "grey92"))
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_requested_scatter, "sample_age_vs_lipoprotein_panel_scatter_facets", 10.0, 9.4)

all_nmr_spearman <- all_nmr_correlations %>%
  filter(method == "spearman") %>%
  mutate(
    y_display = paste0(y_label, " (", y_variable, ")"),
    y_display = factor(y_display, levels = y_display[order(estimate)]),
    sig = if_else(significant_fdr_0_05, "FDR < 0.05", "Not FDR-significant")
  )

sig_all_nmr <- all_nmr_spearman %>%
  filter(significant_fdr_0_05) %>%
  arrange(fdr, p_value)

if (nrow(sig_all_nmr) > 0) {
  p_all_sig <- ggplot(sig_all_nmr, aes(x = estimate, y = reorder(y_display, estimate))) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
    geom_segment(aes(x = 0, xend = estimate, yend = y_display), linewidth = 0.45, color = "grey55") +
    geom_point(aes(size = -log10(p_value)), shape = 21, fill = "#B24745", color = "grey10", stroke = 0.25) +
    scale_size_continuous(range = c(2.2, 5.5), name = expression(-log[10](italic(P)))) +
    labs(
      title = "FDR-significant associations in the all-NMR screen",
      subtitle = "Spearman correlations; P values adjusted across all NMR variables by Benjamini-Hochberg FDR",
      x = "Spearman rho",
      y = NULL
    ) +
    coord_cartesian(xlim = c(-1, 1)) +
    theme_pub(9)
  plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_all_sig, "sample_age_vs_all_nmr_fdr_significant_spearman_lollipop", 7.0, max(2.6, 1.4 + 0.32 * nrow(sig_all_nmr)))
}

cec_spearman <- cec_correlations %>%
  filter(method == "spearman") %>%
  mutate(
    y_label = factor(y_label, levels = y_label[order(estimate)]),
    sig = if_else(significant_fdr_0_05, "FDR < 0.05", "Not FDR-significant")
  )

p_cec_lollipop <- ggplot(cec_spearman, aes(x = estimate, y = y_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_segment(aes(x = 0, xend = estimate, yend = y_label), linewidth = 0.55, color = "grey55") +
  geom_point(aes(fill = sig, size = -log10(p_value)), shape = 21, color = "grey10", stroke = 0.25) +
  scale_fill_manual(values = c("FDR < 0.05" = "#B24745", "Not FDR-significant" = "#4C78A8")) +
  scale_size_continuous(range = c(2.5, 5.5), name = expression(-log[10](italic(P)))) +
  labs(
    title = "Sample age versus CEC index",
    subtitle = "Spearman correlations; P values adjusted across the two CEC index variables",
    x = "Spearman rho",
    y = NULL,
    fill = NULL
  ) +
  coord_cartesian(xlim = c(-1, 1)) +
  theme_pub(9)
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_cec_lollipop, "sample_age_vs_cec_index_spearman_lollipop", 6.0, 2.8)

cec_long <- analysis_df %>%
  select(subject_id, sample_age, all_of(cec_vars)) %>%
  pivot_longer(all_of(cec_vars), names_to = "variable", values_to = "value") %>%
  left_join(variable_manifest %>% select(variable, variable_name), by = "variable") %>%
  filter(is.finite(sample_age), is.finite(value))

p_cec_scatter <- ggplot(cec_long, aes(x = sample_age, y = value)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), alpha = 0.65, size = 1.7, color = "#2F4B5B") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.65, color = "#B24745", fill = "#B24745", alpha = 0.15) +
  facet_wrap(~ variable_name, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = sort(unique(analysis_df$sample_age[is.finite(analysis_df$sample_age)]))) +
  labs(
    title = "Sample age and CEC index",
    x = "Sample age (years)",
    y = "CEC index"
  ) +
  theme_pub(9) +
  theme(panel.grid.major.x = element_line(linewidth = 0.18, color = "grey92"))
plot_manifest[[length(plot_manifest) + 1]] <- save_plot(p_cec_scatter, "sample_age_vs_cec_index_scatter_facets", 6.2, 3.4)

plot_manifest_df <- bind_rows(plot_manifest)
write_csv(plot_manifest_df, file.path(results_dir, "plot_manifest.csv"))

primary_hdl_table <- hdl_spearman %>%
  arrange(fdr, p_value) %>%
  transmute(
    Variable = y_label,
    N = n,
    `Spearman rho` = estimate_label,
    `P value` = p_label,
    FDR = fdr_label
  )

primary_cec_table <- cec_spearman %>%
  arrange(fdr, p_value) %>%
  transmute(
    Variable = y_label,
    N = n,
    `Spearman rho` = estimate_label,
    `P value` = p_label,
    FDR = fdr_label
  )

primary_requested_table <- requested_spearman %>%
  arrange(fdr, p_value) %>%
  transmute(
    Variable = y_label,
    N = n,
    `Spearman rho` = estimate_label,
    `P value` = p_label,
    FDR = fdr_label
  )

primary_all_fdr_table <- sig_all_nmr %>%
  transmute(
    Variable = y_display,
    N = n,
    `Spearman rho` = estimate_label,
    `P value` = p_label,
    FDR = fdr_label
  )

sig_hdl <- hdl_spearman %>% filter(significant_fdr_0_05)
sig_cec <- cec_spearman %>% filter(significant_fdr_0_05)
sig_requested <- requested_spearman %>% filter(significant_fdr_0_05)

html_table <- function(df) {
  HTML(kable(df, format = "html", escape = TRUE, table.attr = "class='data-table'"))
}

img_tag <- function(path, alt) {
  tags$figure(
    tags$img(src = path, alt = alt),
    tags$figcaption(alt)
  )
}

report <- tagList(
  tags$html(
    tags$head(
      tags$meta(charset = "utf-8"),
      tags$title("Sample Age, NMR, HDL, and CEC Correlation Analysis"),
      tags$style(HTML("
        body { font-family: -apple-system, BlinkMacSystemFont, 'Helvetica Neue', Arial, sans-serif; color: #1f2933; margin: 0; background: #ffffff; }
        main { max-width: 1120px; margin: 0 auto; padding: 38px 34px 60px; }
        h1 { font-size: 28px; margin: 0 0 8px; }
        h2 { font-size: 19px; margin: 34px 0 12px; border-bottom: 1px solid #d8dee4; padding-bottom: 6px; }
        h3 { font-size: 15px; margin: 24px 0 8px; }
        p, li { font-size: 14px; line-height: 1.55; }
        .note { color: #4b5563; }
        .callout { border-left: 4px solid #4C78A8; padding: 10px 14px; background: #f7fafc; margin: 16px 0; }
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(430px, 1fr)); gap: 18px; align-items: start; }
        figure { margin: 14px 0 20px; }
        img { max-width: 100%; height: auto; border: 1px solid #e5e7eb; }
        figcaption { font-size: 12px; color: #4b5563; margin-top: 6px; }
        table.data-table { border-collapse: collapse; width: 100%; font-size: 13px; margin: 10px 0 18px; }
        table.data-table th, table.data-table td { border-bottom: 1px solid #e5e7eb; padding: 7px 8px; text-align: left; vertical-align: top; }
        table.data-table th { background: #f3f4f6; font-weight: 700; }
        code { background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }
      "))
    ),
    tags$body(
      tags$main(
        tags$h1("Sample Age, NMR, HDL, and CEC Correlation Analysis"),
        tags$p(class = "note", glue("Generated on {format(Sys.time(), '%Y-%m-%d %H:%M %Z')}.")),
        tags$div(
          class = "callout",
          tags$p("Primary analysis: Spearman correlation between sample age and all NMR variables, HDL-related NMR characteristics, and CEC index. Pearson correlations are included in the CSV outputs as supplemental analyses."),
          tags$p("Multiple testing correction: Benjamini-Hochberg FDR within each analysis family and method.")
        ),
        tags$h2("Input Data and Inclusion"),
        html_table(sample_summary),
        tags$p("Sample age was read from ", tags$code("raw_data/sample_age.xlsx"), ". NMR measurements were read from ", tags$code("clean_data/ASPREE_merged_lipoprofile_clean.csv"), ". CEC indices were read from ", tags$code("outputs/cec_analysis/cec_sample_summary.csv"), ". CEC rows were restricted to ASP sample IDs matching ", tags$code("^ASP-[0-9]{3}$"), ", excluding controls such as PC/QC/NC."),
        img_tag("plots/sample_age_distribution.png", "Distribution of sample age among NMR samples with non-missing sample age."),
        tags$h2("Variable Definition"),
        tags$p("The all-variable screen includes every non-identifier NMR variable listed in the clean NMR codebook. The HDL family shown below includes variables whose codebook labels directly describe HDL particle concentration, HDL subclasses, HDL size, medium/small HDL particles, plus ApoA1 as an HDL-related apolipoprotein. The non-HDL cholesterol variable was explicitly excluded from the HDL-focused family."),
        html_table(variable_manifest %>% filter(analysis_family == "NMR HDL characteristics") %>% select(`Analysis family` = analysis_family, Variable = variable, Label = variable_name, Unit = unit, Note = inclusion_note)),
        tags$h2("FDR-Significant Results from All NMR Variables"),
        tags$p("This section reports only variables significant at FDR < 0.05 after screening all non-identifier NMR variables with Spearman correlation."),
        if (nrow(sig_all_nmr) > 0) {
          tagList(
            img_tag("plots/sample_age_vs_all_nmr_fdr_significant_spearman_lollipop.png", "FDR-significant Spearman correlations from the all-NMR variable screen."),
            html_table(primary_all_fdr_table)
          )
        } else {
          tags$p("No NMR variable was significant at FDR < 0.05 in the all-variable Spearman screen.")
        },
        tags$h2("Sample Age versus NMR HDL Characteristics"),
        tags$div(class = "grid",
          img_tag("plots/sample_age_vs_hdl_spearman_lollipop.png", "Ranked Spearman correlations for sample age versus NMR HDL characteristics."),
          img_tag("plots/sample_age_vs_hdl_spearman_heatmap.png", "Heatmap of Spearman rho and FDR values for sample age versus NMR HDL characteristics.")
        ),
        img_tag("plots/sample_age_vs_hdl_scatter_facets.png", "Facet scatterplots for sample age versus each NMR HDL measurement. Points are horizontally jittered for visibility; red line is a linear fit shown for visual trend only."),
        tags$h3("Primary HDL Correlation Table"),
        html_table(primary_hdl_table),
        tags$p(
          if (nrow(sig_hdl) == 0) {
            "No NMR HDL characteristic was significant at FDR < 0.05 in the primary Spearman analysis."
          } else {
            glue("{nrow(sig_hdl)} NMR HDL characteristic(s) were significant at FDR < 0.05 in the primary Spearman analysis.")
          }
        ),
        tags$h2("Sample Age versus CEC Index"),
        tags$div(class = "grid",
          img_tag("plots/sample_age_vs_cec_index_spearman_lollipop.png", "Spearman correlations for sample age versus CEC indices."),
          img_tag("plots/sample_age_vs_cec_index_scatter_facets.png", "Facet scatterplots for sample age versus CEC indices.")
        ),
        tags$h3("Primary CEC Correlation Table"),
        html_table(primary_cec_table),
        tags$p(
          if (nrow(sig_cec) == 0) {
            "No CEC index was significant at FDR < 0.05 in the primary Spearman analysis."
          } else {
            glue("{nrow(sig_cec)} CEC index variable(s) were significant at FDR < 0.05 in the primary Spearman analysis.")
          }
        ),
        tags$h2("Output Files"),
        tags$ul(
          tags$li(tags$code("analysis_dataset_sample_age_hdl_cec.csv"), ": merged per-sample analysis dataset."),
          tags$li(tags$code("correlation_results_all_nmr_variables.csv"), ": Spearman and Pearson results for all non-identifier NMR variables."),
          tags$li(tags$code("correlation_results_all.csv"), ": all-NMR-variable results plus CEC-index results."),
          tags$li(tags$code("correlation_results_sample_age_vs_hdl.csv"), ": HDL-focused correlation results."),
          tags$li(tags$code("correlation_results_sample_age_vs_lipoprotein_panel.csv"), ": correlation results for the TRL/LDL/HDL lipoprotein panel."),
          tags$li(tags$code("correlation_results_sample_age_vs_cec_index.csv"), ": CEC-index-focused correlation results."),
          tags$li(tags$code("all_nmr_variable_manifest.csv"), ": variables included in the all-NMR screen."),
          tags$li(tags$code("variable_manifest.csv"), ": focused variable families and inclusion notes."),
          tags$li(tags$code("plots/"), ": PNG and PDF versions of all figures.")
        )
      )
    )
  )
)

save_html(report, file = file.path(results_dir, "sample_age_hdl_cec_correlation_report.html"))

message("Analysis complete: ", results_dir)
