# ---- GATA2 – Single & Combination Phenotype Analysis (clean) -----

# ---------- Packages ----------
req <- c("dplyr", "tidyr","data.table","tibble", "stringr","ggplot2","ggVennDiagram","eulerr", "readxl", "ggplot2","purrr","scales")
        
to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org", dependencies = TRUE)
lapply(req, library, character.only = TRUE)

# ---------- Paths & constants ----------
path_plots <- "/Plots/"
if (!dir.exists(path_plots)) dir.create(path_plots, recursive = TRUE, showWarnings = FALSE)

# ------- Counts case & control ----------
N_cases    <- 212
N_controls <- 288600 #296337#After removing GATA2 variants from the cohort##After removing participants with PLP Variants in 227 genes Final Number=296466 #Number=332193 #Removed Withdrawn Participants  # After removing VUS 236,067 Including VUS 330757

# -------- Input files ----------

phenotype_gata2_adult_path <- "phenotype_eid_df_Adults_v8.csv"
phenotype_ukb_path         <- "Phenotype_EID_Wide_v8.xlsx"
select_pheno_full_path     <- "GATA2_Phenotypes_Final_v8.xlsx"

comb_paths <- list(
  Combination_2 = "Combinationsof2_GATA2_022526.csv",
  Combination_3 = "Combinationsof3_GATA2_022526.csv",
  Combination_4 = "Combinationsof4_GATA2_022526.csv"
)

# Single-phenotype counts CSV (Phenotype_Index, Count_GATA2, Count_UKB)
single_pheno_csv_path <- "SinglePheno_GATA2_v8.csv"

# ---------- Helper Functions (stats) ----------
# Vectorized chi-square p for 2x2 (df=1)
chisq_p_2x2 <- function(a,b,c,d) {
  N   <- a + b + c + d
  num <- (a*d - b*c)^2 * N
  den <- (a + b) * (c + d) * (a + c) * (b + d)
  stat <- num / den
  pchisq(stat, df = 1, lower.tail = FALSE)
}

# Compute vectorized RR/OR + CI and p-values from 2x2 table defined by (a, b, c, d)
# a = cases with phenotype, b = cases without, c = controls with, d = controls without
# a =Counts_GATA2,b = N_cases-Counts_GATA2,c=Counts_UKB,d=N_Control-Counts_UKB
compute_2x2_metrics <- function(a, b, c, d, fisher_for_small = TRUE, show_progress = TRUE, pb_step = 1000) {
  # Haldane–Anscombe continuity correction where needed
  any_zero <- (a == 0L) | (b == 0L) | (c == 0L) | (d == 0L)
  aa <- ifelse(any_zero, a + 0.5, a)
  bb <- ifelse(any_zero, b + 0.5, b)
  cc <- ifelse(any_zero, c + 0.5, c)
  dd <- ifelse(any_zero, d + 0.5, d)
  
  # Risks
  p1 <- aa / (aa + bb)  # cases risk (Ratio of people with Phenotype in Cases)
  p0 <- cc / (cc + dd)  # controls risk (Ratio of people with Phenotypes in Control)
  
  # Risk Ratio (Wald CI)
  RR <- p1 / p0
  se_logRR <- sqrt(pmax(0, 1/aa - 1/(aa+bb) + 1/cc - 1/(cc+dd)))
  RR_L <- exp(log(RR) - 1.96 * se_logRR)
  RR_U <- exp(log(RR) + 1.96 * se_logRR)
  logRR<-log10(RR)
  # Odds Ratio (Wald CI)
  OR <- (aa * dd) / (bb * cc)   #(Ratio of Counts_GATA2/N_Cases-Counts_GATA2 )/Ratio of Counts_UKB/N_Controls-Counts_UKB)
  se_logOR <- sqrt(1/aa + 1/bb + 1/cc + 1/dd)
  OR_L <- exp(log(OR) - 1.96 * se_logOR)
  OR_U <- exp(log(OR) + 1.96 * se_logOR)
  logOR<-log10(OR)
  
  # p-values: chi-square by default
  p_value <- chisq_p_2x2(a, b, c, d)
  
  # Fisher only where any cell < 5
  if (isTRUE(fisher_for_small)) {
    small <- (a < 5L) | (b < 5L) | (c < 5L) | (d < 5L)
    if (any(small)) {
      idx <- which(small)
      if (isTRUE(show_progress)) {
        message(sprintf("Running Fisher's exact test on %d rows with small cells...", length(idx)))
        pb <- utils::txtProgressBar(min = 0, max = length(idx), style = 3)
      }
      p_fisher <- numeric(length(idx))
      for (k in seq_along(idx)) {
        i <- idx[k]
        p_fisher[k] <- fisher.test(matrix(c(a[i], b[i], c[i], d[i]), nrow = 2, byrow = TRUE))$p.value
        if (isTRUE(show_progress) && (k %% pb_step == 0L || k == length(idx))) {
          utils::setTxtProgressBar(pb, k)
        }
      }
      if (isTRUE(show_progress)) close(pb)
      p_value[idx] <- p_fisher
    }
  }
  
  # Classification metrics
  Sensitivity <- ifelse((a + b) > 0, a / (a + b), NA_real_)   #TP/TP+FN
  Specificity <- ifelse((c + d) > 0, d / (c + d), NA_real_)   #TN/TN+FP
  PPV         <- ifelse((a + c) > 0, a / (a + c), NA_real_)
  NPV         <- ifelse((b + d) > 0, d / (b + d), NA_real_)
  Accuracy    <- (a + d) / (a + b + c + d)
  PPL         <- ifelse((1 - Specificity) > 0, Sensitivity / (1 - Specificity), Inf)
  NLP         <- ifelse(Specificity > 0, (1 - Sensitivity) / Specificity, Inf)
  
  tibble::tibble(
    RR = RR,logRR=logRR,RR_Lower_CI = RR_L, RR_Upper_CI = RR_U,
    OR = OR,logOR=logOR,OR_Lower_CI = OR_L, OR_Upper_CI = OR_U,
    p_value = p_value,
    Sensitivity = Sensitivity, Specificity = Specificity, PPV = PPV, NPV = NPV,
    Accuracy = Accuracy, PPL = PPL, NLP = NLP
  )
}

# Signed log10 transform for differences
s_log10 <- function(x) sign(x) * log10(1 + abs(x))

# ---------- Generic processor usable for both singles and combinations ----------
process_any_fast <- function(df, id_col, case_col, ctrl_col,
                             N_cases, N_controls,
                             fisher_for_small = TRUE, show_progress = TRUE, pb_step = 1000) {
  # Coerce numeric counts (NA -> 0) and build 2x2 margins
  a <- suppressWarnings(as.numeric(df[[case_col]])); a[is.na(a)] <- 0 #Case #TP
  c <- suppressWarnings(as.numeric(df[[ctrl_col]])); c[is.na(c)] <- 0 #Control #FP
  b <- N_cases - a #FN
  d <- N_controls - c #TN
  
  # Ratios and transforms
  UKB_Ratio   <- c / N_controls
  GATA2_Ratio <- a / N_cases
  Ratio_Difference     <- UKB_Ratio - GATA2_Ratio
  Ratio_of_Proportions <- GATA2_Ratio / pmax(UKB_Ratio, .Machine$double.eps)
  LogRatio_Difference  <- s_log10(Ratio_Difference)
  
  # Vectorized metrics (RR/OR/CI/p)
  m <- compute_2x2_metrics(a, b, c, d, fisher_for_small, show_progress, pb_step)
  
  # Assemble output
  out <- df
  out[["UKB_Ratio"]]            <- UKB_Ratio
  out[["GATA2_Ratio"]]          <- GATA2_Ratio
  out[["Ratio_Difference"]]     <- Ratio_Difference
  out[["Ratio_of_Proportions"]] <- Ratio_of_Proportions
  out[["LogRatio_Difference"]]  <- LogRatio_Difference
  out <- cbind(out, as.data.frame(m))
  
  # Drop trivial and zero-case rows (as before)
  keep <- !((c == 0 & a == 0)) & (a > 0)
  out <- out[keep, , drop = FALSE]
  
  # Multiple-testing correction
  out[["padj"]] <- p.adjust(out$p_value, method = "bonferroni")
  out
}

# ---------- Load data ----------
phenotype_gata2_adult <- read.csv(phenotype_gata2_adult_path, stringsAsFactors = FALSE, check.names = FALSE)
phenotype_ukb         <- readxl::read_xlsx(phenotype_ukb_path)
phenotype_ukb[phenotype_ukb == 0] <- NA

select_pheno_full <- readxl::read_excel(select_pheno_full_path)
select_pheno <- select_pheno_full %>%
  filter(Keep_Adults == "Yes", Present_in_GATA2 == "Yes", Present_in_UKB == "Yes") %>%
  dplyr::select(Phenotype_Index, Phenotype, Phenotype_Description, Category)

select_pheno_combo<-select_pheno_full %>%
  filter(Keep_Adults == "Yes", Present_in_GATA2 == "Yes", Present_in_UKB == "Yes") %>%
  dplyr::select(Phenotype_Index, Phenotype, Phenotype_Description, Category)



# ---------- SINGLE PHENOTYPE SUMMARY ----------
# Read precomputed single-phenotype counts from CSV
# Expect columns: Phenotype_Index, Count_GATA2 (or Counts_GATA2), Count_UKB (or Counts_UKB)
counts <- read.csv(single_pheno_csv_path, stringsAsFactors = FALSE, check.names = FALSE)


single_pheno <- process_any_fast(
  counts,
  id_col   = "Phenotype_Index",
  case_col = "Counts_GATA2",
  ctrl_col = "Counts_UKB",
  N_cases = N_cases,
  N_controls = N_controls,
  fisher_for_small = TRUE,
  show_progress = TRUE
) %>%
  dplyr::left_join(
    select_pheno %>% dplyr::select(Phenotype_Index, Phenotype, Phenotype_Description),
    by = "Phenotype_Index"
  )

# ---------- Forest plot (significant single phenotypes) ----------
single_pheno_sig <- single_pheno %>% filter(padj < 0.01)
pdf(file.path(path_plots, "Significant_Single Phenotypes.pdf"), width = 10, height = 10)
ggplot(single_pheno_sig, aes(x = RR, y = reorder(Phenotype_Description, RR), color = padj < 0.05)) +
  geom_point() +
  geom_errorbar(aes(xmin = RR_Lower_CI, xmax = RR_Upper_CI), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(trans = "log10") +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey"),
                     labels = c("FALSE" = "Not Sig.", "TRUE" = "Bonferroni < 0.05")) +
  labs(x = "Relative Risk (log scale)", y = "Phenotype",
       title = "Forest plot of Relative Risk in Adults",
       color = "Significant") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")
dev.off()

# ---------- Sensitivity/Specificity/RR/OR trend plots ----------
pdf(file.path(path_plots, "Single_Phenotype_Trends.pdf"), width = 10, height = 10)
# Sensitivity
ggplot(single_pheno, aes(x = seq_along(Sensitivity), y = Sensitivity)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Phenotype (index order)", y = "Sensitivity", title = "Sensitivity trend with CI") +
  theme_minimal()
# Specificity
ggplot(single_pheno, aes(x = seq_along(Specificity), y = Specificity)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Phenotype (index order)", y = "Specificity", title = "Specificity trend with CI") +
  theme_minimal()
# RR
single_pheno_non_na <- single_pheno %>% filter(is.finite(RR))
ggplot(single_pheno_non_na, aes(x = seq_along(RR), y = RR)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Phenotype (index order)", y = "Risk Ratio (log scale)", title = "RR trend with CI") +
  theme_minimal()
# OR
single_pheno_non_na <- single_pheno %>% filter(is.finite(OR))
ggplot(single_pheno_non_na, aes(x = seq_along(OR), y = OR)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Phenotype (index order)", y = "Odds Ratio (log scale)", title = "OR trend with CI") +
  theme_minimal()
dev.off()

# ---------- Combination processing ----------


read_and_process_combo <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  process_any_fast(
    df,
    id_col   = "Combination",
    case_col = "Counts_GATA2",
    ctrl_col = "Counts_UKB",
    N_cases = N_cases,
    N_controls = N_controls,
    fisher_for_small = TRUE,
    show_progress = TRUE
  )
}

combo_results <- lapply(comb_paths, read_and_process_combo)
names(combo_results) <- names(comb_paths)

combinations_df_2 <- combo_results$Combination_2
combinations_df_3 <- combo_results$Combination_3
combinations_df_4 <- combo_results$Combination_4


# ---------- Forest plots: top 50 by |RR-1| (avoid tiny deviations) ----------
plot_top50_forest <- function(df, k, title, out_pdf) {
  # Ensure logRR exists; if not, create it
  if (!"logRR" %in% names(df)) {
    df <- df %>% mutate(logRR = log10(RR))
  }
  
  # Also compute log CI if not already present
  if (!"logRR_Lower" %in% names(df) | !"logRR_Upper" %in% names(df)) {
    df <- df %>%
      mutate(
        logRR_Lower = log10(RR_Lower_CI),
        logRR_Upper = log10(RR_Upper_CI)
      )
  }
  
  # Take top k by effect size on log scale
  top_df <- df %>%
    arrange(desc(abs(logRR))) %>%
    slice_head(n = k)
  
  p <- ggplot(top_df, aes(x = logRR, y = reorder(Combination, logRR))) +
    geom_point() +
    geom_errorbar(aes(xmin = logRR_Lower, xmax = logRR_Upper), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # null = logRR = 0
    labs(x = "log10(Relative Risk)", y = "Phenotype Combination", title = title) +
    theme_minimal()
  
  ggsave(file.path(path_plots, out_pdf), p, width = 10, height = 6)
}

plot_top50_forest(combinations_df_2 %>% filter(padj < 0.01), 50,
                  "Forest plot of Relative Risk of top 50 Combinations of 2", "Combination_2_RR.pdf")
plot_top50_forest(combinations_df_3 %>% filter(padj < 0.01), 50,
                  "Forest plot of Relative Risk of top 50 Combinations of 3", "Combination_3_RR.pdf")
plot_top50_forest(combinations_df_4 %>% filter(padj < 0.01), 50,
                  "Forest plot of Relative Risk of top 50 Combinations of 4", "Combination_4_RR.pdf")

# ---------- Combine significant hits across 1–4 ----------
single_pheno_sig    <- single_pheno    %>% mutate(Combination_Number = "Single")        %>% filter(padj < 0.01)
combinations_df_2_sig <- combinations_df_2 %>% mutate(Combination_Number = "Combination_2") %>% filter(padj < 0.01)
combinations_df_3_sig <- combinations_df_3 %>% mutate(Combination_Number = "Combination_3") %>% filter(padj < 0.01)
combinations_df_4_sig <- combinations_df_4 %>% mutate(Combination_Number = "Combination_4") %>% filter(padj < 0.01)

combined_df_AllVar <- bind_rows(single_pheno_sig, combinations_df_2_sig, combinations_df_3_sig, combinations_df_4_sig)
All_combinations   <- bind_rows(combinations_df_2_sig, combinations_df_3_sig, combinations_df_4_sig)

# Save combined significant combinations
readr::write_csv(All_combinations, file.path(path_plots, "All_Significant_Combinations_beforeOverlapremoval_121025.csv"))
readr::write_csv(combined_df_AllVar, file.path(path_plots, "All_Significant_Combinations_And_Single_beforeOverlapRemoval_121025.csv"))


#Matching till here
# ---------- Sensitivity/Specificity by combination size ----------
comb_levels <- c("Single", "Combination_2", "Combination_3", "Combination_4")

combined_df <- combined_df_AllVar %>%
  mutate(Combination = factor(Combination_Number, levels = comb_levels))

summarize_with_ci <- function(df, metric) {
  df %>%
    group_by(Combination) %>%
    summarize(
      Mean = mean(.data[[metric]], na.rm = TRUE),
      SE   = sd(.data[[metric]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[metric]]))),
      Lower_CI = Mean - 1.96 * SE,
      Upper_CI = Mean + 1.96 * SE,
      .groups = "drop"
    )
}

plot_metric_trend <- function(df_sum, metric_label, title, out_pdf) {
  #pdf(file.path(path_plots, out_pdf), width = 8, height = 5)
  ggplot(df_sum, aes(x = Combination, y = Mean, group = 1)) +
    geom_line(linewidth = 1) +
    geom_point(shape = 17, size = 3) +
    geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.2) +
    labs(title = title, x = "", y = metric_label) +
    ylim(0, 1) +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14))
  #dev.off()
}

sens_sum <- summarize_with_ci(combined_df, "Sensitivity")
spec_sum <- summarize_with_ci(combined_df, "Specificity")

plot_metric_trend(sens_sum, "Sensitivity", "Sensitivity of Combination Analysis", "Sensitivity_by_Combination.pdf")
plot_metric_trend(spec_sum, "Specificity", "Specificity of Combination Analysis", "Specificity_by_Combination.pdf")

# ---------- Frequency of phenotype indices among significant combinations ----------
split_phenos <- function(x) trimws(unlist(strsplit(x, "_")))
sig_combo_strings <- c(
  combinations_df_2_sig$Combination,
  combinations_df_3_sig$Combination,
  combinations_df_4_sig$Combination
)
phenotype_list <- split_phenos(sig_combo_strings)
phenotype_freq_df <- as.data.frame(table(phenotype_list)) %>%
  arrange(desc(Freq)) %>%
  rename(Phenotype_Index = phenotype_list)

# add annotations and single-phenotype RR/padj
data_combined_table <- phenotype_freq_df %>%
  left_join(select_pheno_full %>% dplyr::select(Phenotype, Phenotype_Index, Phenotype_Description_Expanded, Keep_Adults),
            by = "Phenotype_Index") %>%
  left_join(single_pheno %>% dplyr::select(Phenotype_Index, RR, padj), by = "Phenotype_Index") %>%
  arrange(desc(Freq))

readr::write_csv(data_combined_table, file.path(path_plots, "Ranked_Phenotypes_Combined_Final_022526.csv"))


# ---------- Combination tables with descriptions ----------
create_combination_table <- function(data_sig, combo_label, select_pheno) {
  # data_sig: df filtered to padj < 0.01 (or not, your choice)
  num <- as.integer(gsub(".*_", "", combo_label))
  tmp <- data_sig %>% filter(Combination_Number == combo_label)
  if (!nrow(tmp)) return(tibble::tibble())
  
  # split into columns Phenotype1..PhenotypeN
  tmp <- tmp %>%
    tidyr::separate(Combination, into = paste0("Phenotype", 1:num), sep = "_", remove = FALSE)
  
  # map indices → descriptions without repeated joins
  map_desc <- setNames(select_pheno$Phenotype_Description, select_pheno$Phenotype_Index)
  for (i in seq_len(num)) {
    idx_col  <- paste0("Phenotype", i)
    desc_col <- paste0("Phenotype", i, "_desc")
    tmp[[desc_col]] <- unname(map_desc[ tmp[[idx_col]] ])
  }
  
  tmp %>%
    dplyr::select(all_of(paste0("Phenotype", 1:num, "_desc")), logRR) %>%
    dplyr::rename_with(~ gsub("_desc$", "", .x)) %>%
    dplyr::arrange(desc(logRR))
}

# Build labeled combined (so the function can filter by label)
combinations_df_2_sig <- combinations_df_2_sig %>% mutate(Combination_Number = "Combination_2")
combinations_df_3_sig <- combinations_df_3_sig %>% mutate(Combination_Number = "Combination_3")
combinations_df_4_sig <- combinations_df_4_sig %>% mutate(Combination_Number = "Combination_4")
All_combinations_labeled <- bind_rows(combinations_df_2_sig, combinations_df_3_sig, combinations_df_4_sig)

table_2 <- create_combination_table(All_combinations_labeled, "Combination_2", select_pheno)
table_3 <- create_combination_table(All_combinations_labeled, "Combination_3", select_pheno)
table_4 <- create_combination_table(All_combinations_labeled, "Combination_4", select_pheno)

readr::write_csv(table_2, file.path(path_plots, "Combinations_2_022526.csv"))
readr::write_csv(table_3, file.path(path_plots, "Combinations_3_022526.csv"))
readr::write_csv(table_4, file.path(path_plots, "Combinations_4_022526.csv"))



# ---------- Relative Risk table for All combinations----------

# Build a long-format dataframe
rr_summary_input <- bind_rows(
  data.frame(Set = "Single Phenotype", RR = single_pheno_sig$RR,logRR=single_pheno_sig$logRR),
  data.frame(Set = "Combinations of 2", RR = combinations_df_2_sig$RR,logRR=combinations_df_2_sig$logRR),
  data.frame(Set = "Combinations of 3", RR = combinations_df_3_sig$RR,logRR=combinations_df_3_sig$logRR),
  data.frame(Set = "Combinations of 4", RR = combinations_df_4_sig$RR,logRR=combinations_df_4_sig$logRR)
)

# Summarize and order by mean
rr_summary <- rr_summary_input %>%
  group_by(Set) %>%
  summarise(
    Count   = n(),
    Min     = min(logRR, na.rm = TRUE),
    Q1      = quantile(logRR, 0.25, na.rm = TRUE),
    Median  = median(logRR, na.rm = TRUE),
    Mean    = mean(logRR, na.rm = TRUE),
    Q3      = quantile(logRR, 0.75, na.rm = TRUE),
    Max     = max(logRR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))   # change to arrange(Mean) for ascending order

rr_summary





# ---------- Odds Ratio table for all combinations  ----------
# Build a long-format dataframe for OR
or_summary_input <- bind_rows(
  data.frame(Set = "Single Phenotype", OR = single_pheno_sig$OR,logOR=single_pheno_sig$logOR),
  data.frame(Set = "Combinations of 2", OR = combinations_df_2_sig$OR,logOR=combinations_df_2_sig$logOR),
  data.frame(Set = "Combinations of 3", OR = combinations_df_3_sig$OR,logOR=combinations_df_3_sig$logOR),
  data.frame(Set = "Combinations of 4", OR = combinations_df_4_sig$OR,logOR=combinations_df_4_sig$logOR)
)

# Summarize and order by mean
or_summary <- or_summary_input %>%
  group_by(Set) %>%
  summarise(
    Count   = n(),
    Min     = min(logOR, na.rm = TRUE),
    Q1      = quantile(logOR, 0.25, na.rm = TRUE),
    Median  = median(logOR, na.rm = TRUE),
    Mean    = mean(logOR, na.rm = TRUE),
    Q3      = quantile(logOR, 0.75, na.rm = TRUE),
    Max     = max(logOR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))   # change to arrange(Mean) for ascending order

or_summary





# ---------- Combination tables with descriptions ----------
rr_summary <- rr_summary_input %>%
  group_by(Set) %>%
  summarise(
    Count     = n(),
    Min       = min(logRR, na.rm = TRUE),
    Q1        = quantile(logRR, 0.25, na.rm = TRUE),
    Median    = median(logRR, na.rm = TRUE),
    Mean      = mean(logRR, na.rm = TRUE),
    Q3        = quantile(logRR, 0.75, na.rm = TRUE),
    Max       = max(logRR, na.rm = TRUE),
    Mean_log  = mean(logRR, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(desc(Mean))  # change to arrange(Mean) for ascending

# 3) OR summary (with Mean logOR)
or_summary_input <- bind_rows(
  data.frame(Set = "Single Phenotype", OR = single_pheno_sig$OR, logOR = single_pheno_sig$logOR),
  data.frame(Set = "Combinations of 2", OR = combinations_df_2_sig$OR, logOR = combinations_df_2_sig$logOR),
  data.frame(Set = "Combinations of 3", OR = combinations_df_3_sig$OR, logOR = combinations_df_3_sig$logOR),
  data.frame(Set = "Combinations of 4", OR = combinations_df_4_sig$OR, logOR = combinations_df_4_sig$logOR)
)

or_summary <- or_summary_input %>%
  group_by(Set) %>%
  summarise(
    Count     = n(),
    Min       = min(logOR, na.rm = TRUE),
    Q1        = quantile(logOR, 0.25, na.rm = TRUE),
    Median    = median(logOR, na.rm = TRUE),
    Mean      = mean(logOR, na.rm = TRUE),
    Q3        = quantile(logOR, 0.75, na.rm = TRUE),
    Max       = max(logOR, na.rm = TRUE),
    Mean_log  = mean(logOR, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(desc(Mean))  # change to arrange(Mean) for ascending

# Optional: print
rr_summary
or_summary





# --- Relative Risk violin ---
rr_means <- rr_summary_input %>%
  group_by(Set) %>%
  summarise(log = mean(logRR, na.rm = TRUE), .groups = "drop")


rr_summary_input <- rr_summary_input %>%
  dplyr::left_join(rr_means, by = "Set") %>%
  dplyr::mutate(Set = reorder(Set, logRR))   # order by mean

p_rr_violin <- ggplot(rr_summary_input, aes(x = logRR, y = Set)) +
  geom_violin(fill = "skyblue", alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    title = "Relative Risk Distribution",
    x = "log10(Relative Risk)", y = NULL
  ) +
  theme_minimal(base_size = 12)

# Print
p_rr_violin +
  scale_x_continuous(
    breaks = seq(-1, 7, by = 1),
  )

# --- Odds Ratio violin ---
or_means <- or_summary_input %>%
  group_by(Set) %>%
  summarise(meanOR = mean(logOR, na.rm = TRUE), .groups = "drop")

or_summary_input <- or_summary_input %>%
  left_join(or_means, by = "Set") %>%
  mutate(Set = reorder(Set, logOR))

p_or_violin <- ggplot(or_summary_input, aes(x = logOR, y = Set)) +
  geom_violin(fill = "salmon", alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.5,outlier.shape=NA) +
  geom_vline(xintercept = 3, linetype = "dashed") +
    labs(
    title = "Odds Ratio Distribution",
    x = "Odds Ratio (log scale)", y = NULL
  ) +
  theme_minimal(base_size = 12)

p_or_violin +
  scale_x_continuous(
    breaks = seq(-1, 7, by = 1),
  )







# --- Sensitivity summary ---
sens_input <- bind_rows(
  data.frame(Set = "Single Phenotype", Sensitivity = single_pheno_sig$Sensitivity),
  data.frame(Set = "Combinations of 2", Sensitivity = combinations_df_2_sig$Sensitivity),
  data.frame(Set = "Combinations of 3", Sensitivity = combinations_df_3_sig$Sensitivity),
  data.frame(Set = "Combinations of 4", Sensitivity = combinations_df_4_sig$Sensitivity)
)

sens_summary <- sens_input %>%
  group_by(Set) %>%
  summarise(
    Count   = n(),
    Min     = min(Sensitivity, na.rm = TRUE),
    Q1      = quantile(Sensitivity, 0.25, na.rm = TRUE),
    Median  = median(Sensitivity, na.rm = TRUE),
    Mean    = mean(Sensitivity, na.rm = TRUE),
    Q3      = quantile(Sensitivity, 0.75, na.rm = TRUE),
    Max     = max(Sensitivity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))   # order by mean sensitivity

# --- Specificity summary ---
spec_input <- bind_rows(
  data.frame(Set = "Single Phenotype", Specificity = single_pheno_sig$Specificity),
  data.frame(Set = "Combinations of 2", Specificity = combinations_df_2_sig$Specificity),
  data.frame(Set = "Combinations of 3", Specificity = combinations_df_3_sig$Specificity),
  data.frame(Set = "Combinations of 4", Specificity = combinations_df_4_sig$Specificity)
)

spec_summary <- spec_input %>%
  group_by(Set) %>%
  summarise(
    Count   = n(),
    Min     = min(Specificity, na.rm = TRUE),
    Q1      = quantile(Specificity, 0.25, na.rm = TRUE),
    Median  = median(Specificity, na.rm = TRUE),
    Mean    = mean(Specificity, na.rm = TRUE),
    Q3      = quantile(Specificity, 0.75, na.rm = TRUE),
    Max     = max(Specificity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean))   # order by mean specificity

# View results
sens_summary
spec_summary


#Table for Percentage of People having the phenotype Before removing overlapping phenotypes


# Helper: IDs present for ≥1 of the given combinations (works for any k)
ids_for_combos <- function(df_wide, combos_df, comb_col = "Combination", sep = "_") {
  if (nrow(combos_df) == 0) return(character(0))
  pheno_long <- df_wide %>%
    pivot_longer(everything(), names_to = "Phenotype", values_to = "ParticipantID") %>%
    filter(!is.na(ParticipantID))
  
  combos <- combos_df %>%
    mutate(!!comb_col := str_replace_all(.data[[comb_col]], "\\s+", "")) %>%
    pull(!!comb_col)
  
  ids <- map(combos, function(x) {
    phs <- str_split(x, sep, simplify = TRUE) %>% as.character()
    pheno_long %>%
      filter(Phenotype %in% phs) %>%
      group_by(ParticipantID) %>%
      summarise(n = n_distinct(Phenotype), .groups = "drop") %>%
      filter(n == length(phs)) %>%
      pull(ParticipantID)
  })
  unique(unlist(ids))
}

# --- Denominator: all unique participants with ≥1 phenotype in the table ---
all_ids <- phenotype_gata2_adult %>% unlist() %>% na.omit() %>% unique()

# --- Sets for each combo size ---
ids2 <- ids_for_combos(phenotype_gata2_adult, combinations_df_2_sig, "Combination", "_")
ids3 <- ids_for_combos(phenotype_gata2_adult, combinations_df_3_sig, "Combination", "_")
ids4 <- ids_for_combos(phenotype_gata2_adult, combinations_df_4_sig, "Combination", "_")

# --- Counts and percents (overlapping: “has at least one …”) ---
n_all <- 212
summary_tbl <- tribble(
  ~Category,            ~Count, ~CountPhenoCombo,    ~Percent,
  "Single phenotype(All)",    length(all_ids), nrow(single_pheno),                100 * length(all_ids) / n_all,        # always 100% of those observed
  "Combinations of 2(Significant Only)",      length(ids2), nrow(combinations_df_2_sig)  ,       100 * length(ids2) / n_all,
  "Combinations of 3(Significant Only)",      length(ids3), nrow(combinations_df_3_sig)   ,      100 * length(ids3) / n_all,
  "Combinations of 4(Significant Only)",      length(ids4), nrow(combinations_df_4_sig)    ,     100 * length(ids4) / n_all
)

summary_tbl


# Get the phenotype columns of interest
pheno_cols <- intersect(colnames(phenotype_gata2_adult), single_pheno_sig$Phenotype_Index)

# Collect unique IDs only from those columns
all_ids <- phenotype_gata2_adult %>%
  dplyr::select(all_of(pheno_cols)) %>%
  unlist() %>%
  na.omit() %>%
  unique()

all_ids
length(all_ids)  # number of unique participants

# --- Counts and percents (Significant Phenotypes only) ---
n_all <- 212
summary_tbl_2 <- tribble(
  ~Category,            ~Count, ~CountPhenoCombo,    ~Percent,
  "Single phenotype(Significant Only)",    length(all_ids),nrow(single_pheno_sig) ,                100 * length(all_ids) / n_all,        # always 100% of those observed
  "Combinations of 2(Significant Only)",      length(ids2), nrow(combinations_df_2_sig)  ,       100 * length(ids2) / n_all,
  "Combinations of 3(Significant Only)",      length(ids3), nrow(combinations_df_3_sig)   ,      100 * length(ids3) / n_all,
  "Combinations of 4(Significant Only)",      length(ids4), nrow(combinations_df_4_sig)    ,     100 * length(ids4) / n_all
)

summary_tbl_2







adult_id_all<-rownames(phenotype_gata2)[phenotype_gata2$Age>=18]
missing<-setdiff(adult_id_all,all_ids)
adult_asymptomatic_ids <- rownames(phenotype_gata2)[ df$Age >= 18 & !df$has_any_pheno ]

variants_asymp<-gata2_survey_data$VariantType[gata2_survey_data$prevalence_patient_deidentifier_lab_family_patient%in%adult_asymptomatic_ids]
variants_asymp

variants_asymp_variant<-gata2_survey_data$prevalence_variant_name[gata2_survey_data$prevalence_patient_deidentifier_lab_family_patient%in%adult_asymptomatic_ids]
variants_asymp_variant


ped_asymptomatic_ids <- rownames(phenotype_gata2)[ df$Age < 18 & !df$has_any_pheno ]

variants_asymp_ped<-gata2_survey_data$VariantType[gata2_survey_data$prevalence_patient_deidentifier_lab_family_patient%in%ped_asymptomatic_ids]
variants_asymp_ped

variants_asymp_variant_ped<-gata2_survey_data$prevalence_variant_name[gata2_survey_data$prevalence_patient_deidentifier_lab_family_patient%in%ped_asymptomatic_ids]
variants_asymp_variant_ped


###LABEL ALL COMBINATIONS

# Inputs expected:
# 1) combos_df: has columns "Combination" (e.g., "P1_P2") and OR (or Odds_Ratio)
# 2) select_pheno: has columns "Phenotype_Index" (e.g., "P1") and "Phenotype" (name)

# --- 0) Make sure the column names are what we expect -------------------------
# Use "OR" if it exists, else try "Odds_Ratio"
combos_df <- All_combinations_labeled %>%
  mutate(OR = OR)

# If select_pheno happens to have slightly different names, try to normalize:
if (!"Phenotype_Index" %in% names(select_pheno)) {
  # fallback if user’s table has 'Phenotype' holding the index
  # and a second column with the name (e.g., Phenotype_Description)
  idx_col <- intersect(names(select_pheno), c("Phenotype_Index","Phenotype","Index"))[1]
  name_col <- intersect(names(select_pheno), c("Phenotype","Phenotype_Description","Name"))[2]
  if (!is.na(idx_col) && !is.na(name_col)) {
    select_pheno <- select_pheno %>%
      transmute(Phenotype_Index = .data[[idx_col]],
                Phenotype       = .data[[name_col]])
  } else {
    stop("`select_pheno` must contain phenotype indices and names.")
  }
}

# Clean any whitespace in lookups
select_pheno <- select_pheno_full %>%
  mutate(Phenotype_Index = str_trim(Phenotype_Index),
         Phenotype       = str_trim(Phenotype))


# --- 2) Combinations_Labels ---------------------------------------------------
# Split "Combination" on "_" → map each index to its Phenotype → build a list column of "Pi-Name"
lookup_name <- select_pheno %>% dplyr::select(Phenotype_Index, Phenotype_Description_Expanded)


# lookup table of phenotype names
lookup_name <- select_pheno %>% dplyr::select(Phenotype_Index, Phenotype_Description_Expanded)

single_df<-combined_df_AllVar[combined_df_AllVar$Combination_Number=="Single",]
single_df_2<-single_df[,c("Phenotype_Index","OR","logOR","Combination_Number")]
colnames(single_df_2)[1]<-"Combination"
combo_df_2<-combos_df[,c("Combination","OR","logOR","Combination_Number")]
combos_df_3<-rbind(combo_df_2,single_df_2)



combos_df <- combos_df_3 %>%
  mutate(Indices = str_split(Combination, "_")) %>%                  # split indices
  unnest_longer(Indices, values_to = "Phenotype_Index") %>%
  left_join(lookup_name, by = "Phenotype_Index") %>%
  mutate(Label = if_else(
    is.na(Phenotype_Description_Expanded),
    Phenotype_Index,
    paste0(Phenotype_Index, "-", Phenotype_Description_Expanded)
  )) %>%
  group_by(Combination) %>%
  summarise(
    OR                  = dplyr::first(OR),               # or unique(OR)[1]
    logOR               = dplyr::first(logOR),
    Combination_Number  = dplyr::first(Combination_Number),
    Combinations_Labels = paste(unique(Label), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(logOR))

# check a few rows
head(combos_df, 5)

write.csv(combos_df,file.path(path_plots, "All_Combinations_Labelled_022526_BeforeFiltering.csv"))

# Result: combos_df has columns: Combination, OR, logOR, Combinations_Labels (a list)
# Example to see a few:
print(combos_df %>% slice_head(n = 5))
# To expand labels as a readable string instead of list:
# combos_df %>% mutate(Combinations_Labels_str = sapply(Combinations_Labels, \(x) paste0("(", paste(x, collapse=", "), ")")))


# ---------- Remove overlapping phenotypes ----------
#Refine_Combinations by removing overlapping phenotypes
overlap_matrix<-read_excel("GATA2_Phenotypes_Final_v8.xlsx",sheet="Sheet1",col_names = FALSE)
concat<-c(overlap_matrix[,"...2"],"_",overlap_matrix[,"...3"])
concat <- paste(as.character(overlap_matrix[["...2"]]),
                as.character(overlap_matrix[["...3"]]),
                sep = "_")

concat <- concat[-c(1:3)]


overlap_matrix_label<-as.matrix.data.frame(overlap_matrix[4:117,4:117])

colnames(overlap_matrix_label)<-concat                           
row.names(overlap_matrix_label)<-concat

overlap_matrix_unlabel<-as.matrix.data.frame(overlap_matrix[4:117,4:117])

colnames(overlap_matrix_unlabel)<-overlap_matrix[3,4:117]
row.names(overlap_matrix_unlabel)<-colnames(overlap_matrix_unlabel)



# ensure matrix form
om <- as.matrix(overlap_matrix_unlabel)

combos_df_filtered <- combos_df %>%
  rowwise() %>%
  mutate(
    has_overlap = {
      ids <- str_split(Combination, "_")[[1]] |> unique()
      # if <2 phenotypes, no pair to check
      if (length(ids) < 2) FALSE else {
        # keep only ids present in the matrix
        ids <- intersect(ids, rownames(om))
        if (length(ids) < 2) FALSE else {
          # check the upper triangle of the submatrix for any 1
          any(om[ids, ids, drop = FALSE][upper.tri(om[ids, ids, drop = FALSE])] == 1)
        }
      }
    }
  ) %>%
  ungroup() %>%
  filter(!has_overlap) %>%     # <-- remove combos with any overlap
  dplyr::select(-has_overlap)

combos_df_filtered$Combination_Number<-gsub("Combination_2","Combinations of 2",combos_df_filtered$Combination_Number)
combos_df_filtered$Combination_Number<-gsub("Combination_3","Combinations of 3",combos_df_filtered$Combination_Number)
combos_df_filtered$Combination_Number<-gsub("Combination_4","Combinations of 4",combos_df_filtered$Combination_Number)



p_or_violin <- ggplot(
  combos_df_filtered %>% mutate(
    Combination_Number = fct_reorder(
      as.factor(Combination_Number),
      logOR,
      .fun = ~ mean(., na.rm = TRUE),
      .desc = FALSE
    )
  ),
  aes(x = logOR, y = Combination_Number)
) +
  geom_violin(fill = "salmon", alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.6,outlier.shape=NA) +
  geom_vline(xintercept = 3, linetype = "dashed",linewidth = 2,color="red") +  # <- if using log10(OR), consider xintercept = 0
  labs(title = "Odds Ratio Distribution Filtered Overlapping Phenotypes",
       x = "log10(Odds Ratio)", y = NULL) +
  theme_minimal(base_size = 12)

p_or_violin +scale_x_continuous(
  breaks = seq(-1, 7, by = 1),
)


combos_df_filtered$Combination_Number <- factor(
  combos_df_filtered$Combination_Number,
  levels = c(
    "Single",
    "Combinations of 2",
    "Combinations of 3",
    "Combinations of 4"
  )
)


#Test2
pdf(file.path(path_plots, "Violin_Plot_Final.pdf"), width = 7, height = 5)
ggplot(combos_df_filtered, aes(x = logOR, y = Combination_Number)) +
  
  geom_violin(
    fill = "#E64B35",
    alpha = 0.35,
    scale = "width",
    color = "black",     # outline color
    linewidth = 0.4      # THICKER edge
  ) +
  
  geom_boxplot(
    width = 0.08,
    outlier.size = 0.5,
    alpha = 0.6,outlier.shape=NA
  ) +
  
  geom_vline(
    xintercept = 3,
    linetype = "dashed",
    linewidth = 1.2,
    color = "red"
  ) +
  
  labs(
    title = "Odds Ratio Distribution",
    x = "log10(Odds Ratio)",
    y = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  scale_x_continuous(breaks = seq(-1, 7, by = 1))

dev.off()


#Updated Summary after filtering
aggregate(logOR ~ Combination_Number, data = combos_df_filtered,
          FUN = function(x) c(n = sum(!is.na(x)),
                              mean = mean(x, na.rm=TRUE),
                              sd = sd(x, na.rm=TRUE),
                              median = median(x, na.rm=TRUE),
                              IQR = IQR(x, na.rm=TRUE),
                              min = min(x, na.rm=TRUE),
                              max = max(x, na.rm=TRUE)))




# --- split singles vs combinations ---
singles <- combined_df_AllVar %>%
  dplyr::select(Phenotype_Index, Counts_GATA2, Counts_UKB) %>%
  rename(Combination = Phenotype_Index)

combinations <- combined_df_AllVar %>%
  dplyr::select(Combination, Counts_GATA2, Counts_UKB)

# --- join depending on Combination_Number ---
combos_df_with_counts <- combos_df_filtered %>%
  left_join(
    bind_rows(singles, combinations),
    by = "Combination"
  )

# assume your mapping table looks like:
lookup <- select_pheno[,c("Phenotype","Phenotype_Description_v2","Phenotype_Index")]
lookup$dict<-paste0(lookup$Phenotype_Index,"-",lookup$Phenotype_Description_v2)
lookup<-lookup[,c("Phenotype","dict")]


# 0) Build a lookup table Phenotype_Index -> New_Category
#    (works whether your column is named Phenotype_Index or Phenotype)
lookup_cat <- select_pheno %>%
  dplyr::select(any_of(c("Phenotype_Index", "Phenotype")), New_Category) %>%
  rename(Phenotype_Index = 1) %>%
  distinct()

lookup_cat_merge<-merge(lookup,lookup_cat)

# 1) Explode each Combination into its component P-indices, map to categories,
#    then collapse back to unique, order-preserving categories per Combination
combination_categories <- combos_df_with_counts %>%
  dplyr::select(Combination) %>%
  distinct() %>%
  mutate(Phenotype_Index_tmp = Combination) %>%
  separate_rows(Phenotype_Index_tmp, sep = "_") %>%
  rename(Phenotype_Index = Phenotype_Index_tmp) %>%
  left_join(lookup_cat, by = "Phenotype_Index") %>%
  group_by(Combination) %>%
  summarise(
    # order-preserving unique
    Category_List = list({ x <- New_Category; x[!is.na(x) & !duplicated(x)] }),
    # handy string version for quick viewing/CSV export
    Category      = paste({ x <- New_Category; x[!is.na(x) & !duplicated(x)] }, collapse = "; "),
    Category_Count = length({ x <- New_Category; x[!is.na(x) & !duplicated(x)] }),
    .groups = "drop"
  )

# 2) Join back to your main table
combos_df_with_counts_and_cat <- combos_df_with_counts %>%
  left_join(combination_categories, by = "Combination")

###Plot Top co-occuring categories Categories
cooc_df <- combos_df_with_counts_and_cat %>%
  filter(!is.na(Category), Category != "") %>%
  mutate(
    Category_Set = Category %>%
      str_split(";") %>%
      map(~ .x |> str_trim() |> discard(~ .x == "") |> unique() |> sort() |> paste(collapse = " + "))
  ) %>%
  mutate(Category_Set = unlist(Category_Set)) %>%
  count(Category_Set, name = "Frequency") %>%
  arrange(desc(Frequency)) %>%
  mutate(
    # how many categories are in this set (for coloring/filtering)
    Set_Size = str_count(Category_Set, fixed(" + ")) + 1,
    # nicer label for hover (line breaks)
    Category_Label = gsub(" \\+ ", "\n• ", paste0("• ", Category_Set))
  )


# Interactive plot of co-occuring category-----
library(plotly)
library(htmlwidgets)

p <- plot_ly(
  cooc_df,
  x = ~Frequency,
  y = ~reorder(Category_Set, Frequency),
  type = "bar",
  color = ~factor(Set_Size),
  hovertemplate = paste(
    "<b>Category set</b>:<br>%{customdata}<br>",
    "<b>Frequency</b>: %{x}<extra></extra>"
  ),
  customdata = ~Category_Set
) %>%
  layout(
    title = "Category Co-occurrence (set-level)",
    xaxis = list(title = "Frequency"),
    yaxis = list(title = "Category set"),
    legend = list(title = list(text = "Set size"))
  )

# Save and open in browser
htmlwidgets::saveWidget(as_widget(p), "category_cooccurrence.html", selfcontained = TRUE)
browseURL("category_cooccurrence.html")


# Get the phenotype columns of interest
pheno_cols <- intersect(colnames(phenotype_gata2_adult), single_pheno_sig$Phenotype_Index[single_pheno_sig$Phenotype_Index%in%combos_df_filtered$Combination])

# Collect unique IDs only from those columns
all_ids <- phenotype_gata2_adult %>%
  dplyr::select(all_of(pheno_cols)) %>%
  unlist() %>%
  na.omit() %>%
  unique()


# --- Sets for each combo size ---
ids2 <- ids_for_combos(phenotype_gata2_adult, combinations_df_2_sig[combinations_df_2_sig$Combination%in%combos_df_filtered$Combination,], "Combination", "_")
ids3 <- ids_for_combos(phenotype_gata2_adult, combinations_df_3_sig[combinations_df_3_sig$Combination%in%combos_df_filtered$Combination,], "Combination", "_")
ids4 <- ids_for_combos(phenotype_gata2_adult, combinations_df_4_sig[combinations_df_4_sig$Combination%in%combos_df_filtered$Combination,], "Combination", "_")



# --- Counts and percents (overlapping: “has at least one …”) ---
n_all <- 212
summary_tbl_final <- tribble(
  ~Category,            ~Count, ~CountPhenoCombo,    ~Percent,
  "Single phenotype(Significant Only)",    length(all_ids), nrow(single_pheno_sig[single_pheno_sig$Phenotype_Index%in%combos_df_filtered$Combination,]),                100 * length(all_ids) / n_all,        # always 100% of those observed
  "Combinations of 2(Significant Only)",      length(ids2), nrow(combinations_df_2_sig[combinations_df_2_sig$Combination%in%combos_df_filtered$Combination,])  ,       100 * length(ids2) / n_all,
  "Combinations of 3(Significant Only)",      length(ids3), nrow(combinations_df_3_sig[combinations_df_3_sig$Combination%in%combos_df_filtered$Combination,])   ,      100 * length(ids3) / n_all,
  "Combinations of 4(Significant Only)",      length(ids4), nrow(combinations_df_4_sig[combinations_df_4_sig$Combination%in%combos_df_filtered$Combination,])    ,     100 * length(ids4) / n_all
)

summary_tbl_final

writexl::write_xlsx(combos_df_with_counts_and_cat[,-8],file.path(path_plots, "All_Combinations_Labelled_Filtered_102125.xlsx"))


combinations_df_2_sig_filt<-combinations_df_2_sig[combinations_df_2_sig$Combination%in%combos_df_filtered$Combination,]
combinations_df_3_sig_filt<-combinations_df_3_sig[combinations_df_3_sig$Combination%in%combos_df_filtered$Combination,]
combinations_df_4_sig_filt<-combinations_df_4_sig[combinations_df_4_sig$Combination%in%combos_df_filtered$Combination,]



#Combination 2

thresholds <- c(-1,0,1,2,3,4,5,6)

# ---- helper: ensure we have an OR column; fall back to exp(logOR) ----
get_or_column <- function(df) {
  nms <- names(df)
  if ("OR" %in% nms) {
    df %>% mutate(OR = as.numeric(OR))
  } else if ("Odds_Ratio" %in% nms) {
    df %>% mutate(OR = as.numeric(Odds_Ratio))
  } else if ("odds_ratio" %in% nms) {
    df %>% mutate(OR = as.numeric(odds_ratio))
  } else if ("logOR" %in% nms) {
    df %>% mutate(OR = 10(as.numeric(logOR)))
  } else {
    stop("No OR or logOR column found in combinations_df_2(_sig).")
  }
}

# ---- helper: map each combination -> vector of ParticipantIDs that meet ALL phenotypes ----
ids_by_combo <- function(df_wide, combos_df, comb_col = "Combination", sep = "_") {
  if (!comb_col %in% names(combos_df)) stop("`comb_col` not in combos_df.")
  # Only pivot columns we actually need
  needed <- combos_df %>%
    dplyr::pull(!!comb_col) %>%
    str_remove_all("\\s+") %>%
    str_split(sep) %>%
    unlist() %>%
    unique()
  needed <- intersect(needed, colnames(df_wide))
  
  pheno_long <- df_wide %>%
    dplyr::select(all_of(needed)) %>%
    pivot_longer(everything(), names_to = "Phenotype", values_to = "ParticipantID") %>%
    filter(!is.na(ParticipantID))
  
  combos_df %>%
    transmute(Combination = str_remove_all(.data[[comb_col]], "\\s+")) %>%
    distinct() %>%
    mutate(
      ids = map(Combination, ~{
        phs <- str_split(.x, sep, simplify = TRUE) %>% as.character()
        phs <- phs[phs != ""]
        pheno_long %>%
          filter(Phenotype %in% phs) %>%
          group_by(ParticipantID) %>%
          summarise(n = n_distinct(Phenotype), .groups = "drop") %>%
          filter(n == length(phs)) %>%
          pull(ParticipantID) %>%
          unique()
      })
    )
}

# ---- inputs: k=2 combinations with OR (or logOR) ----
comb2_df <- get_or_column(combinations_df_2_sig_filt)   # <- change name if needed

# Map each combination to its matching IDs, and attach the OR
ids2_tbl <- ids_by_combo(phenotype_gata2_adult, comb2_df, "Combination") %>%
  left_join(comb2_df %>% dplyr::select(Combination, logOR), by = "Combination")

# Denominator for percent (use your earlier all_ids if available, else 212)
n_total <- 212

# ---- compute coverage for each OR cutoff (≥ cutoff) ----
cov2 <- map_dfr(thresholds, function(t) {
  ids_union <- ids2_tbl %>%
    filter(!is.na(logOR), logOR >= t) %>%
    pull(ids) %>%
    unlist() %>%
    unique()
  tibble(
    OR_cutoff = t,
    n = length(ids_union),
    pct = 100 * length(ids_union) / n_total
  )
})

# ---- plot: counts on y, percent annotated ----
p2_counts <- ggplot(cov2, aes(x = factor(OR_cutoff), y = n)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)),
            vjust = -0.25, size = 3) +
  labs(
    title = "k = 2: Participants with ≥1 combination at OR ≥ cutoff",
    x = "OR cutoff (≥)",
    y = "Number of participants"
  ) +
  theme_minimal(base_size = 12)

p2_counts
# ggsave("k2_coverage_counts.png", p2_counts, width = 7, height = 4.5, dpi = 300)


# ---- inputs: k=3 combinations with OR (or logOR) ----
comb3_df <- get_or_column(combinations_df_3_sig_filt)

# Map each combination to its matching IDs, and attach the OR
ids3_tbl <- ids_by_combo(phenotype_gata2_adult, comb3_df, "Combination") %>%
  dplyr::left_join(comb3_df %>% dplyr::select(Combination, logOR), by = "Combination")

# Denominator for percent (use your earlier all_ids if available, else 212)
n_total <- 212

# ---- compute coverage for each OR cutoff (≥ cutoff) ----
cov3 <- purrr::map_dfr(thresholds, function(t) {
  ids_union <- ids3_tbl %>%
    dplyr::filter(!is.na(logOR), logOR >= t) %>%
    dplyr::pull(ids) %>%
    unlist() %>%
    unique()
  tibble::tibble(
    OR_cutoff = t,
    n = length(ids_union),
    pct = 100 * length(ids_union) / n_total
  )
})

# ---- plot: counts on y, percent annotated ----
p3_counts <- ggplot2::ggplot(cov3, ggplot2::aes(x = factor(OR_cutoff), y = n)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d (%.1f%%)", n, pct)),
                     vjust = -0.25, size = 3) +
  ggplot2::labs(
    title = "k = 3: Participants with ≥1 combination at OR ≥ cutoff",
    x = "OR cutoff (≥)",
    y = "Number of participants"
  ) +
  ggplot2::theme_minimal(base_size = 12)

p3_counts
# ggsave("k3_coverage_counts.png", p3_counts, width = 7, height = 4.5, dpi = 300)

#
# ---- inputs: k=4 combinations with OR (or logOR) ----
comb4_df <- get_or_column(combinations_df_4_sig_filt)

# Map each combination to its matching IDs, and attach the OR
ids4_tbl <- ids_by_combo(phenotype_gata2_adult, comb4_df, "Combination") %>%
  dplyr::left_join(comb4_df %>% dplyr::select(Combination, logOR), by = "Combination")

# Denominator for percent (use your earlier all_ids if available, else 212)
n_total <-212

# ---- compute coverage for each OR cutoff (≥ cutoff) ----
cov4 <- purrr::map_dfr(thresholds, function(t) {
  ids_union <- ids4_tbl %>%
    dplyr::filter(!is.na(logOR), logOR >= t) %>%
    dplyr::pull(ids) %>%
    unlist() %>%
    unique()
  tibble::tibble(
    OR_cutoff = t,
    n = length(ids_union),
    pct = 100 * length(ids_union) / n_total
  )
})

# ---- plot: counts on y, percent annotated ----
p4_counts <- ggplot2::ggplot(cov4, ggplot2::aes(x = factor(OR_cutoff), y = n)) +
  ggplot2::geom_col() +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%d (%.1f%%)", n, pct)),
                     vjust = -0.25, size = 3) +
  ggplot2::labs(
    title = "k = 4: Participants with ≥1 combination at OR ≥ cutoff",
    x = "OR cutoff (≥)",
    y = "Number of participants"
  ) +
  ggplot2::theme_minimal(base_size = 12)

p4_counts
# ggsave("k3_coverage_counts.png", p3_counts, width = 7, height = 4.5, dpi = 300)



#With Combination Number 
# Make sure we have a clean OR and Combination column
ensure_or <- function(df) {
  out <- df
  out %>%
    mutate(Combination = str_replace_all(Combination, "\\s+", "")) %>%
    distinct(Combination, .keep_all = TRUE)
}

# Standardize your dfs
comb2_df <- ensure_or(combinations_df_2_sig_filt)
comb3_df <- ensure_or(combinations_df_3_sig_filt)
comb4_df <- ensure_or(combinations_df_4_sig_filt)


count_combos_cum <- function(df, thresholds, k_label) {
  map_dfr(thresholds, function(t) {
    tibble(
      OR_cutoff = t,
      n_combos  = df %>%
        filter(!is.na(logOR), logOR >= t) %>%
        distinct(Combination) %>%           # <-- get unique combos
        nrow(),                             # <-- count them
      k = k_label
    )
  })
}

combo_counts2_cum <- count_combos_cum(comb2_df, thresholds, "k=2")
combo_counts3_cum <- count_combos_cum(comb3_df, thresholds, "k=3")
combo_counts4_cum <- count_combos_cum(comb4_df, thresholds, "k=4")

combo_counts2_cum
combo_counts3_cum
combo_counts4_cum

combo_counts_cum <- bind_rows(combo_counts2_cum, combo_counts3_cum,combo_counts4_cum)


# --- Plot cumulative counts by cutoff (grouped bars) ---
p_combo_cum <- ggplot(combo_counts_cum, aes(x = factor(OR_cutoff), y = n_combos, fill = k)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = n_combos),
            position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
  labs(title = "Number of combinations meeting OR ≥ cutoff",
       x = "OR cutoff (≥)", y = "# of combinations", fill = NULL) +
  theme_minimal(base_size = 12)
p_combo_cum




# --- reuse your thresholds and n_total ---
n_total <- 212

# --- ids_by_combo must exist from earlier; if not, define it again ---
# ids_by_combo(df_wide, combos_df, comb_col = "Combination", sep = "_") ...

# Build ids tables with OR attached (k=2,3,4) using your comb*_df already created
ids2_tbl <- ids_by_combo(phenotype_gata2_adult, comb2_df, "Combination") %>%
  dplyr::left_join(comb2_df %>% dplyr::select(Combination, logOR), by = "Combination")
ids3_tbl <- ids_by_combo(phenotype_gata2_adult, comb3_df, "Combination") %>%
  dplyr::left_join(comb3_df %>% dplyr::select(Combination, logOR), by = "Combination")
ids4_tbl <- ids_by_combo(phenotype_gata2_adult, comb4_df, "Combination") %>%
  dplyr::left_join(comb4_df %>% dplyr::select(Combination, logOR), by = "Combination")

# Coverage (% covered) at/above each cutoff, per k
cov_k <- function(ids_tbl, k_label) {
  map_dfr(thresholds, function(t) {
    ids_union <- ids_tbl %>% 
      dplyr::filter(!is.na(logOR), logOR >= t) %>%
      dplyr::pull(ids) %>% unlist() %>% unique()
    tibble(OR_cutoff = t, pct = 100 * length(ids_union) / n_total, k = k_label)
  })
}
cov2 <- cov_k(ids2_tbl, "k=2")
cov3 <- cov_k(ids3_tbl, "k=3")
cov4 <- cov_k(ids4_tbl, "k=4")
cov_cum <- bind_rows(cov2, cov3, cov4)

# Join coverage onto your counts
counts_cov <- combo_counts_cum %>%
  dplyr::left_join(cov_cum, by = c("OR_cutoff","k"))

# ---- OPTION A: overlay coverage as a line with a secondary % axis ----
scale_fac <- max(counts_cov$n_combos, na.rm = TRUE) / 100  # map 100% to max bar height

p_combo_cum2 <- ggplot(counts_cov, aes(x = factor(OR_cutoff))) +
  geom_col(aes(y = n_combos, fill = k), position = position_dodge(width = 0.9)) +
  geom_text(aes(y = n_combos, label = n_combos, group = k),
            position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
  geom_line(aes(y = pct * scale_fac, color = k, group = k),
            position = position_dodge(width = 0.9)) +
  geom_point(aes(y = pct * scale_fac, color = k, group = k),
             position = position_dodge(width = 0.9)) +
  scale_y_continuous(
    name = "# of combinations",
    sec.axis = sec_axis(~ . / scale_fac, name = "% covered")
  ) +
  labs(title = "Combinations meeting OR ≥ cutoff (bars) + Coverage (line)",
       x = "OR cutoff (≥)", fill = NULL, color = NULL) +
  theme_minimal(base_size = 12)

p_combo_cum2
# ggsave("combo_counts_with_coverage.png", p_combo_cum2, width = 8, height = 5, dpi = 300)


###Option 2
# keep scale_fac and counts_cov from before
p_combo_cum2 <- ggplot(counts_cov, aes(x = factor(OR_cutoff))) +
  # --- Coverage line (now primary axis) ---
  geom_line(
    aes(y = pct, color = k, group = k),
    size = 0.9, position = position_dodge(width = 0.9)
  ) +
  geom_point(
    aes(y = pct, color = k, group = k),
    size = 2.5, position = position_dodge(width = 0.9)
  ) +
  
  # --- Bars for # of combinations (scaled to secondary axis) ---
  geom_col(
    aes(y = n_combos / max(n_combos) * 100, fill = k),
    position = position_dodge(width = 0.9),
    alpha = 0.6
  ) +
  geom_text(
    aes(y = n_combos / max(n_combos) * 100, label = n_combos, group = k),
    position = position_dodge(width = 0.9),
    vjust = -0.25, size = 3
  ) +
  
  # --- Axis scaling: primary = coverage (%), secondary = # combinations ---
  scale_y_continuous(
    name = "% participants who satisfy criteria",
    sec.axis = sec_axis(
      ~ . * max(counts_cov$n_combos) / 100,
      name = "# of combinations",
      labels = scales::comma
    ),breaks = seq(0, 100, by = 10)
  ) +
  
  labs(
    x = expression("OR cutoff (log" [10] * " scale)"),
    fill = NULL,
    color = NULL
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "black", face = "bold"),
    axis.title.y.right = element_text(color = "grey40"),
    axis.text.y.right  = element_text(color = "grey40"),
    panel.grid.minor = element_blank()
  )

pdf(file.path(path_plots, "LogOR_Cutoff.pdf"), width = 10, height = 8)
p_combo_cum2
dev.off()



#Option 3 
p_combo_cum2 <- ggplot(counts_cov, aes(x = factor(OR_cutoff))) +
  # --- Bars for # of combinations (fully opaque now) ---
  geom_col(
    aes(y = n_combos / max(n_combos) * 100, fill = k),
    position = position_dodge(width = 0.9),
    alpha = 1,            # make bars solid
    color = "black"       # optional thin border around bars
  ) +
  geom_text(
    aes(y = n_combos / max(n_combos) * 100, label = n_combos, group = k),
    position = position_dodge(width = 0.9),
    vjust = -0.25, size = 3
  ) +
  
  # --- Black border (outline) for coverage line ---
  geom_line(
    aes(y = pct, group = k),
    color = "black", linewidth = 2.0,           # thick outline
    position = position_dodge(width = 0.9),
    show.legend = FALSE
  ) +
  geom_point(
    aes(y = pct, group = k),
    color = "black", size = 3.5,
    position = position_dodge(width = 0.9),
    show.legend = FALSE
  ) +
  
  # --- Colored coverage line on top ---
  geom_line(
    aes(y = pct, color = k, group = k),
    linewidth = 1.0,
    position = position_dodge(width = 0.9)
  ) +
  geom_point(
    aes(y = pct, color = k, group = k),
    size = 2.0,
    position = position_dodge(width = 0.9)
  ) +
  
  # --- Axis scaling: left = coverage %, right = # combinations ---
  scale_y_continuous(
    name = "% participants who satisfy criteria",
    breaks = seq(0, 100, by = 10),
    sec.axis = sec_axis(
      ~ . * max(counts_cov$n_combos) / 100,
      name = "# of combinations",
      breaks = seq(0, max(counts_cov$n_combos), by = max(counts_cov$n_combos) / 10),
      labels = scales::comma
    )
  ) +
  labs(
    x = expression("OR cutoff (log" [10] * " scale)"),
    fill = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "black", face = "bold"),
    axis.title.y.right = element_text(color = "grey40"),
    axis.text.y.right  = element_text(color = "grey40"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank()
  )

pdf(file.path(path_plots, "LogOR_Cutoff_Final.pdf"), width = 10, height = 8)
p_combo_cum2
dev.off()



# 1) Build participant sets (unique IDs) for each k
ids2_set <- ids2_tbl %>% pull(ids) %>% unlist(use.names = FALSE) %>% unique()
ids3_set <- ids3_tbl %>% pull(ids) %>% unlist(use.names = FALSE) %>% unique()
ids4_set <- ids4_tbl %>% pull(ids) %>% unlist(use.names = FALSE) %>% unique()

sets <- list(
  `k=2` = ids2_set,
  `k=3` = ids3_set,
  `k=4` = ids4_set
)



# 1) Build participant sets from your tables
ids2_set <- ids2_tbl %>% dplyr::pull(ids) %>% unlist(use.names = FALSE) %>% unique()
ids3_set <- ids3_tbl %>% dplyr::pull(ids) %>% unlist(use.names = FALSE) %>% unique()
ids4_set <- ids4_tbl %>% dplyr::pull(ids) %>% unlist(use.names = FALSE) %>% unique()

sets <- list(`k=2` = ids2_set, `k=3` = ids3_set, `k=4` = ids4_set)
cols <- c(`k=2` = "#1f77b4", `k=3` = "#2ca02c", `k=4` = "#d62728")

# 2) Euler diagram
fit <- euler(sets)

pdf(file.path(path_plots, "Overlap_Combiantions_Euler.pdf"), width = 10, height = 10)
plot(
  fit,
  quantities = TRUE,                                # show counts
  fills      = list(fill = unname(cols[names(sets)]), alpha = 0.35),
  edges      = list(col  = unname(cols[names(sets)]), lwd = 2),
  labels     = list(col = "black"),
  main       = "Participants overlap (k=2, k=3, k=4)"
)
dev.off()



cols <- c(`k=2` = "#1f77b4", `k=3` = "#2ca02c", `k=4` = "#d62728")

pdf(file.path(path_plots, "Overlap_Combiantions.pdf"), width = 10, height = 10)
ggVennDiagram(sets,label_alpha  = 0,set_color=cols) + guides(fill = guide_legend(title = "Title")) +
  theme(legend.title = element_text(color = "red"),
        legend.position = "bottom")
dev.off()


# assume your mapping table looks like:
lookup <- select_pheno[,c("Phenotype","Phenotype_Description_v2","Phenotype_Index")]
lookup$dict<-paste0(lookup$Phenotype_Index,"-",lookup$Phenotype_Description_v2)
lookup<-lookup[,c("Phenotype","dict")]



#Find Participants missing in Combionations of 2
ids2_pid<-unique(unlist(ids2_tbl$ids))
all_ids<-row.names(phenotype_adult)
missing_combo2<-phenotype_adult[setdiff(all_ids,ids2_pid),]


# choose phenotype columns (drop ID cols if present)
pheno_cols <- select_pheno$Phenotype[select_pheno$Present_in_GATA2=="Yes" & select_pheno$Present_in_UKB=="Yes" &select_pheno$Keep_Adults=="Yes"]

df_out <- missing_combo2 %>%
  rowwise() %>%
  mutate(
    Phenotypes_Y = {
      vals <- c_across(all_of(pheno_cols))
      hits <- pheno_cols[which(str_to_upper(as.character(vals)) == "Y")]
      if (length(hits) == 0) "no Phenotype" else paste(hits, collapse = "; ")
    },
    num_pheno = {
      vals <- c_across(all_of(pheno_cols))
      hits <- pheno_cols[which(str_to_upper(as.character(vals)) == "Y")]
      length(hits)
    }
  ) %>%
  ungroup()


missing_combo2_pheno<-cbind(row.names(missing_combo2),df_out$Phenotypes_Y,df_out$num_pheno)



# convert to tibble
missing_combo2_pheno <- as_tibble(missing_combo2_pheno)
colnames(missing_combo2_pheno) <- c("Row_ID", "Phenotypes_Y", "num_pheno")

# create labeled version of Phenotypes_Y using lookup$Phenotype and lookup$dict
missing_combo2_pheno <- missing_combo2_pheno %>%
  rowwise() %>%
  mutate(
    Phenotypes_Label = {
      # split phenotype indices
      idxs <- unlist(str_split(Phenotypes_Y, ";\\s*"))
      idxs <- idxs[idxs != "no Phenotype"]
      
      # map index to dictionary term
      descs <- lookup$dict[match(idxs, lookup$Phenotype)]
      descs <- ifelse(is.na(descs), idxs, descs)
      
      # combine index and mapped term
      labels <- descs
      
      if (length(labels) == 0) "no Phenotype" else paste(labels, collapse = "; ")
    }
  ) %>%
  ungroup()





#Find Participants missing in Combionations of 3
ids3_pid<-unique(unlist(ids3_tbl$ids))
all_ids<-row.names(phenotype_adult)
missing_combo3<-phenotype_adult[setdiff(all_ids,ids3_pid),]


# choose phenotype columns (drop ID cols if present)
pheno_cols <- select_pheno$Phenotype[select_pheno$Present_in_GATA2=="Yes" & select_pheno$Present_in_UKB=="Yes" &select_pheno$Keep_Adults=="Yes"]

df_out <- missing_combo3 %>%
  rowwise() %>%
  mutate(
    Phenotypes_Y = {
      vals <- c_across(all_of(pheno_cols))
      hits <- pheno_cols[which(str_to_upper(as.character(vals)) == "Y")]
      if (length(hits) == 0) "no Phenotype" else paste(hits, collapse = "; ")
    },
    num_pheno = {
      vals <- c_across(all_of(pheno_cols))
      hits <- pheno_cols[which(str_to_upper(as.character(vals)) == "Y")]
      length(hits)
    }
  ) %>%
  ungroup()


missing_combo3_pheno<-cbind(row.names(missing_combo3),df_out$Phenotypes_Y,df_out$num_pheno)


# convert to tibble
missing_combo3_pheno <- as_tibble(missing_combo3_pheno)
colnames(missing_combo3_pheno) <- c("Row_ID", "Phenotypes_Y", "num_pheno")

# create labeled version of Phenotypes_Y using lookup$Phenotype and lookup$dict
missing_combo3_pheno <- missing_combo3_pheno %>%
  rowwise() %>%
  mutate(
    Phenotypes_Label = {
      # split phenotype indices
      idxs <- unlist(str_split(Phenotypes_Y, ";\\s*"))
      idxs <- idxs[idxs != "no Phenotype"]
      
      # map index to dictionary term
      descs <- lookup$dict[match(idxs, lookup$Phenotype)]
      descs <- ifelse(is.na(descs), idxs, descs)
      
      # combine index and mapped term
      labels <- descs
      
      if (length(labels) == 0) "no Phenotype" else paste(labels, collapse = "; ")
    }
  ) %>%
  ungroup()




#Find Participants missing in Combionations of 4
ids4_pid<-unique(unlist(ids4_tbl$ids))
all_ids<-row.names(phenotype_adult)
missing_combo4<-phenotype_adult[setdiff(all_ids,ids4_pid),]


# choose phenotype columns (drop ID cols if present)
pheno_cols <- select_pheno$Phenotype[select_pheno$Present_in_GATA2=="Yes" & select_pheno$Present_in_UKB=="Yes" &select_pheno$Keep_Adults=="Yes"]

df_out <- missing_combo4 %>%
  rowwise() %>%
  mutate(
    Phenotypes_Y = {
      vals <- c_across(all_of(pheno_cols))
      hits <- pheno_cols[which(str_to_upper(as.character(vals)) == "Y")]
      if (length(hits) == 0) "no Phenotype" else paste(hits, collapse = "; ")
    },
    num_pheno = {
      vals <- c_across(all_of(pheno_cols))
      hits <- pheno_cols[which(str_to_upper(as.character(vals)) == "Y")]
      length(hits)
    }
  ) %>%
  ungroup()


missing_combo4_pheno<-cbind(row.names(missing_combo4),df_out$Phenotypes_Y,df_out$num_pheno)



# convert to tibble
missing_combo4_pheno <- as_tibble(missing_combo4_pheno)
colnames(missing_combo4_pheno) <- c("Row_ID", "Phenotypes_Y", "num_pheno")

# create labeled version of Phenotypes_Y using lookup$Phenotype and lookup$dict
missing_combo4_pheno <- missing_combo4_pheno %>%
  rowwise() %>%
  mutate(
    Phenotypes_Label = {
      # split phenotype indices
      idxs <- unlist(str_split(Phenotypes_Y, ";\\s*"))
      idxs <- idxs[idxs != "no Phenotype"]
      
      # map index to dictionary term
      descs <- lookup$dict[match(idxs, lookup$Phenotype)]
      descs <- ifelse(is.na(descs), idxs, descs)
      
      # combine index and mapped term
      labels <- descs
      
      if (length(labels) == 0) "no Phenotype" else paste(labels, collapse = "; ")
    }
  ) %>%
  ungroup()










# --- 1) Build mapping vector: P# -> Phenotype_Description_v2 ---
pheno_map <- select_pheno_full %>%
  dplyr::select(Phenotype_Index, Phenotype_Description_v2,Phenotype) %>%
  distinct() %>%
  mutate(Phenotype = as.character(Phenotype),
    Phenotype_Index = as.character(Phenotype_Index),
    Phenotype_Description_v2 = as.character(Phenotype_Description_v2)
  )
map_vec <- rlang::set_names(pheno_map$Phenotype, pheno_map$Phenotype_Index)

# --- 2) Detect the Y-list column in df_out ---
y_col <- "Phenotypes_Y"
# --- 3) Helper: map 'P1; P5; ...' -> 'P1-Desc; P5-Desc; ...' (or 'no Phenotype') ---
format_pheno_list <- function(s) {
  if (is.na(s) || s == "") return("no Phenotype")
  parts <- str_split(s, ";")[[1]] |> str_trim()
  parts <- parts[parts != ""]
  if (length(parts) == 0) return("no Phenotype")
  if (length(parts) == 1 && tolower(parts) %in% c("no phenotype", "none")) return("no Phenotype")
  labs <- map_chr(parts, function(p) {
    desc <- map_vec[[p]]
    if (!is.null(desc) && !is.na(desc) && desc != "") paste0(p, "-", desc) else p
  })
  paste(labs, collapse = "; ")
}

# --- 4) Get an ID column from rownames (fallback to ParticipantID/eid if needed) ---
df_id <- df_out
if (!is.null(rownames(df_out)) && any(nchar(rownames(df_out)) > 0)) {
  df_id <- rownames_to_column(df_out, var = "ID")
} else if ("ParticipantID" %in% names(df_out)) {
  df_id <- df_out %>% rename(ID = ParticipantID)
} else if ("eid" %in% names(df_out)) {
  df_id <- df_out %>% rename(ID = eid)
} else {
  df_id <- df_out %>% mutate(ID = row_number())
}




writexl::write_xlsx(combos_df_with_counts_and_cat[,-8],file.path(path_plots, "Filtered_Overlap_Phenotypes_102125.xlsx"))




ph_counts_2 <- comb2_df %>%
  filter(!is.na(Combination)) %>%
  mutate(Combination = str_squish(Combination)) %>%
  separate_rows(Combination, sep = "_") %>%         # split P1_P2 -> P1, P2
  rename(Phenotype_Index = Combination) %>%
  count(Phenotype_Index, name = "n") %>%
  arrange(desc(n))

ph_counts_2<-merge(ph_counts_2,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2")],by="Phenotype_Index")
ph_counts_2 <- ph_counts_2[order(-ph_counts_2$n), ]
writexl::write_xlsx(ph_counts_2,file.path(path_plots, "Ranked_list_Combinations_2.xlsx"))




ph_counts_3 <- comb3_df %>%
  filter(!is.na(Combination)) %>%
  mutate(Combination = str_squish(Combination)) %>%
  separate_rows(Combination, sep = "_") %>%         # split P1_P2 -> P1, P2
  rename(Phenotype_Index = Combination) %>%
  count(Phenotype_Index, name = "n") %>%
  arrange(desc(n))

ph_counts_3<-merge(ph_counts_3,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2")],by="Phenotype_Index")
ph_counts_3 <- ph_counts_3[order(-ph_counts_3$n), ]
writexl::write_xlsx(ph_counts_3,file.path(path_plots, "Ranked_list_Combinations_3.xlsx"))

ph_counts_4 <- comb4_df %>%
  filter(!is.na(Combination)) %>%
  mutate(Combination = str_squish(Combination)) %>%
  separate_rows(Combination, sep = "_") %>%         # split P1_P2 -> P1, P2
  rename(Phenotype_Index = Combination) %>%
  count(Phenotype_Index, name = "n") %>%
  arrange(desc(n))

ph_counts_4<-merge(ph_counts_4,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2")],by="Phenotype_Index")
ph_counts_4 <- ph_counts_4[order(-ph_counts_4$n), ]
writexl::write_xlsx(ph_counts_4,file.path(path_plots, "Ranked_list_Combinations_4.xlsx"))


combos_exclude_single<-combos_df_with_counts_and_cat[combos_df_with_counts_and_cat$Combination_Number!="Single",]

ph_counts_all<-combos_exclude_single %>%
  filter(!is.na(Combination)) %>%
  mutate(Combination = str_squish(Combination)) %>%
  separate_rows(Combination, sep = "_") %>%         # split P1_P2 -> P1, P2
  rename(Phenotype_Index = Combination) %>%
  count(Phenotype_Index, name = "n") %>%
  arrange(desc(n))

ph_counts_all<-merge(ph_counts_all,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2")],by="Phenotype_Index")
ph_counts_all <- ph_counts_all[order(-ph_counts_all$n), ]
writexl::write_xlsx(ph_counts_all,file.path(path_plots, "Ranked_list_Combinations_All.xlsx"))

writexl::write_xlsx(combos_df_with_counts_and_cat[,-8],file.path(path_plots, "GATA2_Phenotype_Search.xlsx"))



ggplot(ph_counts_all, aes(x = n, y = reorder(Phenotype_Description_v2, n))) +
  geom_col(fill = "#4C78A8", width = 0.7) +
  labs(
    x = "Count (n)",
    y = "Phenotype",
    title = "Phenotype Frequency"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 14)
  )


pheno_missing<-setdiff(select_pheno_combo$Phenotype_Index,ph_counts_all$Phenotype_Index)
pheno_missing<-select_pheno_full[select_pheno_full$Phenotype_Index%in%pheno_missing,]
write.csv(pheno_missing[,c(1,10,12)],file.path(path_plots, "Phenotype_Missing_from_RankedList.csv"))




# --- helper: canonicalize a combo like "P2_P1" -> "P1_P2"
canon <- function(x) paste(sort(x), collapse = "_")

# --- make all 2-subsets (pairs) from each row's Combination column
make_pair_index <- function(df) {
  df %>%
    filter(!is.na(Combination)) %>%
    mutate(Combination = str_replace_all(Combination, "\\s+", ""),
           toks = str_split(Combination, "_")) %>%
    transmute(
      Combination,
      pair = map(toks, ~{
        x <- sort(.x)
        if (length(x) < 2) character(0)
        else apply(combn(x, 2), 2, canon)
      })
    ) %>%
    unnest(pair)
}

# 1) Canonical set of pairs from combos_2 (P1_P2 == P2_P1)
pairs2 <- comb2_df %>%
  filter(!is.na(Combination)) %>%
  transmute(pair = map_chr(
    str_split(str_replace_all(Combination, "\\s+", ""), "_"),
    canon
  )) %>%
  distinct()

# 2) Build pair→combination index for combos_3 and combos_4
pairs3_idx <- make_pair_index(comb3_df)
pairs4_idx <- make_pair_index(comb4_df)

# 3) Counts per pair
counts3 <- pairs3_idx %>% count(pair, name = "n_in_combos_3")
counts4 <- pairs4_idx %>% count(pair, name = "n_in_combos_4")

# 4) (optional) lists of which combos contain each pair
which3 <- pairs3_idx %>%
  group_by(pair) %>%
  summarise(combos_3_having = list(unique(Combination)), .groups = "drop")

which4 <- pairs4_idx %>%
  group_by(pair) %>%
  summarise(combos_4_having = list(unique(Combination)), .groups = "drop")

# 5) Final table: for every pair from combos_2, how often it appears in combos_3/4
pair_overlap <- pairs2 %>%
  left_join(counts3, by = "pair") %>%
  left_join(counts4, by = "pair") %>%
  left_join(which3, by = "pair") %>%
  left_join(which4, by = "pair") %>%
  mutate(
    n_in_combos_3 = coalesce(n_in_combos_3, 0L),
    n_in_combos_4 = coalesce(n_in_combos_4, 0L)
  ) %>%
  arrange(desc(n_in_combos_3), desc(n_in_combos_4))

pair_overlap

# --- helper: ensure we have a numeric logOR column ---
ensure_logOR <- function(df) {
  out <- df
  out
}

# --- 1) From combos_exclude_single, get logOR for exactly 2-phenotype combos (canonicalized) ---
combos_pairs_or <- combos_exclude_single %>%
  ensure_logOR() %>%
  mutate(
    Combination = str_replace_all(Combination, "\\s+", ""),
    toks = str_split(Combination, "_"),
    k = purrr::map_int(toks, length),
    canon = purrr::map_chr(toks, ~ paste(sort(.x), collapse = "_"))
  ) %>%
  filter(k == 2) %>%
  group_by(canon) %>%
  summarise(
    # if duplicates exist for a pair, average logOR (change to first() if you prefer)
    logOR_pair = if (all(is.na(logOR))) NA_real_ else mean(logOR, na.rm = TRUE),
    .groups = "drop"
  )

# --- 2) Add logOR_pair to your pair_overlap by matching pair <-> canon ---
pair_overlap2 <- pair_overlap %>%
  left_join(combos_pairs_or, by = c("pair" = "canon")) %>%
  mutate(OR_pair = exp(logOR_pair))  # optional: OR on original scale

# --- 3) Add phenotype descriptions from select_pheno ---
desc_map <- select_pheno %>%
  distinct(Phenotype_Index, Phenotype_Description_v2)

pair_overlap2 <- pair_overlap2 %>%
  separate(pair, into = c("P_A", "P_B"), sep = "_", remove = FALSE) %>%
  left_join(desc_map, by = c("P_A" = "Phenotype_Index")) %>%
  rename(P_A_desc = Phenotype_Description_v2) %>%
  left_join(desc_map, by = c("P_B" = "Phenotype_Index")) %>%
  rename(P_B_desc = Phenotype_Description_v2) %>%
  mutate(
    pair_label = paste0(
      P_A, "-", coalesce(P_A_desc, "NA"), "; ",
      P_B, "-", coalesce(P_B_desc, "NA")
    )
  )

# Result:
# pair_overlap2 has:
# - logOR_pair (and OR_pair) for the exact pair from combos_exclude_single
# - P_A_desc / P_B_desc and a combined pair_label
pair_overlap2

writexl::write_xlsx(pair_overlap2[,-c(6,7)],file.path(path_plots, "Overlap_Combinations_102125.xlsx"))




# Build boolean membership per pair: does this k=2 pair appear in any k=3 / k=4 combo?
membership_df <- pair_overlap2 %>%
  transmute(
    pair,
    k3 = n_in_combos_3 > 0,
    k4 = n_in_combos_4 > 0
  )

# Quick redundancy summary (optional)
redundancy <- membership_df %>%
  summarise(
    pairs_total = n(),
    pairs_in_any = sum(k3 | k4),
    pct = 100 * pairs_in_any / pairs_total
  )


#
final_gata2_combo2_3<-combos_exclude_single[combos_exclude_single$Combination_Number%in%c("Combinations of 2","Combinations of 3"),]
final_gata2_combo2_3<-final_gata2_combo2_3[final_gata2_combo2_3$logOR>=3,]
writexl::write_xlsx(final_gata2_combo2_3[,-8],file.path(path_plots, "Final_GATA2_Combinations_atORCutoff3_022526.xlsx"))





# Map each combination to its matching IDs, and attach the OR
ids_tbl <- ids_by_combo(phenotype_gata2_adult, final_gata2_combo2_3, "Combination") %>%
  left_join(final_gata2_combo2_3 %>% dplyr::select(Combination, logOR), by = "Combination")


# Extract all IDs from Match_PID, split by ";", flatten, and count unique
unique_pid <- ids_tbl$ids %>% unlist()%>%
  na.omit() %>%                                  # remove NAs
  str_split(";") %>%                             # split into list of vectors
  unlist() %>%                                   # flatten into single vector
  str_trim() %>%                                 # remove any stray spaces
  unique() 

length(unique_pid)/212


final_gata2_combo2_3$pair<-lapply(X = final_gata2_combo2_3$Combination,canon)
final_gata2_combo2_3$pair<-unlist(final_gata2_combo2_3$pair)
ph_counts_final<-final_gata2_combo2_3 %>%
  filter(!is.na(Combination)) %>%
  mutate(Combination = str_squish(Combination)) %>%
  separate_rows(Combination, sep = "_") %>%         # split P1_P2 -> P1, P2
  rename(Phenotype_Index = Combination) %>%
  count(Phenotype_Index, name = "n") %>%
  arrange(desc(n))

writexl::write_xlsx(ph_counts_final,file.path(path_plots, "Final_Phenotypes_60.xlsx"))




high_counts<-final_gata2_combo2_3[final_gata2_combo2_3$Counts_GATA2>=10,]
ph_counts_highcounts<-high_counts %>%
  filter(!is.na(Combination)) %>%
  mutate(Combination = str_squish(Combination)) %>%
  separate_rows(Combination, sep = "_") %>%         # split P1_P2 -> P1, P2
  rename(Phenotype_Index = Combination) %>%
  count(Phenotype_Index, name = "n") %>%
  arrange(desc(n))
ph_counts_highcounts<-merge(ph_counts_highcounts,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2")],by="Phenotype_Index")
ph_counts_highcounts <- ph_counts_highcounts[order(-ph_counts_highcounts$n), ]

ph_counts_final<-merge(ph_counts_final,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2","New_Category")],by="Phenotype_Index")
ph_counts_final <- ph_counts_final[order(-ph_counts_final$n), ]


writexl::write_xlsx(ph_counts_final,file.path(path_plots, "Ranked_list_Combinations_Final_022526.xlsx"))
writexl::write_xlsx(combos_df_filtered,file.path(path_plots, "Final_All_Combinations_022526.xlsx"))

table(combos_df_filtered$Combination_Number)
combos_exclude_single

#Coverage of Final GATA2 Phenotypes


# 1) Build a lookup: phenotype column -> vector of participant IDs (non-NA, unique)
all_pairs <- na.omit(unique(final_gata2_combo2_3$pair))
all_idxs  <- unique(unlist(strsplit(unlist(all_pairs), "_", fixed = TRUE)))
pheno_cols <- intersect(all_idxs, colnames(phenotype_gata2_adult))

id_list <- lapply(pheno_cols, function(col) {
  x <- phenotype_gata2_adult[[col]]
  ids <- unique(na.omit(as.character(x)))
  ids[ids != ""]
})

names(id_list) <- pheno_cols

# 2) Helper: for "P1_P2" return intersected IDs
match_ids_for_pair <- function(pair_str) {
  if (is.na(pair_str) || !nzchar(pair_str)) return(character(0))
  idxs <- strsplit(pair_str, "_", fixed = TRUE)[[1]]
  # require all columns to exist
  if (!all(idxs %in% names(id_list))) return(character(0))
  Reduce(intersect, id_list[idxs])
}

# 3) Compute matches for each row's 'pair'
matches <- map(final_gata2_combo2_3$pair, match_ids_for_pair)

final_gata2_combo2_3 <- final_gata2_combo2_3 %>%
  mutate(
    Match_PID   = ifelse(lengths(matches) > 0,
                         map_chr(matches, ~ paste(.x, collapse = ";")),
                         NA_character_),
    Match_Count = lengths(matches)
  )

# (optional) quick peek
final_gata2_combo2_3 %>%
  dplyr::select(Combination, pair, Match_Count, Match_PID) %>%
  head()



# Extract all IDs from Match_PID, split by ";", flatten, and count unique
unique_pid <- final_gata2_combo2_3$Match_PID %>%
  na.omit() %>%                                  # remove NAs
  str_split(";") %>%                             # split into list of vectors
  unlist() %>%                                   # flatten into single vector
  str_trim() %>%                                 # remove any stray spaces
  unique()                                       # keep unique IDs



# Count unique participants
length(unique_pid)
length(unique_pid)/N_cases



#Check which phenotypes are each log10 odds ratio
library(ComplexHeatmap)
library(circlize)


# ---- 1. Prepare combos with valid log_OR ----
df <- combos_df_filtered %>%
  mutate(logOR=logOR) %>%
  filter(!is.na(logOR)) %>%
  distinct(Combination, Combination_Number, logOR)

# ---- 2. Expand into phenotype indices ----
long <- df %>%
  transmute(
    Combination,
    logOR,
    Phenotype_Index = str_split(Combination, "_")
  ) %>%
  unnest_longer(Phenotype_Index)

# ---- 3. Define bins and set factor order explicitly (ascending) ----
breaks <- c(-1,0, 1, 2, 3, 4, 5, 6)
labels <- as.character(breaks)

long <- long %>%
  mutate(
    logOR_bin = cut(
      logOR,
      breaks = c(-Inf, breaks[-1], Inf),
      labels = labels,
      right = TRUE
    )
  )

# enforce factor order ascending
long$logOR_bin <- factor(long$logOR_bin, levels = labels, ordered = TRUE)

# ---- 4. Build binary matrix ----
mat_df <- long %>%
  mutate(value = 1L) %>%
  group_by(Phenotype_Index, logOR_bin) %>%
  summarise(value = 1L, .groups = "drop") %>%
  pivot_wider(
    names_from = logOR_bin,
    values_from = value,
    values_fill = 0L
  )

mat <- mat_df %>%
  column_to_rownames("Phenotype_Index") %>%
  as.matrix()

# enforce ordered factor levels
long$logOR_bin <- factor(long$logOR_bin, levels = labels, ordered = TRUE)

long <- long %>%
  left_join(lookup_cat_merge %>% dplyr::select(Phenotype_Index, dict), by = "Phenotype_Index")

# ---- 4. Build binary matrix ----
mat_df <- long %>%
  mutate(value = 1L) %>%
  group_by(dict, logOR_bin) %>%
  summarise(value = 1L, .groups = "drop") %>%
  pivot_wider(
    names_from = logOR_bin,
    values_from = value,
    values_fill = 0L
  )

mat <- mat_df %>%
  column_to_rownames("dict") %>%
  as.matrix()

# ---- FIX: keep only bins that actually exist, but order them ascending ----
existing_bins <- intersect(labels, colnames(mat))
mat <- mat[, existing_bins, drop = FALSE]

# ---- 5. Plot ----
ht_opt$message = FALSE
Heatmap(
  mat,
  name = "Present",
  col = c("0" = "white", "1" = "black"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,   # keep fixed order
  column_title = "log10(OR) bins (ascending)",
  row_title = "Phenotypes",
  row_names_gp = grid::gpar(fontsize = 8),
  use_raster = TRUE
)

# ---- 5. Heatmap ----
ht_opt$message = FALSE
Heatmap(
  mat,
  name = "Present",
  col = c("0" = "white", "1" = "black"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,   # keep logOR order fixed
  column_title = "log10(OR) bins (ascending)",
  row_title = "Phenotypes",
  row_names_gp = grid::gpar(fontsize = 8),
  use_raster = TRUE
)


#Heatmap with Columns as Participant IDs

combos_df_filtered$pair <- sapply(combos_df_filtered$Combination, canon)


# 1) Build a lookup: phenotype column -> vector of participant IDs (non-NA, unique)
all_pairs <- na.omit(unique(combos_df_filtered$pair))
all_idxs  <- unique(unlist(strsplit(all_pairs, "_", fixed = TRUE)))
pheno_cols <- intersect(all_idxs, colnames(phenotype_gata2_adult))

id_list <- lapply(pheno_cols, function(col) {
  x <- phenotype_gata2_adult[[col]]
  ids <- unique(na.omit(as.character(x)))
  ids[ids != ""]
})
names(id_list) <- pheno_cols

# 2) Helper: for "P1_P2" return intersected IDs
match_ids_for_pair <- function(pair_str) {
  if (is.na(pair_str) || !nzchar(pair_str)) return(character(0))
  idxs <- strsplit(pair_str, "_", fixed = TRUE)[[1]]
  # require all columns to exist
  if (!all(idxs %in% names(id_list))) return(character(0))
  Reduce(intersect, id_list[idxs])
}

# 3) Compute matches for each row's 'pair'
matches <- map(combos_df_filtered$pair, match_ids_for_pair)

combos_df_filtered <- combos_df_filtered %>%
  mutate(
    Match_PID   = ifelse(lengths(matches) > 0,
                         map_chr(matches, ~ paste(.x, collapse = ";")),
                         NA_character_),
    Match_Count = lengths(matches)
  )

# (optional) quick peek
combos_df_filtered %>%
  dplyr::select(Combination, pair, Match_Count, Match_PID) %>%
  head()




library(ComplexHeatmap)
library(circlize)

# -----------------------------
# 0) Define log10(OR) bins (ascending)
# -----------------------------
breaks <- c(-1,0, 1, 2, 3, 4, 5, 6,7)
labels  <- as.character(breaks)

# -----------------------------
# 1) Use combos_df_filtered with Match_PID
#    Ensure log_OR exists and compute logOR_bin
# -----------------------------
df <- combos_df_filtered %>%
  mutate(
    logOR = dplyr::coalesce(logOR, if_else(!is.na(OR) & OR > 0, log10(OR), NA_real_))
  ) %>%
  filter(!is.na(logOR)) %>%
  mutate(
    # normalize Match_PID to character (might be NA), and split
    Match_PID = ifelse(is.na(Match_PID), NA_character_, as.character(Match_PID)),
    logOR_bin = cut(
      logOR,
      breaks = c(-Inf, breaks[-1], Inf),
      labels = labels,
      right  = TRUE
    )
  )

df$logOR_bin <- factor(df$logOR_bin, levels = labels, ordered = TRUE)

# helper: split "A;B;C" into vector
split_ids <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  out <- str_split(x, ";", simplify = FALSE)[[1]]
  out <- str_trim(out)
  unique(out[out != ""])
}

# -----------------------------
# 2) Participant → lowest logOR bin (via Match_PID membership)
# -----------------------------
pid_min_bin <- df %>%
  mutate(PID = lapply(Match_PID, split_ids)) %>%
  filter(lengths(PID) > 0) %>%
  dplyr::select(logOR_bin, PID) %>%
  unnest_longer(PID) %>%
  group_by(PID) %>%
  summarise(min_bin_idx = max(match(logOR_bin, labels),na.rm=TRUE ),.groups = "drop")

# -----------------------------
# 3) All participants & column order
#     (include anyone present in phenotype table even if not in combos)
# -----------------------------
# Build phenotype→participant sets from phenotype_gata2_adult (IDs in cells)
# Keep only phenotypes that appear anywhere in combos_df_filtered
phenos_in_combos <- df$Combination %>%
  str_split("_") %>% unlist() %>% unique()

pheno_cols <- intersect(phenos_in_combos, colnames(phenotype_gata2_adult))

id_list <- lapply(pheno_cols, function(col) {
  ids <- unique(na.omit(as.character(phenotype_gata2_adult[[col]])))
  ids[ids != ""]
})
names(id_list) <- pheno_cols

pids_from_pheno <- sort(unique(unlist(id_list)))
pids_from_match <- sort(unique(unlist(lapply(df$Match_PID, split_ids))))

all_pids <- sort(unique(c(pids_from_pheno, pids_from_match)))

# align min-bin index to all_pids (NA for participants not in any matched combo)
min_idx_vec <- setNames(pid_min_bin$min_bin_idx, pid_min_bin$PID)[all_pids]

# order: participants with defined min bin first (ascending), then the rest by ID
ord <- order(is.na(min_idx_vec), min_idx_vec, all_pids, na.last = TRUE)
pid_order <- all_pids[ord]

# top annotation labels (bins) aligned to pid_order
pid_bin_label <- labels[min_idx_vec]
names(pid_bin_label) <- all_pids
pid_bin_label <- pid_bin_label[pid_order]


# -----------------------------
# 4) Build phenotype × participant binary matrix
#     1 if participant has that phenotype (regardless of combo)
# -----------------------------
mat_df <- tibble(Phenotype_Index = pheno_cols) %>%
  mutate(tmp = Phenotype_Index) %>%
  rowwise() %>%
  mutate(PID = list(id_list[[tmp]])) %>%
  ungroup() %>%
  dplyr::select(-tmp) %>%
  unnest_longer(PID, values_to = "ParticipantID") %>%
  mutate(value = 1L) %>%
  pivot_wider(names_from = ParticipantID, values_from = value, values_fill = 0L)
mat_df<-merge(lookup_cat_merge[,c("Phenotype_Index","dict")],mat_df)
mat_df[,-1]
mat <- mat_df %>%
  column_to_rownames("dict") %>%
  as.matrix()

# keep only columns we care about, ordered by pid_order
existing_cols <- intersect(pid_order, colnames(mat))
mat <- mat[, existing_cols, drop = FALSE]


#5) Plot heatmap
#color map for bins present in the top annotation
bins_present <- unique(na.omit(pid_bin_label))
bin_pal <- setNames(
  colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(max(length(bins_present), 3)),
  bins_present
)

ha <- HeatmapAnnotation(
  `log10(OR) bin` = factor(pid_bin_label, levels = labels, ordered = TRUE),
  col = list(`log10(OR) bin` = bin_pal),
  annotation_name_side = "right"
)

ht_opt$message = FALSE

library(ComplexHeatmap)
library(InteractiveComplexHeatmap)

ht <- Heatmap(
  mat,
  name = "Present",
  col = c("0" = "white", "1" = "black"),
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # keep ParticipantID order fixed
  column_title = "Participants (ordered by log10(OR) bin from Match_PID)",
  row_title = "Phenotypes",
  row_names_gp = grid::gpar(fontsize = 8),
  use_raster = TRUE
)

ht <- draw(ht)
htShiny(ht)   # launches an interactive viewer


## For spot check-



Combinations_final_check<-merge(final_gata2_combo2_3,combined_df_AllVar,by="Combination")
Combinations_final_check<-Combinations_final_check[Combinations_final_check$logOR.x>=3,]
Combinations_final_check<-Combinations_final_check[Combinations_final_check$Combination_Number.x!="Single",]
Keep_col<-c("Combinations_Labels","Counts_GATA2.x","Counts_UKB.x","OR.y","logOR.y","OR_Lower_CI","OR_Upper_CI")
Combinations_final_check<-Combinations_final_check[,Keep_col]
Combinations_final_check$CI<-paste0(round(log10(Combinations_final_check$OR_Lower_CI),2),",",round(log10(Combinations_final_check$OR_Upper_CI),2))
Combinations_final_check$logOR.y<-round(Combinations_final_check$logOR.y,2)
colnames(Combinations_final_check)<-c("Combinations_Labels","Counts_GATA2","Counts_UKB","OR","logOR","OR_Lower_CI","OR_Upper_CI","CI")
)



writexl::write_xlsx(Combinations_final_check,file.path(path_plots, "Final_Combinations_withCI.xlsx"))
###


cutoffs <- -1:6

# ---- Prepare data (robust ID parsing) ----
df <- combos_df_filtered %>%
  mutate(
    log10OR = logOR,
    ids = str_split(Match_PID %||% "", ";\\s*")
  ) %>%
  mutate(ids = lapply(ids, function(v) unique(v[nzchar(v)])))  # drop "" and dupes

# >>> IMPORTANT: set your true participant denominator <<<
N_total <- 212   # <- use your known total; if unknown: length(unique(unlist(df$ids)))

# ---- Stacked counts per cutoff × k ----
bar_df <- lapply(cutoffs, function(cut) {
  df %>%
    filter(Combination_Number %in% c("Combinations of 2","Combinations of 3"), log10OR >= cut) %>%
    count(cutoff = cut, k = Combination_Number, name = "n")
}) %>% bind_rows() %>%
  complete(cutoff, k = c("Combinations of 2","Combinations of 3"), fill = list(n = 0)) %>%
  mutate(k = factor(k, levels = c("Combinations of 3","Combinations of 2"),
                    labels = c("Combinations of 3 phenotypes","Combinations of 2 phenotypes")))

totals_df <- bar_df %>%
  group_by(cutoff) %>%
  summarise(total = sum(n), .groups = "drop")

max_bar_height <- max(totals_df$total)

# ---- Coverage (union of IDs hit by k=2 OR k=3) ----
cov_df <- lapply(cutoffs, function(cut) {
  tmp <- df %>% filter(Combination_Number %in% c("Combinations of 2","Combinations of 3"), log10OR >= cut)
  covered_ids <- unique(unlist(tmp$ids))
  cov_pct <- 100 * length(covered_ids) / N_total
  tibble(
    cutoff = cut,
    cover_ids<-length(covered_ids),
    coverage_pct = cov_pct,
    y_plot = cov_pct * max_bar_height / 100
  )
}) %>% bind_rows()


pdf(file.path(path_plots, "Log_Ratio_Cutoff_Coverage_Final_WidthChanged2.pdf"), width = 6, height = 3)
ggplot() +
  geom_col(
    data = bar_df,
    aes(x = factor(cutoff), y = n / max_bar_height * 100, fill = k),
    width = 0.75
  ) +
  geom_text(
    data = bar_df,
    aes(x = factor(cutoff), y = n / max_bar_height * 100, label = n, group = k),
    position = position_stack(vjust = 0.35),
    size = 3, color = "black"
  ) +
  geom_text(
    data = totals_df,
    aes(x = factor(cutoff), y = total / max_bar_height * 100, label = total),
    vjust = -0.4, size = 3.2, fontface = "bold"
  ) +
  geom_line(
    data = cov_df,
    aes(x = factor(cutoff), y = coverage_pct, group = 1),
    linewidth = 1
  ) +
  geom_point(
    data = cov_df,
    aes(x = factor(cutoff), y = coverage_pct),
    size = 2
  ) +
  scale_fill_manual(values = c("#00BA38", "#F8766D")) +
  scale_y_continuous(
    name = "% participants who satisfy criteria",
    breaks = seq(0, 100, by = 10),
    sec.axis = sec_axis(~ . * max_bar_height / 100,
                        name = "# of combinations",
                        labels = scales::comma)
  ) +
  expand_limits(y = 110) +
  labs(
    x = expression("OR cutoff (log" [10] * " scale)"),
    fill = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),   # remove major grid lines
    panel.grid.minor = element_blank(),   # remove minor grid lines
    panel.background = element_blank(), 
    legend.position = c(1, 0.8),       # ⬅️ move legend inside (x=85%, y=85%)
    legend.justification = c(1, 0.8),  # remove background shading
    axis.title.y.left = element_text(color = "black", size = 10),
    axis.title.y.right = element_text(color = "black", size = 10)
  )

dev.off()


###Ranked List of Phenotypes

ggplot(ph_counts_final, aes(x = reorder(Phenotype_Description_v2, -n), y = n)) +
  geom_col(fill = "steelblue") +   # steel blue bars
  geom_text(aes(label = n),vjust=-0.7,hjust=-0.1, size = 3,angle=45) +  # show counts above bars
  labs(
    x = "Phenotype",
    y = "Number of times in all the selected Combinations",
    title = "Phenotype Frequency in Selected Combinations"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # tilt labels
    panel.grid = element_blank(),               # clean up grid
  )

ggsave(file.path(path_plots, "Ranked_List_of_phenotypes_60.pdf"), width = 15, height = 10)

# Recursive Phenotype ----

#Create a matrix 


#1) Helper: columns-of-IDs  ->  binary matrix per participant
# df: wide, columns are phenotypes (P1..P103 etc), cells are participant IDs (or NA)
# label: "Case" or "Control"

ids_to_binary <- function(df, label) {
  pheno_cols <- colnames(df)
  
  long <- df %>%
    pivot_longer(
      cols = all_of(pheno_cols),
      names_to = "Phenotype_Index",
      values_to = "ParticipantID"
    ) %>%
    filter(!is.na(ParticipantID), ParticipantID != "") %>%
    mutate(ParticipantID = as.character(ParticipantID)) %>%
    distinct(ParticipantID, Phenotype_Index) %>%   # dedupe if same ID listed twice
    mutate(value = 1L)
  
  wide <- long %>%
    pivot_wider(
      id_cols = ParticipantID,
      names_from = Phenotype_Index,
      values_from = value,
      values_fill = 0
    ) %>%
    arrange(ParticipantID) %>%
    mutate(Label = label, ID = ParticipantID, .before = 1)
  
  # ensure all original phenotype columns exist (even if empty)
  missing <- setdiff(pheno_cols, colnames(wide))
  if (length(missing) > 0) wide[missing] <- 0L
  
  # reorder: ID, Label, phenotypes
  wide %>% select(ID, Label, all_of(sort(pheno_cols)))
}

# 2) Build 0/1 matrices for cases and controls
# EXPECTS: phenotype_gata2_adult and phenotype_ukb in your workspace,
# each shaped like your screenshot (columns=P*, cells = EIDs or NA).

case_bin <- ids_to_binary(phenotype_gata2_adult, "Case")
setdiff(row.names(phenotype_gata2_adult_binary),case_bin$ID)


phenotype_ukb <- phenotype_ukb %>% dplyr::mutate(across(everything(), as.character))
ctrl_bin <- ids_to_binary(phenotype_ukb, "Control")

# CSV with a column named 'eid'
eid_vec <- readr::read_csv("Control_GATA2_EID_Final.csv", show_col_types = FALSE) |>
  dplyr::pull(eid) |> as.character() |> unique()


# If ctrl_bin_all doesn't exist yet, start it from ctrl_bin
if (!exists("ctrl_bin_all")) ctrl_bin_all <- ctrl_bin


# 2) Find EIDs that are not already in ctrl_bin
have_ids  <- as.character(ctrl_bin$ID)
new_eids  <- setdiff(eid_vec, have_ids)

# Nothing to add? stop here.
if (length(new_eids) == 0L) {
  message("No new EIDs to add; all are already present.")
} else {
  # 3) Build zero rows for all phenotype columns
  pheno_cols <- setdiff(colnames(ctrl_bin), c("ID", "Label"))
  zero_mat   <- matrix(0L, nrow = length(new_eids), ncol = length(pheno_cols),
                       dimnames = list(NULL, pheno_cols))
  
  new_rows <- tibble(ID = new_eids, Label = "Control") |>
    bind_cols(as_tibble(zero_mat))
  
  # 4) Append to ctrl_bin_all, keep column order, and (optionally) sort
  ctrl_bin_all <- bind_rows(ctrl_bin_all, new_rows) |>
    arrange(ID)
  
  message(sprintf("Added %d new controls. ctrl_bin_all now has %d rows and %d columns.",
                  length(new_eids), nrow(ctrl_bin_all), ncol(ctrl_bin_all)))
}


final_control<-read.csv("Control_GATA2_EID_Final.csv")
ctrl_bin_final<-ctrl_bin_all[ctrl_bin_all$ID%in%final_control$eid,]

# Keep only phenotypes present in BOTH (for a clean common design)
common_pheno <- intersect(
  setdiff(colnames(case_bin), c("ID","Label")),
  setdiff(colnames(ctrl_bin_final), c("ID","Label"))
)




#Add missing case id

# Helper: add missing IDs (as all-zero rows) to a binary matrix
add_missing_ids <- function(bin_df, all_ids, label = c("Case","Control")) {
  label <- match.arg(label)
  all_ids <- as.character(unique(all_ids[!is.na(all_ids) & all_ids != ""]))
  
  have_ids  <- as.character(bin_df$ID)
  new_ids   <- setdiff(all_ids, have_ids)
  if (length(new_ids) == 0L) return(bin_df)
  
  pheno_cols <- setdiff(colnames(bin_df), c("ID", "Label"))
  zero_mat   <- matrix(0L, nrow = length(new_ids), ncol = length(pheno_cols),
                       dimnames = list(NULL, pheno_cols))
  
  new_rows <- tibble(ID = new_ids, Label = label) %>% bind_cols(as_tibble(zero_mat))
  bind_rows(bin_df, new_rows) %>% arrange(ID)
}

# 1) You already have:
# case_bin <- ids_to_binary(phenotype_gata2_adult, "Case")

# 2) Provide a master list of all case EIDs (one per line, or from a data frame)
# Example: read from a text file 'Cases_for_GATA2_eid' (one EID per line)
all_case_ids <-row.names(phenotype_adult)

# 3) Add the missing IDs as all-zero rows
case_bin_all <- add_missing_ids(case_bin, all_case_ids, label = "Case")

model_base <- bind_rows(
  case_bin_all %>% select(ID, Label, all_of(common_pheno)),
  ctrl_bin_final %>% select(ID, Label, all_of(common_pheno))
)

model_base_og<-model_base
model_base<-model_base_og
#model_base<-model_base[, !grepl("P38", names(model_base))]



#R Part Packages 
library(rpart)
library(rpart.plot)
library(pROC)

# model_base: rows = ParticipantID, columns = P1..P103 (0/1), and Label (Case/Control or 1/0)
# Ensure Label is a factor with 2 levels
model_base<-model_base[,-1]
model_base$Label <- factor(model_base$Label, levels = c("Control","Case"))  # adjust if 0/1
pheno_map <- setNames(
  select_pheno_107$Phenotype_Description,
  select_pheno_107$Phenotype_Index
)
old_names <- colnames(model_base)

new_names <- ifelse(
  old_names %in% names(pheno_map),
  pheno_map[old_names],
  old_names
)

colnames(model_base) <- new_names

model_base[["MDS|Pancytopenia"]] <-
  as.integer(
    (model_base$MDS == 1) |
      (model_base$Pancytopenia == 1) 
  )

model_base <- model_base[, !names(model_base) %in% c("MDS","Pancytopenia")]
# ── 1) Train/Test Split ────────────────────────────────────────────────────────
set.seed(42)
#idx  <- sample(seq_len(nrow(model_base)), size = floor(0.8*nrow(model_base)))
train <- model_base
test  <- model_base

# ── 2) Fit rpart (with x-val for cp selection) ────────────────────────────────
fit0 <- rpart(
  Label ~ .,
  data   = train,
  method = "class",
  control = rpart.control(cp = 0.001, minsplit = 5,minbucket=2, maxdepth = 30, xval = 0),
  parms   = list(split = "gini")  # or "information"
)

# ── 3) Pick best cp and prune ─────────────────────────────────────────────────
#View(fit0$cptable)                        # inspect table
bestcp <- fit0$cptable[which.min(fit0$cptable[,"rel error"]), "CP"]
fit    <- prune(fit0, cp = bestcp)

# ── 4) Visualize tree ─────────────────────────────────────────────────────────
pdf(file.path(path_plots, "Decision_Tree_Final_Label_Depth8_.pdf"), width = 8, height = 8)
rpart.plot(fit, type = 1, extra = 1, under = TRUE, fallen.leaves = TRUE,
           main = "Phenotype → Case/Control (rpart)",cex = 0.5)
dev.off()



# ── 5) Variable importance (quick feature readout) ────────────────────────────
imp_variable<-as.data.frame(sort(fit$variable.importance, decreasing = TRUE))  # top 20
imp_variable<-rownames_to_column(imp_variable,"Phenotype")
colnames(imp_variable)[2]<-"Phenotype_Importance"
#imp_variable<-merge(imp_variable,select_pheno[,c("Phenotype_Index","Phenotype_Description_v2")])
imp_variable <- imp_variable[order(-imp_variable$Phenotype_Importance), ]


importance_df <- data.frame(
  Phenotype_Description = names(fit$variable.importance),
  Importance = fit$variable.importance
)
importance_df<-merge(importance_df,select_pheno_combo[,c("Phenotype_Index","Phenotype_Description")])
pdf(file.path(path_plots, "BarPlot_Rpartion.pdf"), width = 8, height = 6)
ggplot(imp_variable, aes(x = reorder(Phenotype, Phenotype_Importance), y = Phenotype_Importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  
  labs(x = "Phenotype", y = "Importance Score",
       title = "Variable Importance in Decision Tree") +
  theme_minimal()+
  theme(panel.grid=element_blank())

dev.off()


#Venn Diagram
set1 <- unique(ph_counts_final$Phenotype_Index)
set2 <- unique(importance_df$Phenotype_Index)

venn_list <- list(
  "ph_count_final" = set1,
  "imp_variable"   = set2
)

ggvenn(
  venn_list,
  fill_color = c("#4DAF4A", "#377EB8"),
  stroke_size = 0.5,
  set_name_size = 4
)


venn.diagram(
  x = list(
    ph_count_final = set1,
    imp_variable   = set2
  ),
  filename = NULL,
  fill = c("darkgreen", "steelblue"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.2
)

length(intersect(set1, set2))   # overlap
length(setdiff(set1, set2))     # only in ph_count_final
setdiff(set2, set1)     # only in imp_variable



# 1) Unique phenotype sets
A <- unique(na.omit(ph_counts_final$Phenotype_Index))
B <- unique(na.omit(importance_df$Phenotype_Index))

# 2) Fit Euler diagram (areas proportional to counts)
fit <- euler(list(
  Combination_Analysis = A,
  R.Part   = B
))

# 3) Plot
plot(
  fit,
  fills = c("#4DAF4A", "#377EB8"),
  edges = TRUE,
  quantities = TRUE,
  labels = TRUE
)


low_freq_pheno <- counts %>%
  filter(Counts_GATA2 < 5) %>%
  pull(Phenotype_Index) %>%
  as.character() %>%
  trimws()

enriched_combos <- Combinations_final_check %>%
  filter(
    logOR > 2,
    Counts_GATA2 > 5
  )

combo_low_freq_any <- enriched_combos %>%
  rowwise() %>%
  mutate(
    pheno_ids = list(
      strsplit(Combinations_Labels, ";")[[1]] %>%
        trimws() %>%
        sub("-.*", "", .)
    ),
    any_low_freq = any(unlist(pheno_ids) %in% low_freq_pheno)
  ) %>%
  ungroup() %>%
  filter(any_low_freq)

  
# your 20 phenotype indices (character vector)
# e.g. phenos20 <- c("P38","P37","P67", ... )
# make sure this exists:

# 1) Extract the phenotype indices from your importance table
phenos_use <- c(importance_df$Phenotype_Index,"P38","P32") 

# 1) keep only the 20 phenotypes that actually exist in the data
cols_use <- intersect(phenos_use, colnames(phenotype_gata2_adult))

# 2) gather all participant IDs across those columns and dedupe
unique_ids <- phenotype_gata2_adult %>%
  select(all_of(cols_use)) %>%
  pivot_longer(everything(), names_to = "Phenotype_Index", values_to = "ParticipantID") %>%
  filter(!is.na(ParticipantID), ParticipantID != "") %>%
  mutate(ParticipantID = as.character(ParticipantID)) %>%
  distinct(ParticipantID) %>%
  pull(ParticipantID)

# 3) how many people have ≥1 of the 20 phenotypes?
n_with_any_of_20 <- length(unique_ids)
n_with_any_of_20/212


# ── 6) Evaluate on test set ───────────────────────────────────────────────────
# Class predictions
pred_class <- predict(fit, newdata = test, type = "class")
tab <- table(Truth = test$Label, Pred = pred_class); tab

acc  <- sum(diag(tab)) / sum(tab)
acc

# Probabilities for ROC/AUC
pred_prob <- predict(fit, newdata = test, type = "prob")[, "Case"]
roc_obj   <- roc(response = test$Label, predictor = pred_prob, levels = c("Control","Case"))
auc(roc_obj)

# Optional: Youden-optimal threshold (instead of 0.5)
coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"))




#GATA2 Genotypes ----
gata2_survey_data$Consequence_Class <- dplyr::case_when(
  
  gata2_survey_data$VariantType %in% c(
    "Nonsense_Mutation",
    "Frameshift_Deletion",
    "Frameshift_Insertion",
    "Large_Deletion"
  ) ~ "LOF",
  
  gata2_survey_data$VariantType %in% c(
    "Missense_Mutation",
    "Inframe_Deletion",
    "Inframe_Insertion"
  ) ~ "Missense",
  
  gata2_survey_data$VariantType %in% c(
    "Noncoding_Mutation"
  ) ~ "Noncoding",
  
  gata2_survey_data$VariantType %in% c(
    "Synonymous"
  ) ~ "Synonymous",
  
  TRUE ~ "Other"
)


library(dplyr)

patient_genotype <- gata2_survey_data %>%
  group_by(prevalence_patient_deidentifier_lab_family_patient) %>%
  summarise(
    Age = first(prevalence_age), 
    Sex= first(prevalence_gender),
    has_LOF        = any(Consequence_Class == "LOF"),
    has_Missense   = any(Consequence_Class == "Missense"),
    has_Noncoding  = any(Consequence_Class == "Noncoding"),
    n_variants     = n()
  )




#GATA2 Zinc finger boundaries ----
ZF1 <- c(294, 344)  # ZF1 aa 294–344
ZF2 <- c(349, 398)  # ZF2 aa 349–398

#Identify CNV/structural descriptions so we don't parse exon numbers as AA positions
is_structural <- function(x) {
  if (is.na(x)) return(FALSE)
  s <- str_to_lower(str_trim(x))
  str_detect(s, "exon\\s*\\d") || str_detect(s, "large\\s*deletion") || str_detect(s, "deletion") && !str_detect(s, "^p\\.")
}

# ---- Extract AA positions from strings like:
# p.Arg396Trp, p.Ala341=, p.Arg362_Asn365del, p.Ala345_Gly346ins..., p.Thr358Asn/p.Leu359Val,
# and delta-style "∆340-381" or "Δ340-381"
extract_aa_positions <- function(x) {
  if (is.na(x)) return(integer(0))
  s <- str_trim(x)
  
  if (s %in% c("", "na", "NA", "?", "n/a")) return(integer(0))
  if (is_structural(s)) return(integer(0))  # handle separately as CNV/structural
  
  # unify delta symbols and split multi-annotations
  s <- str_replace_all(s, "∆", "Δ")
  parts <- unlist(str_split(s, "\\s*/\\s*"))
  
  pos <- integer(0)
  for (p in parts) {
    p <- str_trim(p)
    
    # Δ340-381 style
    if (str_detect(p, "^Δ\\s*\\d+\\s*-\\s*\\d+")) {
      nums <- as.integer(str_extract_all(p, "\\d+")[[1]])
      pos <- c(pos, nums)
      next
    }
    
    # protein HGVS-like strings: grab all numbers (works for del/ins/fs too)
    # examples: p.Arg362_Asn365del -> 362,365 ; p.Ala103GlnfsTer16 -> 103,16 (we'll use min/max)
    if (str_detect(p, "^(p\\.|p\\()")) {
      nums <- as.integer(str_extract_all(p, "\\d+")[[1]])
      pos <- c(pos, nums)
    }
  }
  
  # remove obvious "Ter" artifact numbers when they are much smaller than AA position:
  # e.g., p.Arg384GlyfsTer3 gives 384 and 3; we only need the AA region.
  # Keep numbers >= 20 as a simple heuristic (you can relax/tighten).
  pos <- pos[pos >= 20]
  
  unique(pos)
}

# ---- Assign domain based on overlap with ZF1/ZF2 intervals ----
assign_zf_domain <- function(pos_vec, original_string = NA_character_) {
  if (!is.na(original_string) && is_structural(original_string)) return("CNV_or_Structural")
  if (length(pos_vec) == 0) return("Unknown")
  
  lo <- min(pos_vec, na.rm = TRUE)
  hi <- max(pos_vec, na.rm = TRUE)
  
  overlap <- function(a1, a2, b1, b2) !(a2 < b1 || b2 < a1)
  
  in_zf1 <- overlap(lo, hi, ZF1[1], ZF1[2])
  in_zf2 <- overlap(lo, hi, ZF2[1], ZF2[2])
  
  if (in_zf1 && in_zf2) "ZF1_and_ZF2"
  else if (in_zf1) "ZF1"
  else if (in_zf2) "ZF2"
  else "Non_ZF"
}


# change this to your actual protein-change column name:
protein_col <- "prevalence_protein_change"   # e.g. "HGVSp" or whatever you're using

gata2_survey_data <- gata2_survey_data %>%
  mutate(
    aa_pos_vec = map(.data[[protein_col]], extract_aa_positions),
    aa_pos_min = map_int(aa_pos_vec, ~ if (length(.x)==0) NA_integer_ else min(.x)),
    aa_pos_max = map_int(aa_pos_vec, ~ if (length(.x)==0) NA_integer_ else max(.x)),
    ZF_domain  = map2_chr(aa_pos_vec, .data[[protein_col]], assign_zf_domain)
  )

table(gata2_survey_data$ZF_domain, useNA = "ifany")

  patient_zf <- gata2_survey_data %>%
  group_by(prevalence_patient_deidentifier_lab_family_patient) %>%
  summarise(
    Age = first(prevalence_age), 
    Sex= first(prevalence_gender),
    has_ZF1 = any(ZF_domain == "ZF1"),
    has_ZF2 = any(ZF_domain == "ZF2"),
    has_any_ZF = any(ZF_domain %in% c("ZF1","ZF2","ZF1_and_ZF2")),
    has_CNV = any(ZF_domain == "CNV_or_Structural"),
    n_variants = n(),
    .groups = "drop"
  )



library()

# Inputs:
# combos_exclude_single: has column "Combination" like "P36_P37_P34"
# case_bin_all: 212 rows, columns: .ID (or ID) + P1..P103 (0/1)

combo_col <- "Combination"  # change if needed
id_col    <- ".ID"          # change if needed

DTc <- as.data.table(case_bin_all)
DTx <- as.data.table(combos_exclude_single)

# phenotype columns available in case_bin_all
pheno_cols <- grep("^P\\d+$", names(DTc), value = TRUE)

# parse combos into list of phenotype names
combo_list <- DTx[[combo_col]]
combo_parts <- strsplit(as.character(combo_list), "_", fixed = TRUE)

# pre-allocate output matrix (212 x 4981)
out_mat <- matrix(0L, nrow = nrow(DTc), ncol = length(combo_list))
colnames(out_mat) <- combo_list

# for quick column lookup
col_index <- setNames(seq_along(pheno_cols), pheno_cols)

# numeric phenotype matrix for fast row-wise ops
X <- as.matrix(DTc[, ..pheno_cols])
storage.mode(X) <- "integer"

missing_combos <- list()

for (j in seq_along(combo_parts)) {
  parts <- combo_parts[[j]]
  
  # check if all phenotype columns exist
  miss <- parts[!parts %in% pheno_cols]
  if (length(miss) > 0) {
    missing_combos[[length(missing_combos) + 1]] <- list(combo = combo_list[j], missing = miss)
    next
  }
  
  idx <- col_index[parts]
  # combo present if all required phenotypes == 1  (row-wise min == AND for 0/1)
  out_mat[, j] <- apply(X[, idx, drop = FALSE], 1, min)
}

combo_df <- data.frame(
  DTc[[id_col]],
  out_mat,
  check.names = FALSE
)
names(combo_df)[1] <- id_col

dim(combo_df)  # (212, 4982) including ID column

# Optional: see combos where phenotype columns were missing
if (length(missing_combos) > 0) {
  print(sprintf("[WARN] %d combos had missing phenotype columns in case_bin_all.", length(missing_combos)))
  print(missing_combos[[1]])
}


colnames(patient_genotype)[1]<-".ID"
colnames(patient_zf)[1]<-".ID"

library(dplyr)

id_col <- ".ID"   # change if your ID column differs

analysis_df <- combo_df %>%
  left_join(patient_zf,       by = id_col) %>%
  left_join(patient_genotype, by = id_col)

table(analysis_df$Sex.x)



analysis_df <- analysis_df %>%
  mutate(
    Sex = ifelse(Sex.x %in% c("M","F"), Sex.x, NA)
  )


# Make sure predictors are correct types
str(analysis_df[, c("has_ZF2","has_ZF1","has_CNV","Sex.x","Age.x")])
analysis_df$Sex.x <- factor(analysis_df$Sex.x, levels = c("F","M"))

library(logistf)
fit <- logistf(
  P36_P37_P34 ~ has_LOF+ Age.x + Sex.x,
  data = analysis_df
)
summary(fit)








library(dplyr)
library(logistf)
library(tidyr)

age_col <- "Age.x"
sex_col <- "Sex"

# Clean sex (drop U if present)
analysis_df <- analysis_df %>%
  mutate(
    Sex_clean = ifelse(.data[[sex_col]] %in% c("M","F"), .data[[sex_col]], NA),
    Sex_clean = factor(Sex_clean, levels = c("F","M"))
  )

fit_one <- function(dat, combo, geno_term, min_events = 1) {
  y <- dat[[combo]]
  if (all(is.na(y)) || length(unique(y[!is.na(y)])) < 2) return(NULL)
  events <- sum(y == 1, na.rm = TRUE)
  if (events < min_events) return(NULL)
  
  # skip if predictor is constant
  x <- dat[[geno_term]]
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) < 2) return(NULL)
  print(as.formula(paste0(combo, " ~ ", geno_term, " + ", age_col, " + Sex_clean")))
  f <- as.formula(paste0(combo, " ~ ", geno_term, " + ", age_col, " + Sex_clean"))
  
  fit <- tryCatch(logistf(f, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  beta <- coef(fit)[geno_term]
  ci   <- suppressWarnings(confint(fit))[geno_term, ]
  pval <- fit$prob[geno_term]
  
  tibble(
    Combination = combo,
    Events = events,
    Predictor = geno_term,
    OR = exp(beta),
    CI_low = exp(ci[1]),
    CI_high = exp(ci[2]),
    p_value = pval
  )
}


# combo columns (use combo_df to avoid grabbing .Label etc.)
id_col <-c( "prevalence_patient_deidentifier_lab_family_patient",".ID")
combo_cols <- setdiff(names(combo_df), id_col)

geno_terms <- c("has_LOF", "has_Missense", "has_Noncoding")  # add "has_ZF2" etc. if you want

res_long <- bind_rows(lapply(combo_cols, function(combo) {
  bind_rows(lapply(geno_terms, function(g) fit_one(analysis_df, combo, g, min_events = 2)))
}))

# FDR correction within each genotype term
res_long <- res_long %>%
  group_by(Predictor) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

res_long %>% arrange(FDR) %>% head(20)


res_wide <- res_long %>%
  mutate(across(c(OR, CI_low, CI_high, p_value, FDR), ~ round(.x, 4))) %>%
  pivot_wider(
    id_cols = c(Combination, Events),
    names_from = Predictor,
    values_from = c(OR, CI_low, CI_high, p_value, FDR),
    names_glue = "{.value}_{Predictor}"
  )

res_wide %>% head()





# 0A) Are combo columns correct?
length(combo_cols)
head(combo_cols, 5)
tail(combo_cols, 5)

# 0B) Do these columns exist in analysis_df?
all(c("has_LOF","has_Missense","has_Noncoding","Age.x","Sex_clean") %in% names(analysis_df))

# 0C) Is Sex_clean actually created?
table(analysis_df$Sex.x ,useNA = "ifany")

# 0D) Try ONE combo, ONE predictor (should return a 1-row tibble)
test1 <- fit_one(analysis_df, combo_cols[1], "has_LOF", min_events = 1)
test1

# 0E) Is res_long empty after running?
nrow(res_long)










library(tibble)
library(logistf)

fit_one <- function(dat, combo, geno_term, min_events = 5, age_col = "Age.x") {
  
  y <- dat[[combo]]
  if (all(is.na(y)) || length(unique(y[!is.na(y)])) < 2) return(NULL)
  
  events <- sum(y == 1, na.rm = TRUE)
  if (events < min_events) return(NULL)
  
  x <- dat[[geno_term]]
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) < 2) return(NULL)
  
  f <- as.formula(paste0("`", combo, "` ~ ", geno_term, " + `", age_col, "` + Sex_clean"))
  fit <- tryCatch(logistf(f, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # coef name is often has_LOFTRUE / has_MissenseTRUE / etc.
  term_in_fit <- grep(paste0("^", geno_term), names(coef(fit)), value = TRUE)
  if (length(term_in_fit) == 0) return(NULL)
  term_in_fit <- term_in_fit[1]
  
  beta <- coef(fit)[term_in_fit]
  ci   <- suppressWarnings(confint(fit))[term_in_fit, ]
  pval <- fit$prob[term_in_fit]
  
  tibble(
    Combination = combo,
    Events = events,
    Predictor = geno_term,
    Term = term_in_fit,
    OR = exp(beta),
    CI_low = exp(ci[1]),
    CI_high = exp(ci[2]),
    p_value = pval
  )
}

fit_one(analysis_df, "P36_P37_P34", "has_LOF", min_events = 1)
fit_one(analysis_df, "P36_P37_P34", "has_Noncoding", min_events = 1)





library(dplyr)

geno_terms <- c("has_LOF", "has_Missense", "has_Noncoding")

res_long <- bind_rows(lapply(combo_cols, function(cc) {
  bind_rows(lapply(geno_terms, function(g) {
    fit_one(analysis_df, cc, g, min_events = 5)
  }))
}))

# FDR within each genotype class
res_long <- res_long %>%
  group_by(Predictor) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  ungroup()


res_long_sig<-res_long[res_long$p_value<0.05,]



#Unadjusted comparison
library(dplyr)
library(logistf)
library(tibble)

fit_one_unadj <- function(dat, combo, geno_term, min_events = 5) {
  
  # outcome checks
  y <- dat[[combo]]
  if (all(is.na(y)) || length(unique(y[!is.na(y)])) < 2) return(NULL)
  
  events <- sum(y == 1, na.rm = TRUE)
  if (events < min_events) return(NULL)
  
  # predictor checks
  x <- dat[[geno_term]]
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) < 2) return(NULL)
  
  # unadjusted model
  f <- as.formula(paste0("`", combo, "` ~ ", geno_term))
  
  fit <- tryCatch(logistf(f, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  # correct coefficient name (e.g. has_LOFTRUE)
  term_in_fit <- grep(paste0("^", geno_term), names(coef(fit)), value = TRUE)
  if (length(term_in_fit) == 0) return(NULL)
  term_in_fit <- term_in_fit[1]
  
  beta <- coef(fit)[term_in_fit]
  ci   <- suppressWarnings(confint(fit))[term_in_fit, ]
  pval <- fit$prob[term_in_fit]
  
  tibble(
    Combination = combo,
    Events = events,
    Predictor = geno_term,
    OR = exp(beta),
    CI_low = exp(ci[1]),
    CI_high = exp(ci[2]),
    p_value = pval
  )
}

geno_terms <- c("has_LOF", "has_Missense", "has_Noncoding")

res_long_unadj <- bind_rows(lapply(combo_cols, function(cc) {
  bind_rows(lapply(geno_terms, function(g) {
    fit_one_unadj(analysis_df, cc, g, min_events = 5)
  }))
}))
res_long_unadj_sig<-res_long_unadj[res_long_unadj$p_value<0.05,]









#Obtain Phenotypes for the top frequent vairants.
#Variant1 c.1017+572C>T
gata2_v1<-gata2_survey_data[gata2_survey_data$prevalence_cdna_change=="c.1017+572C>T",]
gata2_v1_check<-phenotype_gata2[rownames(phenotype_gata2)%in%gata2_v1$prevalence_patient_deidentifier_lab_family_patient,]
pheno_cols <- colnames(gata2_v1_check)
pheno_counts <- gata2_v1_check %>%
  summarise(across(
    all_of(pheno_cols),
    ~ sum(. == "Y", na.rm = TRUE)
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Phenotype",
    values_to = "Count_Y"
  )
pheno_counts_labeled <- pheno_counts %>%
  left_join(
    select_pheno_full %>%
      select(Phenotype, Phenotype_Description_Expanded),
    by = "Phenotype"
  )
total_n <- nrow(gata2_v1_check)

final_pheno_summary_v1 <- pheno_counts_labeled %>%
  mutate(
    Phenotype_Label = ifelse(
      is.na(Phenotype_Description_Expanded),
      Phenotype,
      Phenotype_Description_Expanded
    ),
    Percent = round((Count_Y / total_n) * 100, 1)
  ) %>%
  select(
    Phenotype,
    Phenotype_Label,
    Count_Y,
    Percent
  ) %>%
  arrange(desc(Count_Y))


#v2 R398W,R396W,R396C

gata2_v2<-gata2_survey_data[gata2_survey_data$prevalence_cdna_change%in%c("c.1192C>T","c.1186C>T","c.1187G>A"),]
gata2_v2_check<-phenotype_gata2[rownames(phenotype_gata2)%in%gata2_v2$prevalence_patient_deidentifier_lab_family_patient,]

pheno_cols <- colnames(gata2_v2_check)
pheno_counts <- gata2_v2_check %>%
  summarise(across(
    all_of(pheno_cols),
    ~ sum(. == "Y", na.rm = TRUE)
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Phenotype",
    values_to = "Count_Y"
  )
pheno_counts_labeled <- pheno_counts %>%
  left_join(
    select_pheno_full %>%
      select(Phenotype, Phenotype_Description_Expanded),
    by = "Phenotype"
  )
total_n <- nrow(gata2_v2_check)

final_pheno_summary_v2 <- pheno_counts_labeled %>%
  mutate(
    Phenotype_Label = ifelse(
      is.na(Phenotype_Description_Expanded),
      Phenotype,
      Phenotype_Description_Expanded
    ),
    Percent = round((Count_Y / total_n) * 100, 1)
  ) %>%
  select(
    Phenotype,
    Phenotype_Label,
    Count_Y,
    Percent
  ) %>%
  arrange(desc(Count_Y))

