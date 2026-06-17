# =============================================================================
# Individual Gene — Survival Correlation (Univariate Cox + KM)
# Adapted for: genes as COLUMNS in clinical dataframe
# =============================================================================

library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)

# --- 1. LOAD DATA -------------------------------------------------------------
clinical <- read.csv("clinical.csv", stringsAsFactors = FALSE)

# Clean column names and survival variables
clinical <- clinical %>%
  rename(
    OS_months  = Overall.Survival..Months.,
    PTK2       = PTK2.x
  ) %>%
  select(-PTK2.y) %>%
  mutate(event = as.numeric(status)) %>%
  filter(!is.na(OS_months))

# --- 2. DEFINE GENES ----------------------------------------------------------
sodium_genes <- c("ATP1A1", "SCN8A",  "SLC4A4", "SLC4A7",
                  "SLC5A6", "ATP1A3", "SLC12A3", "SLC9A1",
                  "SLC4A8", "FXYD5",  "SCN1B",  "SLC13A3")

all_genes <- c("PTK2", sodium_genes)

# Confirm all genes are present as columns
all_genes <- all_genes[all_genes %in% colnames(clinical)]
cat("Genes found in dataframe:", paste(all_genes, collapse = ", "), "\n\n")

# =============================================================================
# PART A — UNIVARIATE COX TABLE (continuous expression)
# =============================================================================
cox_results <- lapply(all_genes, function(g) {
  df      <- clinical
  df$expr <- df[[g]]                          # extract gene column by name
  fit     <- coxph(Surv(OS_months, event) ~ expr, data = df)
  s       <- summary(fit)
  data.frame(
    Gene      = g,
    HR        = round(s$conf.int[, "exp(coef)"],  3),
    CI_lo     = round(s$conf.int[, "lower .95"],  3),
    CI_hi     = round(s$conf.int[, "upper .95"],  3),
    P_value   = round(s$coefficients[, "Pr(>|z|)"], 4),
    Sig       = ifelse(s$coefficients[, "Pr(>|z|)"] < 0.05, "*", ""),
    Direction = ifelse(s$conf.int[, "exp(coef)"] > 1,
                       "Risk  (High expr = worse OS)",
                       "Protective (High expr = better OS)")
  )
})

cox_df <- do.call(rbind, cox_results)
rownames(cox_df) <- NULL

# Sort by p-value
cox_df <- cox_df %>% arrange(P_value)

cat("========== Univariate Cox Results ==========\n")
print(cox_df)
write.csv(cox_df, "Gene_Cox_results.csv", row.names = FALSE)

# =============================================================================
# PART B — FOREST PLOT of all genes
# =============================================================================
cox_df$Gene <- factor(cox_df$Gene, levels = rev(cox_df$Gene))

forest_plot <- ggplot(cox_df, aes(x = HR, y = Gene, color = Direction)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), height = 0.25, linewidth = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(values = c(
    "Risk  (High expr = worse OS)"       = "#D73027",
    "Protective (High expr = better OS)" = "#4575B4"
  )) +
  labs(
    title    = "Univariate Cox Regression — Each Gene vs Overall Survival",
    subtitle = "TCGA-SKCM | HR > 1: risk gene | HR < 1: protective gene",
    x        = "Hazard Ratio (95% CI)",
    y        = NULL,
    color    = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  ) +
  # Add p-value labels
  geom_text(aes(x = max(CI_hi) + 0.05,
                label = ifelse(P_value < 0.05,
                               paste0("p=", P_value, " *"),
                               paste0("p=", P_value))),
            hjust = 0, size = 3, color = "black")

print(forest_plot)
ggsave("Forest_plot_all_genes.pdf",
       plot   = forest_plot,
       width  = 10,
       height = 6)

# =============================================================================
# PART C — KAPLAN-MEIER FOR EACH GENE (median split)
# =============================================================================
km_plots <- list()

for (g in all_genes) {

  df           <- clinical
  gene_median  <- median(df[[g]], na.rm = TRUE)
  df$gene_group <- factor(
    ifelse(df[[g]] >= gene_median, "High", "Low"),
    levels = c("High", "Low")
  )

  km_fit <- survfit(Surv(OS_months, event) ~ gene_group, data = df)

  p <- ggsurvplot(
    km_fit,
    data            = df,
    pval            = TRUE,
    pval.method     = TRUE,
    conf.int        = TRUE,
    risk.table      = TRUE,
    risk.table.col  = "strata",
    palette         = c("#D73027", "#4575B4"),
    legend.labs     = c(paste(g, "High"), paste(g, "Low")),
    legend.title    = "Expression",
    xlab            = "Time (Months)",
    ylab            = "Overall Survival Probability",
    title           = paste(g, "Expression — TCGA-SKCM Overall Survival"),
    ggtheme         = theme_classic2(base_size = 12),
    surv.median.line = "hv",
    tables.theme    = theme_cleantable()
  )

  km_plots[[g]] <- p
  ggsave(paste0("KM_", g, ".pdf"),
         plot  = print(p),
         width = 8,
         height = 7)
  cat("KM plot saved for:", g, "\n")
}

# =============================================================================
# PART D — COMBINED KM PANEL (all genes in one PDF)
# =============================================================================
pdf("KM_all_genes_panel.pdf", width = 16, height = 9)
arrange_ggsurvplots(km_plots,
                    ncol  = 4,
                    nrow  = ceiling(length(all_genes) / 4),
                    title = "Individual Gene KM — TCGA-SKCM")
dev.off()

# =============================================================================
# PART E — SUMMARY: WHICH GENES TO INCLUDE IN SIGNATURE
# =============================================================================
cat("\n========== Signature Building Recommendation ==========\n")
cat("\nRisk genes (HR > 1) — include as POSITIVE in signature:\n")
print(cox_df %>% filter(HR > 1) %>% select(Gene, HR, P_value, Sig))

cat("\nProtective genes (HR < 1) — include as NEGATIVE in signature:\n")
print(cox_df %>% filter(HR < 1) %>% select(Gene, HR, P_value, Sig))

cat("\nSignificant genes only (p < 0.05):\n")
print(cox_df %>% filter(P_value < 0.05) %>% select(Gene, HR, CI_lo, CI_hi, P_value, Direction))

cat("\nDone. Check Gene_Cox_results.csv and Forest_plot_all_genes.pdf\n")
