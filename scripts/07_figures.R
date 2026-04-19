# =============================================================================
# Script 07: Publication Figures
# =============================================================================

source("scripts/utils.R")
library(ggplot2)

# ============================================================
# Figure 1: Combined forest plot (all three targets)
# ============================================================

mr <- fread(file.path(dir_tables, "mr_ivw_corrected.csv"))
mr[, OR := exp(b)]
mr[, OR_lo := exp(b - 1.96 * se)]
mr[, OR_hi := exp(b + 1.96 * se)]

# Clean outcome labels
mr[, outcome_short := outcome]
mr[grepl("gram-neg", outcome_short), outcome_short := "Other septicaemia"]

# Order: bacterial top, viral middle, control bottom
mr[, category := factor(category,
                        levels = c("negative_control", "viral", "bacterial"))]
mr[, outcome_short := reorder(outcome_short, as.numeric(category))]

p1 <- ggplot(mr, aes(x = OR, y = outcome_short, colour = category)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0.25, linewidth = 0.4) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  facet_wrap(~target, ncol = 3) +
  scale_colour_manual(
    values = c("bacterial" = "#C0392B", "viral" = "#2980B9",
               "negative_control" = "#7F8C8D"),
    labels = c("Bacterial", "Viral", "Negative control"),
    name = NULL
  ) +
  scale_x_continuous(trans = "log2", breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
  labs(x = "Odds ratio per 1 SD decrease in LDL-C (95% CI, log scale)",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8))

ggsave(file.path(dir_figures, "figure1_forest_combined.pdf"),
       p1, width = 10, height = 4.5, dpi = 300)
ggsave(file.path(dir_figures, "figure1_forest_combined.png"),
       p1, width = 10, height = 4.5, dpi = 300)
log_msg("Figure 1 saved")

# ============================================================
# Figure 2: Required sample size roadmap
# ============================================================

pow <- fread(file.path(dir_tables, "power_calculations_enhanced.csv"))

outcome_labels <- c(
  AB1_OTHER_SEPSIS    = "Other septicaemia",
  AB1_STREPTO_SEPSIS  = "Streptococcal septicaemia",
  J10_PNEUMOPNEUMO    = "Pneumococcal pneumonia",
  J10_INFLUENZA       = "Influenza",
  AB1_HERPES_SIMPLEX  = "Herpes simplex",
  AB1_ZOSTER          = "Herpes zoster",
  ST19_FRACT_FOREA    = "Forearm fracture"
)
pow[, outcome_short := outcome_labels[outcome]]

pow_long <- melt(pow, id.vars = c("target", "outcome_short"),
                 measure.vars = c("n_total", "required_N_bonf"),
                 variable.name = "type", value.name = "N")
pow_long[, type := ifelse(type == "n_total",
                          "Current (FinnGen R12)",
                          "Required (OR=0.85, 80% power)")]

# Order outcomes by required N for PCSK9
outcome_order <- pow[target == "PCSK9"][order(required_N_bonf), outcome_short]
pow_long[, outcome_short := factor(outcome_short, levels = outcome_order)]

p2 <- ggplot(pow_long, aes(x = N, y = outcome_short, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~target, ncol = 3) +
  scale_x_log10(
    labels = function(x) ifelse(x >= 1e6, paste0(round(x / 1e6, 1), "M"),
                          ifelse(x >= 1e3, paste0(round(x / 1e3), "K"), x)),
    breaks = c(1e5, 1e6, 1e7, 1e8, 1e9)
  ) +
  scale_fill_manual(
    values = c("Current (FinnGen R12)" = "#3498DB",
               "Required (OR=0.85, 80% power)" = "#E74C3C"),
    name = NULL
  ) +
  labs(x = "Outcome GWAS sample size (log scale)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8))

ggsave(file.path(dir_figures, "figure2_power_gap.pdf"),
       p2, width = 10, height = 5, dpi = 300)
ggsave(file.path(dir_figures, "figure2_power_gap.png"),
       p2, width = 10, height = 5, dpi = 300)
log_msg("Figure 2 saved")

# ============================================================
# Figure 3: Power simulation by biobank combination
# ============================================================

simulate_power <- function(r2, case_prop, n_seq, alpha) {
  sapply(n_seq, function(n) {
    var_beta <- 1 / (n * r2 * case_prop * (1 - case_prop))
    se_val <- sqrt(var_beta)
    z_alpha <- qnorm(1 - alpha / 2)
    pnorm(abs(log(0.85)) / se_val - z_alpha)
  })
}

n_seq <- 10^seq(5, 8, by = 0.05)
bonf_alpha <- 0.05 / 21
r2_pcsk9 <- 0.0017

outcomes_sim <- data.frame(
  outcome = c("Other septicaemia", "Influenza", "Herpes zoster",
              "Streptococcal septicaemia", "Pneumococcal pneumonia"),
  case_prop = c(17133 / 456181, 11558 / 427096, 7132 / 487448,
                3239 / 483555, 1625 / 417163)
)

sim_data <- list()
for (i in seq_len(nrow(outcomes_sim))) {
  pow_vals <- simulate_power(r2_pcsk9, outcomes_sim$case_prop[i], n_seq, bonf_alpha)
  sim_data[[i]] <- data.frame(N = n_seq, power = pow_vals,
                               outcome = outcomes_sim$outcome[i])
}
sim_df <- do.call(rbind, sim_data)

biobanks <- data.frame(
  label = c("FinnGen R12\n(0.5M)", "FinnGen+UKBB\n(1M)",
            "FinnGen+UKBB+MVP\n(2M)", "GBMI target\n(5M)",
            "Future\n(10M)"),
  N = c(5e5, 1e6, 2e6, 5e6, 1e7)
)

p3 <- ggplot(sim_df, aes(x = N, y = power * 100, colour = outcome)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey40",
             linewidth = 0.3) +
  geom_vline(data = biobanks, aes(xintercept = N), linetype = "dotted",
             colour = "grey60", linewidth = 0.3) +
  geom_text(data = biobanks, aes(x = N, y = 95, label = label),
            inherit.aes = FALSE, size = 2.3, colour = "grey40", hjust = 0.5) +
  scale_x_log10(
    labels = function(x) ifelse(x >= 1e6, paste0(x / 1e6, "M"),
                                paste0(x / 1e3, "K")),
    breaks = c(1e5, 5e5, 1e6, 2e6, 5e6, 1e7, 5e7, 1e8)
  ) +
  scale_colour_brewer(palette = "Set1", name = NULL) +
  labs(title = "PCSK9-proxied LDL-lowering: power by outcome GWAS sample size",
       subtitle = "OR = 0.85 per 1 SD LDL-C decrease, Bonferroni-corrected alpha",
       x = "Outcome GWAS total sample size (log scale)",
       y = "Statistical power (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(dir_figures, "figure3_power_simulation.pdf"),
       p3, width = 9, height = 5, dpi = 300)
ggsave(file.path(dir_figures, "figure3_power_simulation.png"),
       p3, width = 9, height = 5, dpi = 300)
log_msg("Figure 3 saved")

log_msg("Script 07 complete.")
