# DEPRIORITIZED (moved from pixel_approach/src/07_regression.R): this whole
# script regresses the CSD predictors against Actual_Drop_15_16, which is
# no longer the trusted ground truth -- real-world validation (MTBS fire,
# LCMAP mortality, USFS ADS beetle -- 09/10/08) is, and Actual_Drop doesn't
# agree with it on cluster severity ordering. Not auto-invoked by any
# master pipeline script, so nothing breaks by leaving it here. Run paths
# below are relative to the repo root (as invoked via
# `Rscript pixel_approach/src/deprecated/07_regression.R`), not to this
# file's own location, so moving it didn't require any path changes.

library(tidyverse)
library(car)

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

cat("Loading ndvi_pixel_tau.csv...\n")
d <- read_csv("pixel_approach/data/processed/ndvi_pixel_tau.csv",
              show_col_types = FALSE) |>
  filter(!is.na(Actual_Drop_15_16), !is.na(tau_r1),
         !is.na(tau_var), !is.na(tau_skew))

cat(sprintf("✓ %d pixels with complete data\n\n", nrow(d)))

# --------------------------------------------------------------------------
# 1. Two-predictor model: tau_r1 + tau_var
# --------------------------------------------------------------------------
m2 <- lm(Actual_Drop_15_16 ~ tau_r1 + tau_var, data = d)
cat("--- Model 1: tau_r1 + tau_var ---\n")
print(summary(m2)$coefficients)
cat(sprintf("R-squared: %.5f\n\n", summary(m2)$r.squared))

# --------------------------------------------------------------------------
# 2. Three-predictor model: tau_r1 + tau_var + tau_skew
# --------------------------------------------------------------------------
m3 <- lm(Actual_Drop_15_16 ~ tau_r1 + tau_var + tau_skew, data = d)
cat("--- Model 2: tau_r1 + tau_var + tau_skew ---\n")
print(summary(m3)$coefficients)
cat(sprintf("R-squared: %.5f\n\n", summary(m3)$r.squared))

# --------------------------------------------------------------------------
# 3. VIF — check multicollinearity before accepting 3-predictor model
# --------------------------------------------------------------------------
cat("--- VIF (3-predictor model) ---\n")
print(vif(m3))
cat("(VIF near 1 = no multicollinearity; > 5 = concern)\n\n")

# --------------------------------------------------------------------------
# 4. Model comparison
# --------------------------------------------------------------------------
cat("--- ANOVA (model comparison) ---\n")
print(anova(m2, m3))

# --------------------------------------------------------------------------
# 4. Coefficient plot
# --------------------------------------------------------------------------
coefs <- broom::tidy(m3, conf.int = TRUE) |>
  filter(term != "(Intercept)")

p <- ggplot(coefs, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, linewidth = 0.6) +
  geom_point(size = 3) +
  labs(
    title    = "Regression Coefficients: CSD Predictors of NDVI Collapse",
    subtitle = sprintf("Outcome: Actual_Drop_15_16 | R² = %.3f | n = %s",
                       summary(m3)$r.squared, format(nrow(d), big.mark = ",")),
    x = "Coefficient estimate (95% CI)",
    y = NULL
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

ggsave("pixel_approach/charts/regression_coefs.pdf",
       plot = p, width = 7, height = 4, dpi = 150)
cat("✓ Saved: pixel_approach/charts/regression_coefs.pdf\n")
