library(tidyverse)
library(car)

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

cat("Loading pixel_tau.csv...\n")
d <- read_csv("pixel_approach/data/processed/pixel_tau.csv",
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
