# Knowledge Base


## STL Decomposition

**STL** = Seasonal and Trend decomposition using Loess.

Splits a raw time series into three components:

$$Y_t = T_t + S_t + \epsilon_t$$

| Component | Symbol | What it captures |
|-----------|--------|-----------------|
| Trend | $T_t$ | Slow multi-year direction |
| Seasonal | $S_t$ | Repeating annual cycle |
| Residual | $\epsilon_t$ | Month-by-month anomalies — the signal we use |

---

### Why we use it

Raw NDVI has a strong annual cycle — vegetation greens up every spring and browns every fall. This seasonal pattern is so large it drowns out the slow structural decline happening underneath during the collapse window. STL strips out the seasonal cycle and long-term trend, leaving only $\epsilon_t$: anomalies relative to what is expected for that time of year.

**The residual answers:** *Is this month's NDVI lower than it should be for this time of year?*

VPD and ppt (used as regressors in climate correction) are not deseasonalized — their anomalies are large enough to be detectable in raw values.

---

### The 4 Panels

When you run `plot(decomp)` in R, you see 4 panels stacked vertically. Here is what each one is:

---

#### Panel 1 — Data

The raw `ndvi` column, plotted against time. No transformation applied. This is exactly what you would get from `plot(date, ndvi)` directly from the CSV.

---

#### Panel 2 — Seasonal ($S_t$)

##### `s.window` parameter

Controls how the seasonal component is estimated:

| Setting | Behavior |
|---------|----------|
| `"periodic"` | All years get equal weight. One fixed seasonal shape learned from all years. |
| A number (e.g., `7`) | Only nearby years contribute. Seasonal pattern can slowly evolve over time. |

---

Because `s.window = "periodic"`, the seasonal component is the **mean of each month slot across all years**:

- Typical January = mean of all Januaries
- Typical February = mean of all Februaries
- ...
- Typical December = mean of all Decembers

This gives 12 fixed values. STL then stamps the same pattern across all months in the series — every year looks identical in this panel. The wave never changes shape because `"periodic"` forces one fixed seasonal baseline for the entire 2000–2014 period.

Example — "typical July":

```
Jul 2000: 0.xxx
Jul 2001: 0.xxx
Jul 2002: 0.xxx
...
Jul 2014: 0.xxx

Mean of all N Julies = 0.xxx  ← S_t for every July in the chart
```

---

#### Panel 3 — Trend ($T_t$)

STL first computes the deseasonalized series:

$$Y_t - S_t$$

Then applies Loess to this series using **time** as the x-axis:

- x = month index (1 = Feb 2000, 2 = Mar 2000, ... N = Dec 2014)
- y = deseasonalized NDVI ($Y_t - S_t$)
- Window: `t.window = 19` (R's default for this data — 19 months)
- Distance unit: **months**, not years

For each month $t_0$, Loess looks at the 19 nearest months in time, applies the tricubic weight:

$$W(d_i) = (1 - |d_i|^3)^3 \quad \text{where } d_i = \frac{|t_i - t_0|}{h},\ h = 9$$

Fits a local weighted line, and takes the intercept as $T_{t_0}$. Repeated for all months → a smooth slow-moving curve. See **Weighted Least Squares — Simple Example** below for the full numeric calculation.

---

### Weighted Least Squares — Simple Example

For each month $t_0$, Loess fits a weighted line through the neighboring points. The intercept of that line is the trend value $T_{t_0}$.

Here is the full calculation with 3 simple points centered at $t_0 = 0$:

| $t_i$ | $y_i$ | $w_i$ |
|--------|--------|--------|
| $-1$ | $2$ | $0.5$ |
| $0$ | $4$ | $1.0$ |
| $+1$ | $3$ | $0.5$ |

**Design matrix, weight matrix, and response vector:**

$$X = \begin{bmatrix} 1 & -1 \\ 1 & 0 \\ 1 & 1 \end{bmatrix} \qquad W = \begin{bmatrix} 0.5 & 0 & 0 \\ 0 & 1.0 & 0 \\ 0 & 0 & 0.5 \end{bmatrix} \qquad y = \begin{bmatrix} 2 \\ 4 \\ 3 \end{bmatrix}$$

The first column of $X$ is all ones (for the intercept). The second column is the centered time index $t_i - t_0$.

**Compute $X^TWX$:**

$$X^TWX = \begin{bmatrix} \sum w_i & \sum t_i w_i \\ \sum t_i w_i & \sum t_i^2 w_i \end{bmatrix}$$

$$= \begin{bmatrix} 0.5+1.0+0.5 & (-1)(0.5)+(0)(1.0)+(1)(0.5) \\ (-1)(0.5)+(0)(1.0)+(1)(0.5) & (-1)^2(0.5)+(0)^2(1.0)+(1)^2(0.5) \end{bmatrix} = \begin{bmatrix} 2.0 & 0 \\ 0 & 1.0 \end{bmatrix}$$

**Compute $X^TWy$:**

$$X^TWy = \begin{bmatrix} \sum w_i y_i \\ \sum t_i w_i y_i \end{bmatrix} = \begin{bmatrix} (0.5)(2)+(1.0)(4)+(0.5)(3) \\ (-1)(0.5)(2)+(0)(1.0)(4)+(1)(0.5)(3) \end{bmatrix} = \begin{bmatrix} 6.5 \\ 0.5 \end{bmatrix}$$

**Invert $X^TWX$:**

Because $X^TWX$ is diagonal, the inverse just flips each diagonal value:

$$(X^TWX)^{-1} = \begin{bmatrix} 1/2.0 & 0 \\ 0 & 1/1.0 \end{bmatrix} = \begin{bmatrix} 0.5 & 0 \\ 0 & 1.0 \end{bmatrix}$$

**Solve for $\hat{\beta}$:**

$$\hat{\beta} = (X^TWX)^{-1} X^TWy = \begin{bmatrix} 0.5 & 0 \\ 0 & 1.0 \end{bmatrix} \begin{bmatrix} 6.5 \\ 0.5 \end{bmatrix} = \begin{bmatrix} 3.25 \\ 0.50 \end{bmatrix}$$

**The trend value at $t_0$:**

$$T_{t_0} = \hat{\beta}_0 = 3.25$$

The slope $\hat{\beta}_1 = 0.50$ describes the local direction but is discarded. Only the intercept is kept as the trend value.

Note: the result (3.25) is pulled toward the center point (4) more than a plain average $\frac{2+4+3}{3} = 3.0$ — because the center weight $w = 1.0$ is double the edge weights $w = 0.5$.

---

### Why the matrix collapses to a weighted mean

When the window is symmetric around $t_0$ (e.g., $-1, 0, +1$), the time points balance out: $\sum t_i w_i = 0$. This makes $X^TWX$ diagonal with zeros in the off-diagonal, which decouples $\hat{\beta}_0$ from $\hat{\beta}_1$. The full matrix formula then collapses into two simple expressions.

**The trend (intercept)** is just the weighted mean of the data points:

$$\hat{\beta}_0 = \frac{\sum w_i y_i}{\sum w_i}$$

$$= \frac{(0.5 \times 2) + (1.0 \times 4) + (0.5 \times 3)}{0.5 + 1.0 + 0.5} = \frac{6.5}{2.0} = 3.25$$

**The slope** is a weighted ratio of how the data changes over time:

$$\hat{\beta}_1 = \frac{\sum w_i t_i y_i}{\sum w_i t_i^2}$$

$$= \frac{(-1 \times 0.5 \times 2) + (0 \times 1.0 \times 4) + (1 \times 0.5 \times 3)}{(1^2 \times 0.5) + (0^2 \times 1.0) + (1^2 \times 0.5)} = \frac{0.5}{1.0} = 0.5$$

The slope is computed but discarded — only $\hat{\beta}_0$ is kept as $T_{t_0}$.

The matrix form is what R runs under the hood. The simplified form above is what the math actually reduces to for a symmetric centered window.

---

### Loess: how weights are assigned

Loess fits many small local weighted regressions across overlapping windows rather than one global curve. The weight assigned to each neighboring year follows the **tricubic kernel** (Cleveland, 1979):

$$W(d) = \begin{cases} \left(1 - |d|^3\right)^3 & \text{for } |d| < 1 \\ 0 & \text{for } |d| \geq 1 \end{cases}$$

Where $d$ is the normalized distance from the target year (actual distance ÷ window half-width $h$).

- $|d| = 0$ (same year): weight = 1.0 — full weight
- $|d| = 0.5$ (halfway to edge): weight = 0.67
- $|d| = 1.0$ (at window edge): weight = 0 — dropped completely
- $|d| > 1$ (outside window): weight = 0

---

**Reference:** Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990). STL: A seasonal-trend decomposition procedure based on loess. *Journal of Official Statistics*, 6(1), 3–73.

---

#### Panel 4 — Remainder ($\epsilon_t$)

What is left after removing both seasonal and trend — simple subtraction:

$$\epsilon_t = Y_t - S_t - T_t$$

One value per month, centered near zero. It answers: *Is the value lower than expected for this time of year, after accounting for the long-term trend?* The output has the same shape as the input series, with the seasonal cycle and trend removed.

---

### R implementation

```r
# convert to monthly ts object
x_ts <- ts(x, start = c(2000, 2), frequency = 12)

# run STL
decomp <- stl(x_ts, s.window = "periodic", robust = TRUE)

# plot all 4 panels: raw, seasonal, trend, remainder
plot(decomp)

# extract remainder (Panel 4)
residual <- decomp$time.series[, "remainder"]
```

---

## CSD Rolling Statistics

Input: climate-corrected residual (`clim_resid`), growing-season months only (Mar–Nov), per spatial unit (2000–2014).

**Climate correction** removes shared climate forcing (VPD and ppt) so that `clim_resid` reflects only the biological/structural response:

1. Fit a linear model on the pre-drought baseline (2000–2011) using `lm()` — R's ordinary least squares function, which finds the coefficients $\hat{\beta}$ that minimize the sum of squared errors between observed and predicted values. This estimates $\hat{\beta}_0, \hat{\beta}_1, \hat{\beta}_2$ that best explain NDVI variation from climate during the pre-drought period:

$$\widehat{ndvi\_deseasoned}_t = \hat{\beta}_0 + \hat{\beta}_1 \cdot vpd_t + \hat{\beta}_2 \cdot ppt_t$$

**Illustrative example** — 4 months, 2 predictors (vpd and ppt):

OLS is only meaningful when observations > unknowns. Here we have 4 months and 3 unknowns ($\hat{\beta}_0, \hat{\beta}_1, \hat{\beta}_2$) — no perfect fit exists, so `lm()` finds the coefficients that minimize the sum of squared residuals. In the actual project we have ~108 months (12 years × 9 growing-season months) and only 3 unknowns.

| Month | $vpd$ | $ppt$ | $ndvi\_deseasoned$ ($y$) |
|-------|-------|-------|--------------------------|
| 1 | 0.8 | 120 | 0.45 |
| 2 | 1.2 | 80  | 0.42 |
| 3 | 1.6 | 60  | 0.38 |
| 4 | 1.0 | 100 | 0.44 |

**Step 1 — Build the design matrix X and response vector y:**

$$X = \begin{bmatrix} 1 & 0.8 & 120 \\ 1 & 1.2 & 80 \\ 1 & 1.6 & 60 \\ 1 & 1.0 & 100 \end{bmatrix} \qquad y = \begin{bmatrix} 0.45 \\ 0.42 \\ 0.38 \\ 0.44 \end{bmatrix}$$

The first column of X is all ones (for the intercept $\hat{\beta}_0$).

**Step 2 — Compute $X^TX$ and $X^Ty$:**

$$X^TX = \begin{bmatrix} 4 & 4.6 & 360 \\ 4.6 & 5.64 & 388 \\ 360 & 388 & 34400 \end{bmatrix} \qquad X^Ty = \begin{bmatrix} 1.69 \\ 1.912 \\ 154.4 \end{bmatrix}$$

**Step 3 — Solve the normal equations $(X^TX)\hat{\beta} = X^Ty$:**

$$4\hat{\beta}_0 + 4.6\hat{\beta}_1 + 360\hat{\beta}_2 = 1.69 \quad \text{...(I)}$$
$$4.6\hat{\beta}_0 + 5.64\hat{\beta}_1 + 388\hat{\beta}_2 = 1.912 \quad \text{...(II)}$$
$$360\hat{\beta}_0 + 388\hat{\beta}_1 + 34400\hat{\beta}_2 = 154.4 \quad \text{...(III)}$$

**Eliminate $\hat{\beta}_0$ from (II) using (I):**

Multiply (I) by $\frac{4.6}{4} = 1.15$ and subtract from (II):

$$\underbrace{(4.6 - 1.15 \times 4)}_0\hat{\beta}_0 + \underbrace{(5.64 - 1.15 \times 4.6)}_{0.35}\hat{\beta}_1 + \underbrace{(388 - 1.15 \times 360)}_{-26}\hat{\beta}_2 = \underbrace{1.912 - 1.15 \times 1.69}_{-0.0315}$$

$$0.35\hat{\beta}_1 - 26\hat{\beta}_2 = -0.0315 \quad \text{...(IV)}$$

**Eliminate $\hat{\beta}_0$ from (III) using (I):**

Multiply (I) by $\frac{360}{4} = 90$ and subtract from (III):

$$\underbrace{(360 - 90 \times 4)}_0\hat{\beta}_0 + \underbrace{(388 - 90 \times 4.6)}_{-26}\hat{\beta}_1 + \underbrace{(34400 - 90 \times 360)}_{2000}\hat{\beta}_2 = \underbrace{154.4 - 90 \times 1.69}_{2.3}$$

$$-26\hat{\beta}_1 + 2000\hat{\beta}_2 = 2.3 \quad \text{...(V)}$$

**Eliminate $\hat{\beta}_1$ — solve for $\hat{\beta}_2$:**

From (IV): $\hat{\beta}_1 = \dfrac{-0.0315 + 26\hat{\beta}_2}{0.35}$. Substitute into (V):

$$-26 \times \frac{-0.0315 + 26\hat{\beta}_2}{0.35} + 2000\hat{\beta}_2 = 2.3$$

$$\frac{0.819}{0.35} - \frac{676}{0.35}\hat{\beta}_2 + 2000\hat{\beta}_2 = 2.3$$

$$2.34 - 1931.43\hat{\beta}_2 + 2000\hat{\beta}_2 = 2.3$$

$$68.57\hat{\beta}_2 = -0.04 \implies \hat{\beta}_2 = \frac{-0.04}{68.57} = -0.00058$$

**Back-substitute $\hat{\beta}_2$ into (IV):**

$$0.35\hat{\beta}_1 - 26(-0.00058) = -0.0315$$
$$0.35\hat{\beta}_1 = -0.0315 - 0.0151 = -0.0466$$
$$\hat{\beta}_1 = \frac{-0.0466}{0.35} = -0.133$$

**Back-substitute $\hat{\beta}_1$ and $\hat{\beta}_2$ into (I):**

$$4\hat{\beta}_0 + 4.6(-0.133) + 360(-0.00058) = 1.69$$
$$4\hat{\beta}_0 - 0.612 - 0.209 = 1.69$$
$$4\hat{\beta}_0 = 2.511 \implies \hat{\beta}_0 = 0.628$$

$$\hat{\beta}_0 = 0.628 \qquad \hat{\beta}_1 = -0.133 \qquad \hat{\beta}_2 = -0.00058$$

**Fitted model:** $\widehat{ndvi} = 0.628 - 0.133 \cdot vpd - 0.00058 \cdot ppt$

**Predictions and residuals:**

| Month | $\hat{y}$ | $y$ | residual |
|-------|-----------|-----|----------|
| 1 | $0.628 - 0.133(0.8) - 0.00058(120) = 0.452$ | 0.45 | −0.002 |
| 2 | $0.628 - 0.133(1.2) - 0.00058(80)  = 0.422$ | 0.42 | −0.002 |
| 3 | $0.628 - 0.133(1.6) - 0.00058(60)  = 0.380$ | 0.38 | −0.000 |
| 4 | $0.628 - 0.133(1.0) - 0.00058(100) = 0.437$ | 0.44 | +0.003 |

No month fits perfectly — `lm()` chose these coefficients because they minimize $(-0.002)^2 + (-0.002)^2 + (0.000)^2 + (0.003)^2 = 0.000017$. Any other coefficients would give a larger total.

2. Using those same coefficients from step 1, plug in each month's observed VPD and ppt for the full 2000–2014 series to get the climate-predicted NDVI for each month.

3. Subtract to isolate the residual:

$$clim\_resid_t = ndvi\_deseasoned_t - \hat{\beta}_0 - \hat{\beta}_1 \cdot vpd_t - \hat{\beta}_2 \cdot ppt_t$$

Training on 2000–2011 ensures the coefficients capture the pre-drought climate–NDVI relationship and are not contaminated by the collapse signal.

A rolling window slides across the series one month at a time. At each position, three statistics are computed on the values inside the window.

---

### 1. Lag-1 Autocorrelation ($r_1$)

**Question:** does this month's residual predict next month's?

$$r_1 = \frac{\sum(\epsilon_t - \mu)(\epsilon_{t+1} - \mu)}{\sum(\epsilon_t - \mu)^2}$$

**Illustrative example** — hypothetical 5-month window:

| Month | $\epsilon_t$ | $\epsilon_t - \mu$ |
|-------|-------------|-------------------|
| 1   | −0.030 | −0.017 |
| 2   | −0.050 | −0.037 |
| 3   | +0.005 | +0.018 |
| 4   | +0.015 | +0.028 |
| 5   | −0.005 | +0.008 |

$$\mu = \frac{-0.030 + (-0.050) + 0.005 + 0.015 + (-0.005)}{5} = \frac{-0.065}{5} = -0.013$$

Check: deviations must sum to zero: $-0.017 - 0.037 + 0.018 + 0.028 + 0.008 = 0$ ✓

**Numerator** — multiply each consecutive pair of deviations:

$$(-0.017)(-0.037) + (-0.037)(+0.018) + (+0.018)(+0.028) + (+0.028)(+0.008)$$
$$= 0.000629 - 0.000666 + 0.000504 + 0.000224 = +0.000691$$

**Denominator** — sum of all squared deviations:

$$(-0.017)^2 + (-0.037)^2 + (0.018)^2 + (0.028)^2 + (0.008)^2$$
$$= 0.000289 + 0.001369 + 0.000324 + 0.000784 + 0.000064 = 0.002830$$

$$r_1 = \frac{0.000691}{0.002830} \approx \mathbf{0.24}$$

This value is assigned to **Month 5** — the last month of the window.

When the system is healthy, a below-normal month is not predictive of the next — $r_1$ stays low. As the collapse deepens, $r_1$ rises closer to 1.

| $r_1$ value | Meaning |
|-------------|---------|
| Close to 0 | Each month is independent — system recovers quickly from shocks |
| Close to 1 | Low month stays low next month — system recovering slowly |

**Rising $r_1$** = system losing resilience = critical slowing down signal.

**In the real computation:** window size = 54 growing-season months (~40% of N). Same formula, 54 values instead of 5. Each $r_1$ value is assigned to the **last month of its window** — the first valid value is **Nov 2005** (end of the first complete 54-month window starting Mar 2000). The window slides one month at a time → 82 valid $r_1$ values output. "Consecutive" here means consecutive in the growing-season series: the month after Nov is Mar of the next year (Dec–Feb are excluded).

---

### 2. Rolling Variance ($\sigma^2$)

**Question:** how much do residuals fluctuate inside the window?

$$\sigma^2 = \frac{1}{w-1}\sum_{t=1}^{w}(\epsilon_t - \mu)^2$$

**Illustrative example** — same hypothetical 5-month window ($\mu = -0.013$):

| Month | $\epsilon_t$ | $\epsilon_t - \mu$ | $(\epsilon_t - \mu)^2$ |
|-------|-------------|-------------------|----------------------|
| 1   | −0.030 | −0.017 | 0.000289 |
| 2   | −0.050 | −0.037 | 0.001369 |
| 3   | +0.005 | +0.018 | 0.000324 |
| 4   | +0.015 | +0.028 | 0.000784 |
| 5   | −0.005 | +0.008 | 0.000064 |

$$\sigma^2 = \frac{0.000289 + 0.001369 + 0.000324 + 0.000784 + 0.000064}{5-1} = \frac{0.002830}{4} = \mathbf{0.000708}$$

$\sigma^2 = 0.000708$ is in squared NDVI units, which is hard to interpret directly. Taking the square root gives the standard deviation $\sigma = \sqrt{0.000708} \approx 0.027$ NDVI units — a typical month in this window sits about 0.027 NDVI units above or below the window mean. As the collapse deepens and residuals swing further from the mean, $\sigma^2$ rises.

**Why $w-1$ not $w$:** dividing by $w-1$ (Bessel's correction) gives an unbiased estimate of true variance from a sample.

$\sigma^2$ is assigned to the **last month of the window**. Same rule as $r_1$.

**Rising $\sigma^2$** = residuals swinging wider = system less stable, being pushed further from equilibrium.

---

### 3. Skewness ($\gamma$)

**Question:** is the distribution of residuals symmetric, or lopsided to one side?

$$\gamma = -\frac{\frac{1}{w}\sum(\epsilon_t - \mu)^3}{\left(\frac{1}{w}\sum(\epsilon_t - \mu)^2\right)^{3/2}}$$

**Illustrative example** — same hypothetical 5-month window ($\mu = -0.013$):

| Month | $\epsilon_t - \mu$ | $(\epsilon_t - \mu)^3$ |
|-------|-------------------|----------------------|
| 1   | −0.017 | −0.0000049 |
| 2   | −0.037 | −0.0000507 |
| 3   | +0.018 | +0.0000058 |
| 4   | +0.028 | +0.0000220 |
| 5   | +0.008 | +0.0000005 |

**Numerator:**
$$\frac{1}{5} \times (-0.0000049 - 0.0000507 + 0.0000058 + 0.0000220 + 0.0000005) = \frac{-0.0000273}{5} = -0.00000546$$

**Denominator** — $\left(\frac{1}{w}\sum(\epsilon_t - \mu)^2\right)^{3/2}$, using the sum of squares from above ($0.002830$):
$$\left(\frac{0.002830}{5}\right)^{3/2} = (0.000566)^{3/2} = 0.000566 \times \sqrt{0.000566} = 0.000566 \times 0.02379 \approx 0.0000135$$

$$\gamma = -\frac{-0.00000546}{0.0000135} \approx \mathbf{+0.40}$$

The raw skewness is negative (the −0.037 deviation in month 2 dominates the cubed sum, pulling it left). The sign is flipped by convention, so that rising $\gamma$ corresponds to residuals becoming more left-skewed. (Note: this is a project-specific convention, not something the `earlywarnings` package itself does — it uses `abs(skewness(...))` rather than a sign flip.) Assigned to the **last month of the window**, same rule as $r_1$ and $\sigma^2$.

| $\gamma$ value | Meaning |
|----------------|---------|
| $\approx 0$ | Symmetric — equally likely to be above or below normal |
| Shifting positive or negative | System spending more time on one side — flickering near tipping point |

---

### What the output looks like

Each statistic produces a time series of valid values once the first full window is complete. CSD theory predicts a signature pattern as a system approaches a tipping point:

- $r_1$ trends toward 1: the system's memory is increasing — it's losing resilience and struggling to bounce back to its normal state after a perturbation.
- $\sigma^2$ trends upward: the system fluctuates more wildly as the equilibrium it sits in becomes a shallower, flatter "valley" with less restoring force.
- $\gamma$ shifts significantly (direction not universal): as the system nears the tipping point, its distribution becomes asymmetrical because it starts leaning toward the alternative state it's about to flip into. Unlike $r_1$ and $\sigma^2$, whether $\gamma$ rises or falls depends on which direction that alternative state lies in.
- Stable systems: none of the three show a consistent trend

---

## Kendall's Tau ($\tau$)

**Question:** does a CSD statistic (r₁, variance, skewness) rise monotonically over time?

Kendall's τ looks at all possible **pairs** of observations. For n time steps, that is $\frac{n(n-1)}{2}$ pairs.

For each pair $(i, j)$ where $date_i < date_j$:
- **Concordant** if $value_j > value_i$ — later value is higher, consistent with upward trend
- **Discordant** if $value_j < value_i$ — later value is lower, inconsistent with upward trend

$$\tau = \frac{C - D}{C + D}$$

where $C$ = number of concordant pairs, $D$ = number of discordant pairs.

| $\tau$ | Meaning |
|--------|---------|
| $+1$ | Every later observation is higher — perfect upward trend |
| $-1$ | Every later observation is lower — perfect downward trend |
| $0$ | No monotonic trend |

**Illustrative example** — 4 time steps of r₁:

| Time step | r₁ |
|-----------|-----|
| 1 | 0.10 |
| 2 | 0.18 |
| 3 | 0.15 |
| 4 | 0.25 |

All pairs (6 total):

| Pair | Later > Earlier? | Concordant? |
|------|-----------------|-------------|
| (1,2): 0.18 vs 0.10 | yes | C |
| (1,3): 0.15 vs 0.10 | yes | C |
| (1,4): 0.25 vs 0.10 | yes | C |
| (2,3): 0.15 vs 0.18 | no  | D |
| (2,4): 0.25 vs 0.18 | yes | C |
| (3,4): 0.25 vs 0.15 | yes | C |

$C = 5,\ D = 1$

$$\tau = \frac{5 - 1}{5 + 1} = \frac{4}{6} \approx +0.67$$

Strong upward trend — most later observations are higher.

**Key advantage over Pearson correlation:** does not assume the trend is linear, only that it is monotonic — later observations tend to be higher. Appropriate here since CSD statistics don't need to rise at a constant rate, just consistently over time.

**In R:**
```r
cor.test(as.numeric(date), r1, method = "kendall")
```
`as.numeric(date)` converts dates to days since 1970-01-01. The p-value tests whether $\tau$ is significantly different from 0.

