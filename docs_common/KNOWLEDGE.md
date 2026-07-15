# Knowledge Base

## STL Decomposition

**STL** = Seasonal and Trend decomposition using Loess.

Splits a raw time series into three components:

$$Y_t = T_t + S_t + \epsilon_t$$

| Component | Symbol | What it captures |
|-----------|--------|-----------------|
| Trend | $T_t$ | Slow multi-year direction |
| Seasonal | $S_t$ | Repeating annual cycle |
| Residual | $\epsilon_t$ | Month-by-month anomalies |

---

### Why we use it

A raw monthly time series often has a strong annual cycle (e.g., a vegetation index that greens up every spring and browns every fall). This seasonal pattern is so large it can drown out slower, multi-year structural change happening underneath. STL strips out the seasonal cycle and long-term trend, leaving only $\epsilon_t$: anomalies relative to what is expected for that time of year.

**The residual answers:** *Is this month's value lower than it should be for this time of year?*

---

### The 4 Panels

When you run `plot(decomp)` in R, you see 4 panels stacked vertically. Here is what each one is:

---

#### Panel 1 — Data

The raw value column, plotted against time. No transformation applied. This is exactly what you would get from `plot(date, x)` directly from the data.

---

#### Panel 2 — Seasonal ($S_t$)

##### `s.window` parameter

Controls how the seasonal component is estimated:

| Setting | Behavior |
|---------|----------|
| `"periodic"` | All years get equal weight. One fixed seasonal shape learned from all years. |
| A number (e.g., `3`) | Only nearby years contribute. Seasonal pattern can slowly evolve over time. |

---

Because `s.window = "periodic"`, the seasonal component is the **mean of each month slot across all years**:

- Typical January = mean of all Januaries
- Typical February = mean of all Februaries
- ...
- Typical December = mean of all Decembers

This gives 12 fixed values. STL then stamps the same pattern across all months in the series — every year looks identical in this panel. The wave never changes shape because `"periodic"` forces one fixed seasonal baseline for the entire series.

Illustrative example — "typical July" over a 15-year record (values here are made up for the math, not measured data):

| Year | Y1 | Y2 | Y3 | Y4 | Y5 | Y6 | Y7 | Y8 | Y9 | Y10 | Y11 | Y12 | Y13 | Y14 | Y15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| July value | 0.60 | 0.62 | 0.61 | 0.63 | 0.65 | 0.64 | 0.63 | 0.62 | 0.60 | 0.58 | 0.59 | 0.57 | 0.55 | 0.53 | 0.50 |

$$\text{Mean of all 15 Julies} = 0.595 \quad \leftarrow S_t \text{ for } \textbf{every} \text{ July in the chart, Y1 through Y15 alike}$$

---

##### The math behind a numeric `s.window`

With a number instead of `"periodic"`, the seasonal estimate for each year is built from a **local** loess smooth of that year's subseries (e.g., just the Julys) instead of a plain mean over all years — so nearby years count more than distant ones, using the same tricubic weighting idea as the trend panel below.

That local smooth is not the final answer, though: STL then runs it through an additional **low-pass filter** (two moving averages the length of one full year, then a length-3 moving average, then one more loess pass) and subtracts that filtered result back off, before the seasonal component is finalized. This extra step exists to strip out any long-term trend that leaked into the seasonal estimate — but it also means the *net* influence of a neighboring year on the final $S_t$ is not simply the tricubic weight from the first step alone.

So: a numeric `s.window` does let the seasonal baseline drift across years (unlike `"periodic"`'s single fixed value), but because of the extra low-pass filtering step above, the exact amount of drift contributed by any one neighboring year is not simply its tricubic weight from the first step — there isn't a simple hand-computable formula for the net effect the way there is for the trend panel below.

(`s.window` should be odd: an odd number splits evenly into "N years before, itself, N years after" — e.g. 3 = 1 before + self + 1 after. An even number can't center on the target year this way. R's `?stl` docs recommend `s.window` be odd and at least 7, per Cleveland et al.'s original STL paper — not a hard rule R enforces, just what they found reliable. Paper: [Cleveland, R.B., Cleveland, W.S., McRae, J.E., Terpenning, I. (1990). STL: A Seasonal-Trend Decomposition Procedure Based on Loess. *Journal of Official Statistics*, 6(1), 3-73.](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/stl-a-seasonal-trend-decomposition-procedure-based-on-loess.pdf))

---

#### Panel 3 — Trend ($T_t$)

STL first computes the deseasonalized series:

$$Y_t - S_t$$

Then applies Loess to this series using **time** as the x-axis:

- x = month index (1 = first month of series, 2 = second month, ... N = last month)
- y = deseasonalized value ($Y_t - S_t$)
- Window: `t.window` — must be odd; a larger window means a smoother, slower-moving trend. Illustrated below with `t.window = 19` purely as an example.
- Distance unit: **months**, not years

For each month $t_0$, Loess looks at the `t.window` nearest months in time, applies the tricubic weight:

$$W(d_i) = (1 - |d_i|^3)^3 \quad \text{where } d_i = \frac{|t_i - t_0|}{h},\ h = \frac{t.window - 1}{2}$$

(with the illustrative `t.window = 19` above, $h = 9$)

Fits a local weighted line, and takes the intercept as $T_{t_0}$. Repeated for all months → a smooth slow-moving curve. See **Weighted Least Squares — Simple Example** below for the full numeric calculation.

---

### Weighted Least Squares — Simple Example

For each month $t_0$, Loess fits a weighted line through the neighboring points. The intercept of that line is the trend value $T_{t_0}$.

Here is the full calculation with 3 simple points centered at $t_0 = 0$:

| $t_i$ | $y_i$ | $w_i$ |
|--------|--------|--------|
| $-1$ | $2$ | $0.67$ |
| $0$ | $4$ | $1.0$ |
| $+1$ | $3$ | $0.67$ |

(These are the tricubic weights for a neighbor at half the window's half-width: $W(0.5) = (1-0.5^3)^3 \approx 0.67$ — see the general tricubic formula below.)

**Design matrix, weight matrix, and response vector:**

$$X = \begin{bmatrix} 1 & -1 \\ 1 & 0 \\ 1 & 1 \end{bmatrix} \qquad W = \begin{bmatrix} 0.67 & 0 & 0 \\ 0 & 1.0 & 0 \\ 0 & 0 & 0.67 \end{bmatrix} \qquad y = \begin{bmatrix} 2 \\ 4 \\ 3 \end{bmatrix}$$

The first column of $X$ is all ones (for the intercept). The second column is the centered time index $t_i - t_0$.

**Compute $X^TWX$:**

$$X^TWX = \begin{bmatrix} \sum w_i & \sum t_i w_i \\ \sum t_i w_i & \sum t_i^2 w_i \end{bmatrix}$$

$$= \begin{bmatrix} 0.67+1.0+0.67 & (-1)(0.67)+(0)(1.0)+(1)(0.67) \\ (-1)(0.67)+(0)(1.0)+(1)(0.67) & (-1)^2(0.67)+(0)^2(1.0)+(1)^2(0.67) \end{bmatrix} = \begin{bmatrix} 2.34 & 0 \\ 0 & 1.34 \end{bmatrix}$$

**Compute $X^TWy$:**

$$X^TWy = \begin{bmatrix} \sum w_i y_i \\ \sum t_i w_i y_i \end{bmatrix} = \begin{bmatrix} (0.67)(2)+(1.0)(4)+(0.67)(3) \\ (-1)(0.67)(2)+(0)(1.0)(4)+(1)(0.67)(3) \end{bmatrix} = \begin{bmatrix} 7.35 \\ 0.67 \end{bmatrix}$$

**Invert $X^TWX$:**

Because $X^TWX$ is diagonal, the inverse just flips each diagonal value:

$$(X^TWX)^{-1} = \begin{bmatrix} 1/2.34 & 0 \\ 0 & 1/1.34 \end{bmatrix} = \begin{bmatrix} 0.4274 & 0 \\ 0 & 0.7463 \end{bmatrix}$$

**Solve for $\hat{\beta}$:**

$$\hat{\beta} = (X^TWX)^{-1} X^TWy = \begin{bmatrix} 0.4274 & 0 \\ 0 & 0.7463 \end{bmatrix} \begin{bmatrix} 7.35 \\ 0.67 \end{bmatrix} = \begin{bmatrix} 3.14 \\ 0.50 \end{bmatrix}$$

**The trend value at $t_0$:**

$$T_{t_0} = \hat{\beta}_0 = 3.14$$

The slope $\hat{\beta}_1 = 0.50$ describes the local direction but is discarded. Only the intercept is kept as the trend value.

Note: the result (3.14) is pulled toward the center point (4) more than a plain average $\frac{2+4+3}{3} = 3.0$ — because the center weight $w = 1.0$ is higher than the edge weights $w = 0.67$.

---

### Why the matrix collapses to a weighted mean

When the window is symmetric around $t_0$ (e.g., $-1, 0, +1$), the time points balance out: $\sum t_i w_i = 0$. This makes $X^TWX$ diagonal with zeros in the off-diagonal, which decouples $\hat{\beta}_0$ from $\hat{\beta}_1$. The full matrix formula then collapses into two simple expressions.

**The trend (intercept)** is just the weighted mean of the data points:

$$\hat{\beta}_0 = \frac{\sum w_i y_i}{\sum w_i}$$

$$= \frac{(0.67 \times 2) + (1.0 \times 4) + (0.67 \times 3)}{0.67 + 1.0 + 0.67} = \frac{7.35}{2.34} = 3.14$$

**The slope** is a weighted ratio of how the data changes over time:

$$\hat{\beta}_1 = \frac{\sum w_i t_i y_i}{\sum w_i t_i^2}$$

$$= \frac{(-1 \times 0.67 \times 2) + (0 \times 1.0 \times 4) + (1 \times 0.67 \times 3)}{(1^2 \times 0.67) + (0^2 \times 1.0) + (1^2 \times 0.67)} = \frac{0.67}{1.34} = 0.5$$

The slope is computed but discarded — only $\hat{\beta}_0$ is kept as $T_{t_0}$.

The matrix form is the general definition of weighted least squares. The simplified form above is what that math reduces to for a symmetric centered window — but the *actual* numerical method a given implementation uses to solve it (e.g. QR decomposition rather than explicit matrix inversion, for numerical stability) may differ from either form shown here while still producing the same fitted values.

---

### Loess: how weights are assigned

Loess fits many small local weighted regressions across overlapping windows rather than one global curve. The weight assigned to each neighboring year follows the **tricubic kernel** (Cleveland, 1979):

$$W(d) = \begin{cases} \left(1 - |d|^3\right)^3 & \text{for } |d| < 1 \\ 0 & \text{for } |d| \geq 1 \end{cases}$$

Where $d$ is the normalized distance from the target year (actual distance ÷ window half-width $h$).

- $|d| = 0$ (same year): weight = 1.0 — full weight
- $|d| = 0.5$ (halfway to edge): weight = 0.67
- $|d| = 1.0$ (at window edge): weight = 0 — dropped completely
- $|d| > 1$ (outside window): weight = 0

For the trend panel, this tricubic weight directly determines each neighboring month's influence on $T_{t_0}$ (illustrated above with `t.window = 19`, i.e. $h=9$: full weight at the target month, decaying to zero by 9 months out, since $d=9/9=1.0$ hits the formula's own cutoff). For the seasonal panel with a numeric `s.window`, this same tricubic weight is only the *first* of several steps — see the caveat under "The math behind a numeric `s.window`" above.

---

**Reference:** [Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990). STL: A seasonal-trend decomposition procedure based on loess. *Journal of Official Statistics*, 6(1), 3–73.](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/stl-a-seasonal-trend-decomposition-procedure-based-on-loess.pdf)

---

#### Panel 4 — Remainder ($\epsilon_t$)

What is left after removing both seasonal and trend — simple subtraction:

$$\epsilon_t = Y_t - S_t - T_t$$

One value per month, centered near zero. It answers: *Is the value lower than expected for this time of year, after accounting for the long-term trend?* The output has the same shape as the input series, with the seasonal cycle and trend removed.

---

### R implementation

```r
# convert to monthly ts object (start = c(year, month) of the first observation)
x_ts <- ts(x, start = c(1, 1), frequency = 12)

# run STL
decomp <- stl(x_ts, s.window = "periodic", robust = TRUE)

# plot all 4 panels: raw, seasonal, trend, remainder
plot(decomp)

# extract remainder (Panel 4)
residual <- decomp$time.series[, "remainder"]
```

Per R's own documentation (`?stl`), `stl()` is a direct port of Cleveland et al.'s original Fortran STL code — it does **not** call R's general-purpose `loess()` function internally. So the tricubic-kernel description above matches the actual published algorithm, not just a generic approximation of "some kind of local regression." That said, the exact *numerical* method used to solve each local weighted fit (e.g. for stability) isn't something this doc verifies — only that the textbook formulas above match the algorithm's published definition and, where checked against real R output, its results.

---

## CSD (Critical Slowing Down) Rolling Statistics

**CSD** = Critical Slowing Down — a theory from dynamical systems that a system losing resilience (approaching a tipping point, e.g. an ecosystem collapse or regime shift) tends to show detectable warning signs beforehand: it recovers more slowly from small perturbations, fluctuates more, and its fluctuations become asymmetric. The three statistics below (lag-1 autocorrelation, variance, skewness), tracked over a rolling window, are the standard way to measure those warning signs.

**Input:** a climate-corrected residual series ($\epsilon_t$) — one value per time step, computed separately for each spatial or observational unit. For vegetation indices, this is sometimes restricted to growing-season months only, since dormant-season values carry little information about vegetation stress.

**Climate correction** removes shared climate forcing (e.g. two climate covariates $x_1, x_2$) so that $\epsilon_t$ reflects only the biological/structural response:

1. Fit a linear model on an early baseline period using `lm()` — R's ordinary least squares function, which finds the coefficients $\hat{\beta}$ that minimize the sum of squared errors between observed and predicted values. This estimates $\hat{\beta}_0, \hat{\beta}_1, \hat{\beta}_2$ that best explain $y$'s variation from climate during that baseline period:

$$\hat{y}_t = \hat{\beta}_0 + \hat{\beta}_1 \cdot x_{1,t} + \hat{\beta}_2 \cdot x_{2,t}$$

**Illustrative example** — 4 months, 2 predictors ($x_1$ and $x_2$):

OLS (Ordinary Least Squares) is only meaningful when observations > unknowns. Here we have 4 months and 3 unknowns ($\hat{\beta}_0, \hat{\beta}_1, \hat{\beta}_2$) — no perfect fit exists, so `lm()` finds the coefficients that minimize the sum of squared residuals. A real application would use far more observations than unknowns; 4 is used here purely to keep the arithmetic tractable by hand.

| Month | $x_1$ | $x_2$ | $y$ |
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

**Fitted model:** $\hat{y} = 0.628 - 0.133 \cdot x_1 - 0.00058 \cdot x_2$

**Predictions and residuals:**

| Month | $\hat{y}$ | $y$ | residual |
|-------|-----------|-----|----------|
| 1 | $0.628 - 0.133(0.8) - 0.00058(120) = 0.452$ | 0.45 | −0.002 |
| 2 | $0.628 - 0.133(1.2) - 0.00058(80)  = 0.422$ | 0.42 | −0.002 |
| 3 | $0.628 - 0.133(1.6) - 0.00058(60)  = 0.380$ | 0.38 | −0.000 |
| 4 | $0.628 - 0.133(1.0) - 0.00058(100) = 0.437$ | 0.44 | +0.003 |

No month fits perfectly — `lm()` chose these coefficients because they minimize $(-0.002)^2 + (-0.002)^2 + (0.000)^2 + (0.003)^2 = 0.000017$. Any other coefficients would give a larger total.

2. Using those same coefficients from step 1, plug in each time step's observed $x_1$ and $x_2$ to get the climate-predicted $y$ for every time step in the full series.

3. Subtract to isolate the residual:

$$\epsilon_t = y_t - \hat{\beta}_0 - \hat{\beta}_1 \cdot x_{1,t} - \hat{\beta}_2 \cdot x_{2,t}$$

Training only on the early baseline period ensures the coefficients capture the normal climate–response relationship, uncontaminated by whatever anomaly the analysis is trying to detect.

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

$\sigma^2 = 0.000708$ is in squared units of the original series, which is hard to interpret directly. Taking the square root gives the standard deviation $\sigma = \sqrt{0.000708} \approx 0.027$ — a typical month in this window sits about 0.027 units above or below the window mean. As the collapse deepens and residuals swing further from the mean, $\sigma^2$ rises.

**Why $w-1$ not $w$:** dividing by $w-1$ (Bessel's correction) gives an unbiased estimate of true variance from a sample.

$\sigma^2$ is assigned to the **last month of the window**. Same rule as $r_1$.

**Rising $\sigma^2$** = residuals swinging wider = system less stable, being pushed further from equilibrium.

---

### 3. Skewness ($\gamma$)

**Question:** is the distribution of residuals symmetric, or lopsided to one side?

$$\gamma = \frac{\frac{1}{w}\sum(\epsilon_t - \mu)^3}{\left(\frac{1}{w}\sum(\epsilon_t - \mu)^2\right)^{3/2}}$$

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

$$\gamma = \frac{-0.00000546}{0.0000135} \approx \mathbf{-0.40}$$

Negative skew: the most extreme deviation in this window (month 2, −0.037) sits on the low side, pulling the distribution's tail to the left — a minority of unusually low readings against an otherwise fairly symmetric set of values. Assigned to the **last month of the window**, same rule as $r_1$ and $\sigma^2$.

| $\gamma$ value | Meaning |
|----------------|---------|
| $\approx 0$ | Symmetric — no side dominates |
| Negative | Left-skewed — a minority of unusually low readings pull the tail down |
| Positive | Right-skewed — a minority of unusually high readings pull the tail up |

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

---

## One-Way ANOVA

**Question:** do three or more groups have genuinely different means, or could the differences we see just be random noise?

ANOVA (Analysis of Variance) compares variation **between** groups to variation **within** groups:

$$F = \frac{\text{between-group variance}}{\text{within-group variance}} = \frac{MS_{between}}{MS_{within}}$$

($MS$ = Mean Square — a sum of squares divided by its degrees of freedom, i.e. an average squared deviation.)

If the groups really differ, the average distance between group means should be large relative to how spread out each group is internally. If the "groups" are really just random subsets of the same population, between-group and within-group variance should be about the same size, giving $F \approx 1$.

**Illustrative example** — 3 groups, 3 observations each:

| Group | Values | Mean |
|---|---|---|
| 1 | 0.10, 0.12, 0.14 | 0.12 |
| 2 | 0.14, 0.15, 0.16 | 0.15 |
| 3 | 0.20, 0.22, 0.24 | 0.22 |

Grand mean = (0.12+0.15+0.22)/3 = 0.163

**Between-group sum of squares** ($SS$ = Sum of Squares; each group's deviation from the grand mean, weighted by group size $n_k=3$):

$$SS_{between} = 3(0.12-0.163)^2 + 3(0.15-0.163)^2 + 3(0.22-0.163)^2 = 0.00563+0.00053+0.00963 = 0.0158$$

Degrees of freedom between = $k-1 = 3-1=2$, so $MS_{between} = 0.0158/2 = 0.0079$

**Within-group sum of squares** (spread inside each group around its own mean — same squared-deviation idea as $SS_{between}$, just measured from each group's *own* mean instead of the grand mean):

Group 1 (values 0.10, 0.12, 0.14; mean 0.12): $(0.10-0.12)^2+(0.12-0.12)^2+(0.14-0.12)^2 = 0.0004+0+0.0004 = 0.0008$

Group 2 (values 0.14, 0.15, 0.16; mean 0.15): $(0.14-0.15)^2+(0.15-0.15)^2+(0.16-0.15)^2 = 0.0001+0+0.0001 = 0.0002$

Group 3 (values 0.20, 0.22, 0.24; mean 0.22): $(0.20-0.22)^2+(0.22-0.22)^2+(0.24-0.22)^2 = 0.0004+0+0.0004 = 0.0008$

$$SS_{within} = \underbrace{0.0008}_{\text{group 1}} + \underbrace{0.0002}_{\text{group 2}} + \underbrace{0.0008}_{\text{group 3}} = 0.0018$$

Degrees of freedom within = $N-k = 9-3=6$ (9 total observations, minus 1 "used up" per group to compute that group's own mean), so $MS_{within} = 0.0018/6 = 0.0003$

**F-statistic:**

$$F = \frac{0.0079}{0.0003} \approx 26.3$$

**Where the p-value comes from:** F follows a known probability distribution when there's truly no group difference (random noise only), shaped by two numbers, $df_{between}$ and $df_{within}$ — here 2 and 6. The p-value is the probability of landing at or beyond the observed F under that distribution.

When $df_{between}=2$ (as it is here — always true for 3 groups, since $k-1=2$), that probability has a clean closed form:

$$p = \left(\frac{df_{within}}{df_{within}+2F}\right)^{df_{within}/2}$$

$$= \left(\frac{6}{6+2(26.33)}\right)^{3} = \left(\frac{6}{58.66}\right)^3 = (0.1023)^3 \approx 0.00107$$

Verified: this hand formula gives 0.001070108, matching both R's `pf(26.33, 2, 6, lower.tail=FALSE)` and `aov()`'s reported p-value exactly.

A large F means the between-group spread dwarfs the within-group spread — strong evidence the groups really differ. The p-value asks: *if there were actually no difference between groups, how likely is an F this large by chance alone?* A small p-value (< 0.05, conventionally) means: reject the "no difference" hypothesis.

| F value | Interpretation |
|---|---|
| Close to 1 | Between-group variation ≈ within-group variation — groups probably don't really differ |
| Much greater than 1 | Between-group variation dwarfs within-group — groups likely genuinely differ |

**Caveat:** ANOVA's F-test only tells you *some* pair of groups differs — not *which* pair. That's what Tukey HSD (below) is for.

---

## Tukey HSD (Honest Significant Difference)

**Question:** ANOVA said *some* pair of groups differs — but which pair(s), specifically?

If you had 3 groups and just ran 3 separate t-tests (1 vs 2, 2 vs 3, 1 vs 3), each test has its own ~5% chance of a false positive. Run enough pairwise tests and you're likely to find "significant" differences purely by chance — the **multiple comparisons problem**. Tukey HSD corrects for this: it controls the **family-wise error rate** (the chance of *any* false positive across *all* pairwise comparisons combined), not just each comparison individually.

For each pair of group means, Tukey HSD computes:

$$q = \frac{|\bar{y}_i - \bar{y}_j|}{SE}, \qquad SE = \sqrt{\frac{MS_{within}}{2}\left(\frac{1}{n_i}+\frac{1}{n_j}\right)}$$

and compares $q$ to the **studentized range distribution** — a distribution calibrated for "the largest gap among $k$ group means," which is inherently wider than a single pairwise comparison's distribution, precisely because there are more chances to see a big gap by luck when comparing many pairs at once.

**Illustrative example** — using the ANOVA numbers above ($MS_{within}=0.0003$, $n_i=n_j=3$ for every group):

$$SE = \sqrt{\frac{0.0003}{2}\left(\frac{1}{3}+\frac{1}{3}\right)} = \sqrt{0.0003 \times \tfrac{1}{3}} = \sqrt{0.0001} = 0.01$$

Group 1 vs Group 3: $|0.12-0.22| = 0.10 \to q = 0.10/0.01 = 10.0$ — compared against the studentized range critical value for $k=3$ groups and $df_{within}=6$ ($q_{crit} \approx 4.34$ at $\alpha=0.05$). Since $10.0 > 4.34$, this pair is significantly different.

Tukey HSD returns, for every pair: the mean difference, a confidence interval, and an adjusted p-value — already corrected for testing multiple pairs, so the p-value can be read directly without further correction.

---

## ROC Curve and AUC

**Question:** using a continuous score to predict a yes/no outcome, how good is that score at separating the two groups — across every possible cutoff, not just one?

Any continuous score can be turned into a binary classifier by picking a threshold: "predict positive if score ≥ threshold." Different thresholds trade off two error types:

- **True Positive Rate (TPR / Sensitivity):** of the observations that are truly positive, what fraction did we correctly flag? $TPR = TP/(TP+FN)$
- **False Positive Rate (FPR):** of the observations that are truly negative, what fraction did we incorrectly flag anyway? $FPR = FP/(FP+TN)$

The ROC curve plots TPR (y-axis) against FPR (x-axis) as the threshold sweeps from strict (flag almost nothing) to loose (flag almost everything). AUC (Area Under the Curve) condenses that whole curve into one number: the probability that a randomly-chosen truly-positive observation gets a higher score than a randomly-chosen truly-negative observation.

**Illustrative example** — 4 observations, sorted by score, descending:

| Observation | Score | Truly positive? |
|---|---|---|
| A | 0.9 | yes |
| B | 0.7 | no |
| C | 0.5 | yes |
| D | 0.2 | no |

2 positive observations (A, C), 2 negative (B, D). Sweeping the threshold down through each score and recomputing *cumulative* TPR/FPR (cumulative counts as each observation gets included, in descending score order, no synthetic (0,0) anchor point):

| Observation included | TP | FP | TPR | FPR |
|---|---|---|---|---|
| A | 1 | 0 | 0.5 | 0.0 |
| A, B | 1 | 1 | 0.5 | 0.5 |
| A, B, C | 2 | 1 | 1.0 | 0.5 |
| A, B, C, D | 2 | 2 | 1.0 | 1.0 |

AUC via the trapezoidal rule over these 4 (FPR, TPR) points:

$$AUC = \sum \Delta FPR \times \frac{TPR_i + TPR_{i+1}}{2}$$

$$= (0.5-0)\times\frac{0.5+0.5}{2} + (0.5-0.5)\times\frac{0.5+1.0}{2} + (1.0-0.5)\times\frac{1.0+1.0}{2}$$

$$= 0.25 + 0 + 0.5 = \mathbf{0.75}$$

**Cross-check via the probability definition:** of all 4 (positive, negative) pairs — (A,B), (A,D), (C,B), (C,D) — the positive observation scores higher in 3 of them (only C=0.5 loses to B=0.7). $3/4 = 0.75$ — matches exactly, as it must (the trapezoidal AUC and the "probability of correct ranking" are mathematically the same quantity).

| AUC | Meaning |
|---|---|
| 0.5 | No better than a coin flip — the score carries no information about the outcome |
| 1.0 | Perfect separation — every positive observation scores higher than every negative observation |
| 0.5–0.7 | Weak-to-modest discrimination |
| 0.7+ | Reasonably strong discrimination |

---

## Choosing K: Elbow Method and Silhouette Score

K-means requires deciding the number of clusters, K, in advance. Two diagnostic tools help pick a reasonable value by running K-means repeatedly across a range of K (e.g., K = 2 through 8) and comparing how each K performs.

### Elbow Method (WCSS)

**Question:** as K increases, how much does within-cluster spread keep shrinking — and where do the improvements stop being worth the added complexity?

WCSS = Within-Cluster Sum of Squares — for each point, the squared distance to its own cluster's center, summed over all points:

$$WCSS = \sum_{k=1}^{K}\sum_{i \in \text{cluster } k} \lVert x_i - c_k \rVert^2$$

WCSS *always* decreases as K increases (more clusters can only fit the data at least as well — in the extreme, K = N gives WCSS = 0, one cluster per point). The question isn't "which K minimizes WCSS" — more K always wins that trivially — it's "where does adding another cluster stop buying much." Plotted against K, WCSS drops steeply at first, then flattens into a shallow, diminishing-returns tail. The bend where the steep drop meets the flat tail is the "elbow" — a visual, somewhat subjective cue for a reasonable K.

**Illustrative example** — 6 points with an obvious two-cluster structure: $\{1, 2, 3\}$ and $\{8, 9, 10\}$.

$K=1$ (one cluster, center = mean of all 6 = 5.5):
$$WCSS = (1-5.5)^2+(2-5.5)^2+(3-5.5)^2+(8-5.5)^2+(9-5.5)^2+(10-5.5)^2 = 77.5$$

$K=2$ (natural split, centers 2 and 9):
$$WCSS = \underbrace{[(1-2)^2+(2-2)^2+(3-2)^2]}_{2} + \underbrace{[(8-9)^2+(9-9)^2+(10-9)^2]}_{2} = 4.0$$

$K=3$ (e.g. $\{1,2\}, \{3\}, \{8,9,10\}$, centers 1.5, 3, 9):
$$WCSS = \underbrace{[(1-1.5)^2+(2-1.5)^2]}_{0.5} + \underbrace{0}_{\{3\}} + \underbrace{[(8-9)^2+(9-9)^2+(10-9)^2]}_{2} = 2.5$$

| K | WCSS | Drop from previous K |
|---|---|---|
| 1 | 77.5 | — |
| 2 | 4.0 | −73.5 (steep) |
| 3 | 2.5 | −1.5 (small) |

The huge drop from $K=1\to2$ (77.5 → 4.0) reflects the data's real two-cluster structure; the much smaller drop from $K=2\to3$ (4.0 → 2.5) is just fitting noise inside already-good clusters. The elbow is the bend at $K=2$, where the curve flattens.

### Silhouette Score

**Question:** for a given K, how well-separated are the resulting clusters — is each point clearly closer to its own cluster than to the nearest other one?

For each point $i$, define:
- $a(i)$ = average distance from $i$ to every other point in its **own** cluster (how tightly it fits where it was placed)
- $b(i)$ = average distance from $i$ to every point in the **nearest other** cluster (how close the next-best option is)

$$s(i) = \frac{b(i)-a(i)}{\max(a(i),b(i))}$$

The **silhouette score** for a given K is the mean of $s(i)$ across all points.

| $s(i)$ | Meaning |
|---|---|
| Close to +1 | Point fits its own cluster well, is far from the next-nearest cluster — good clustering |
| Close to 0 | Point sits right on the boundary between two clusters — ambiguous |
| Close to −1 | Point is, on average, closer to a *different* cluster than its own — likely mis-clustered |

**Illustrative example** — 1 point, 2 candidate clusters:

Point $i$ is at distance 1, 1, 2 from the other 3 points in its own cluster (own-cluster distances), and distance 4, 5 from the 2 points in the nearest other cluster:

$$a(i) = \frac{1+1+2}{3} = 1.33 \qquad b(i) = \frac{4+5}{2} = 4.5$$

$$s(i) = \frac{4.5-1.33}{\max(1.33,4.5)} = \frac{3.17}{4.5} \approx 0.70$$

A silhouette of 0.70 for this point indicates it's well-placed — solidly closer to its own cluster than to the next-nearest one.

**Why both, and why they can disagree:** WCSS keeps dropping no matter what (favors bigger K), while silhouette explicitly penalizes clusters that are close together / overlapping (can favor smaller K). When they disagree, that's real information — the "statistically optimal" K by silhouette isn't automatically the most *useful* or *interpretable* K. For a large dataset, silhouette is often computed on a random subsample rather than every point, since it requires a full pairwise distance matrix ($O(n^2)$).

