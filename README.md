# homa2calc

An R package for computing **HOMA2 %B**, **HOMA2 %S**, and **HOMA2-IR** from the
Oxford Diabetes Trials Unit model of steady-state beta-cell function and insulin
resistance.

## Background

HOMA2 (Homeostasis Model Assessment 2) is a mathematical model of fasting
glucose-insulin homeostasis developed by the Diabetes Trials Unit, University
of Oxford. This package implements HOMA2 via bilinear interpolation over
reference tables extracted from the **HOMA2 Calculator v2.2.4 validation
dataset**. Grid-point values match the official DTU calculator exactly.

## Installation

```r
# Install from GitHub
remotes::install_github("ZanettaT/homa2calc")
```

## Usage

### Three assay modes

| Function | Hormone input | Units | Valid range |
|---|---|---|---|
| `homa2_insulin()` | Non-specific insulin (cross-reacts with proinsulin) | pmol/L | 20–400 |
| `homa2_specific_insulin()` | Proinsulin-free (specific) insulin | pmol/L | 20–300 |
| `homa2_cpeptide()` | Fasting C-peptide | nmol/L | 0.2–3.5 |

All functions take **glucose in mmol/L** (valid range: 3.0–25.0).

### Unit conversions

| From | To (pmol/L or nmol/L) | Formula |
|---|---|---|
| Insulin µU/mL (conventional RIA) | pmol/L | × 6.0 |
| Insulin µU/mL (proinsulin-free assay) | pmol/L | × 6.0 → use `homa2_specific_insulin()` |
| Glucose mg/dL | mmol/L | ÷ 18.0 |
| C-peptide ng/mL | nmol/L | ÷ 3.02 |
| C-peptide pmol/L | nmol/L | ÷ 1000 |

### Examples

```r
library(homa2calc)

# Single observation
homa2_insulin(glucose = 5.0, insulin = 60)
#   homa2_b homa2_s homa2_ir
# 1   105.6    90.1    1.109

# i000 is in µU/mL; convert to pmol/L by × 6 for standard RIA
dpp <- dpp %>%
  mutate(
    ins_pmol = i000 * 6.0,
    homa2    = homa2_insulin(g000_mmol, ins_pmol)
  ) %>%
  mutate(
    homa2_ir    = homa2$homa2_ir,
    homa2_b_ins = homa2$homa2_b
  )

# With C-peptide (if available; convert from ng/mL to nmol/L)
dpp <- dpp %>%
  mutate(
    c0_nmol   = c0_ng_ml / 3.02,
    homa2_cpep = homa2_cpeptide(g000_mmol, c0_nmol)
  ) %>%
  mutate(homa2_b_cpep = homa2_cpep$homa2_b)
```

### Output columns

| Column | Description |
|---|---|
| `homa2_b` | Beta-cell function (%)|
| `homa2_s` | Insulin sensitivity (%)|
| `homa2_ir` | Insulin resistance index = 100 / %S|

Values outside the valid input range are returned as `NA`.

## Validation

At grid-point values, this package reproduces the official DTU calculator output
exactly. Interpolated values are accurate to within rounding error (~0.1 units).

```r
# Exact grid-point match examples (HOMA2 v2.2.4 expected values)
homa2_insulin(3.0, 20)   # %B=139.8, %S=308.9, IR=0.3237
homa2_insulin(14.0, 96)  # %B=21.7,  %S=42.7,  IR=2.3419
homa2_cpeptide(3.0, 0.2) # %B=151.8, %S=272.6
```

## References
- HOMA2 Calculator v2.2.4, Diabetes Trials Unit, University of Oxford.
  <https://www.dtu.ox.ac.uk/homacalculator/>
