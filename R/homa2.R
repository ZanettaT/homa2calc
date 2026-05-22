#' HOMA2 Calculator
#'
#' Computes HOMA2 %B (beta-cell function), HOMA2 %S (insulin sensitivity),
#' and HOMA2-IR (insulin resistance) using bilinear interpolation over
#' reference tables derived from the HOMA2 Calculator v2.2.4 validation
#' dataset (Diabetes Trials Unit, University of Oxford).
#'
#' @references
#' Levy JC, Matthews DR, Hermans MP. Correct homeostasis model assessment
#' (HOMA) evaluation uses the computer program. Diabetes Care. 1998;21(12):2191-2.
#'
#' Wallace TM, Levy JC, Matthews DR. Use and abuse of HOMA modeling.
#' Diabetes Care. 2004;27(6):1487-95.
#'
#' HOMA2 Calculator v2.2.4, Diabetes Trials Unit, University of Oxford.
#' https://www.dtu.ox.ac.uk/homacalculator/

# ── Internal bilinear interpolator ─────────────────────────────────────────

.homa2_interp_single <- function(glucose, hormone, gluc_axis, horm_axis,
                                  b_mat, s_mat) {
  n_g <- length(gluc_axis)
  n_h <- length(horm_axis)

  # Out-of-range → NA
  if (is.na(glucose) || is.na(hormone) ||
      glucose < gluc_axis[1]  || glucose > gluc_axis[n_g] ||
      hormone < horm_axis[1]  || hormone > horm_axis[n_h]) {
    return(c(pct_b = NA_real_, pct_s = NA_real_, ir = NA_real_))
  }

  # Bounding indices (glucose axis)
  g1 <- max(which(gluc_axis <= glucose))
  g2 <- min(which(gluc_axis >= glucose))
  # Bounding indices (hormone axis)
  h1 <- max(which(horm_axis <= hormone))
  h2 <- min(which(horm_axis >= hormone))

  # Fractional positions
  tx <- if (g2 == g1) 0 else (glucose - gluc_axis[g1]) / (gluc_axis[g2] - gluc_axis[g1])
  ty <- if (h2 == h1) 0 else (hormone - horm_axis[h1]) / (horm_axis[h2] - horm_axis[h1])

  # Bilinear interpolation
  bilerp <- function(mat) {
    (1 - tx) * (1 - ty) * mat[g1, h1] +
      tx       * (1 - ty) * mat[g2, h1] +
      (1 - tx) *      ty  * mat[g1, h2] +
      tx       *      ty  * mat[g2, h2]
  }

  pct_b <- bilerp(b_mat)
  pct_s <- bilerp(s_mat)
  ir    <- 100 / pct_s

  c(pct_b = round(pct_b, 1), pct_s = round(pct_s, 1), ir = ir)
}


# ── Vectorised wrapper ──────────────────────────────────────────────────────

.homa2_vectorised <- function(glucose, hormone, gluc_axis, horm_axis,
                               b_mat, s_mat) {
  n <- length(glucose)
  out_b  <- numeric(n)
  out_s  <- numeric(n)
  out_ir <- numeric(n)

  for (i in seq_len(n)) {
    res <- .homa2_interp_single(glucose[i], hormone[i],
                                gluc_axis, horm_axis, b_mat, s_mat)
    out_b[i]  <- res["pct_b"]
    out_s[i]  <- res["pct_s"]
    out_ir[i] <- res["ir"]
  }

  data.frame(homa2_b = out_b, homa2_s = out_s, homa2_ir = out_ir)
}


# ── Public API ──────────────────────────────────────────────────────────────

#' HOMA2 using non-specific insulin (cross-reacts with proinsulin)
#'
#' @param glucose Fasting plasma glucose in **mmol/L** (valid range: 3.0–25.0).
#' @param insulin Fasting plasma insulin in **pmol/L** (valid range: 20–400).
#'   To convert from µU/mL multiply by 6.0 (conventional RIA) or 6.945
#'   (proinsulin-free assay — use \code{homa2_specific_insulin()} instead).
#' @return A data frame with columns \code{homa2_b} (%B), \code{homa2_s} (%S),
#'   and \code{homa2_ir} (IR = 100/%S). Values outside the valid range are NA.
#' @export
#' @examples
#' # Single observation
#' homa2_insulin(glucose = 5.0, insulin = 60)
#'
#' # Vectorised over a data frame (µU/mL → pmol/L, multiply by 6)
#' # dpp$homa2 <- homa2_insulin(dpp$g000_mmol, dpp$i000 * 6)
homa2_insulin <- function(glucose, insulin) {
  glucose <- as.numeric(glucose)
  insulin <- as.numeric(insulin)
  stopifnot(length(glucose) == length(insulin))
  .homa2_vectorised(glucose, insulin,
                    homa2_ins_gluc_axis, homa2_ins_horm_axis,
                    homa2_ins_b_mat,     homa2_ins_s_mat)
}


#' HOMA2 using specific (proinsulin-free) insulin
#'
#' @param glucose Fasting plasma glucose in **mmol/L** (valid range: 3.0–25.0).
#' @param specific_insulin Fasting specific insulin in **pmol/L**
#'   (valid range: 20–300).
#'   To convert from µU/mL with a specific assay multiply by 6.0.
#' @inheritParams homa2_insulin
#' @export
homa2_specific_insulin <- function(glucose, specific_insulin) {
  glucose          <- as.numeric(glucose)
  specific_insulin <- as.numeric(specific_insulin)
  stopifnot(length(glucose) == length(specific_insulin))
  .homa2_vectorised(glucose, specific_insulin,
                    homa2_spec_gluc_axis, homa2_spec_horm_axis,
                    homa2_spec_b_mat,     homa2_spec_s_mat)
}


#' HOMA2 using fasting C-peptide
#'
#' @param glucose Fasting plasma glucose in **mmol/L** (valid range: 3.0–25.0).
#' @param cpeptide Fasting C-peptide in **nmol/L** (valid range: 0.2–3.5).
#'   Common unit conversions:
#'   \itemize{
#'     \item ng/mL → nmol/L: divide by 3.02
#'     \item pmol/L → nmol/L: divide by 1000
#'   }
#' @inheritParams homa2_insulin
#' @export
homa2_cpeptide <- function(glucose, cpeptide) {
  glucose  <- as.numeric(glucose)
  cpeptide <- as.numeric(cpeptide)
  stopifnot(length(glucose) == length(cpeptide))
  .homa2_vectorised(glucose, cpeptide,
                    homa2_cpep_gluc_axis, homa2_cpep_horm_axis,
                    homa2_cpep_b_mat,     homa2_cpep_s_mat)
}
