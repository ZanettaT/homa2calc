library(testthat)
source("../../R/tables.R")
source("../../R/homa2.R")

# ── Insulin mode: exact grid-point checks ──────────────────────────────────
test_that("homa2_insulin matches reference at grid points", {
  r <- homa2_insulin(3.0, 20)
  expect_equal(r$homa2_b, 139.8)
  expect_equal(r$homa2_s, 308.9)

  r2 <- homa2_insulin(5.2, 58)
  expect_equal(r2$homa2_b, 93.0)
  expect_equal(r2$homa2_s, 91.4)

  r3 <- homa2_insulin(14.0, 96)
  expect_equal(r3$homa2_b, 21.7)
  expect_equal(r3$homa2_s, 42.7)
})

# ── C-peptide mode: exact grid-point checks ────────────────────────────────
test_that("homa2_cpeptide matches reference at grid points", {
  r <- homa2_cpeptide(3.0, 0.2)
  expect_equal(r$homa2_b, 151.8)
  expect_equal(r$homa2_s, 272.6)

  r2 <- homa2_cpeptide(14.0, 3.5)
  expect_equal(r2$homa2_b, 97.2)
  expect_equal(r2$homa2_s, 8.9)
})

# ── Specific insulin mode ──────────────────────────────────────────────────
test_that("homa2_specific_insulin matches reference at grid points", {
  r <- homa2_specific_insulin(3.0, 20)
  expect_equal(r$homa2_b, 150.3)
  expect_equal(r$homa2_s, 276.9)
})

# ── Out-of-range returns NA ────────────────────────────────────────────────
test_that("out-of-range inputs return NA", {
  r <- homa2_insulin(2.0, 50)
  expect_true(is.na(r$homa2_ir))

  r2 <- homa2_insulin(5.0, 500)
  expect_true(is.na(r2$homa2_ir))
})

# ── Vectorised output ──────────────────────────────────────────────────────
test_that("vectorised output has correct structure", {
  res <- homa2_insulin(c(4.0, 6.0, 8.5), c(60, 80, 120))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_named(res, c("homa2_b", "homa2_s", "homa2_ir"))
})

# ── IR = 100 / %S ──────────────────────────────────────────────────────────
test_that("homa2_ir = 100 / homa2_s", {
  r <- homa2_insulin(6.0, 80)
  expect_equal(r$homa2_ir, 100 / r$homa2_s, tolerance = 1e-4)
})
