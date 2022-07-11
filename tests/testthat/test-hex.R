test_that("equilibrium frequencies same as in hwep", {
  yww <- c(29/120, 21/120, 17/120, 10/120, 10/120, 10/120, 23/120)
  rho <- sum(yww * 0:6) / 6
  cgeno <- hex_recursive(yww = yww, niter = 30, alpha = 0.1)
  hgeno <- hwep::hwefreq(r = rho, alpha = 0.1, ploidy = 6)
  expect_equal(cgeno, hgeno, ignore_attr = TRUE)
})


test_that("EM estimates are same as in hwep", {
  nvec <- c(29, 21, 17, 10, 10, 10, 23)
  yww <- nvec / sum(nvec)
  hout <- hex_em(yww = yww, niter = 1000)
  lout <- hwep::rmlike(nvec = nvec, thresh = 0)
  expect_equal(lout$p, hout$p, tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(stats::convolve(hout$p, rev(hout$p), type = "open"), hout$q)
})

