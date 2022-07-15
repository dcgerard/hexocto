test_that("octoploid segregation frequencies are same as in Wang et al (2021)", {
p <- c(1, 0, 0, 0, 0)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(1, 0, 0, 0, 0, 0, 0, 0, 0), alpha = 0)
)

p <- c(9/16, 3/8, 1/16, 0, 0)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 1, 0, 0, 0, 0, 0, 0, 0), alpha = 0)
)

p <- c(225/784, 45/98, 87/392, 3/98, 1/784)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 1, 0, 0, 0, 0, 0, 0), alpha = 0)
)

p <- c(25/196, 75/196, 285/784, 45/392, 9/784)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 0, 1, 0, 0, 0, 0, 0), alpha = 0)
)

p <- c(9/196, 12/49, 41/98, 12/49, 9/196)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 0, 0, 1, 0, 0, 0, 0), alpha = 0)
)

p <- c(9/784, 45/392, 285/784, 75/196, 25/196)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 0, 0, 0, 1, 0, 0, 0), alpha = 0)
)

p <- c(1/784, 3/98, 87/392, 45/98, 225/784)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 0, 0, 0, 0, 1, 0, 0), alpha = 0)
)

p <- c(0, 0, 1/16, 3/8, 9/16)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 0, 0, 0, 0, 0, 1, 0), alpha = 0)
)

p <- c(0, 0, 0, 0, 1)
expect_equal(
  convolve(p, rev(p), type = "open"),
  octo_onegen(yww = c(0, 0, 0, 0, 0, 0, 0, 0, 1), alpha = 0)
)

})
