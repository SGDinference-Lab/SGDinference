test_that("sgdi_lm rss", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  sgdi.out = sgdi_lm(x,y)
  out = sgdi_lm(x,y,inference="rss")
  expect_identical(length(out$V_hat_sub), as.integer(1))
})
