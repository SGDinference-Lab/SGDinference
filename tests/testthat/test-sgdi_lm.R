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

test_that("sgdi_lm gamma_0 matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  out1 = sgdi_lm(x,y,gamma_0=0.1)
  out2 = sgdi_lm(x,y,gamma_0=10)
  check = max(abs(out1$beta_hat - out2$beta_hat))
  expect_false(check==0)
})

test_that("sgdi_lm vs. sgd_lm", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  out1 = sgd_lm(x,y)
  out2 = sgdi_lm(x,y)
  check = max(abs(out1$beta_hat - out2$beta_hat))
  expect_true(check==0)
})

