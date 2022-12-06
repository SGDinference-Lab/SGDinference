test_that("sgdi_qr gamma_0 matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  out = sgdi_qr(x,y,gamma_0=0.1)
  check = max(abs(out$beta_hat - bt0))
  expect_true(check > 1)
})

test_that("sgdi_qr vs. sgd_qr", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  out1 = sgd_qr(x,y)
  out2 = sgdi_qr(x,y)
  check = max(abs(out1$beta_hat - out2$beta_hat))
  expect_true(check==0)
})

