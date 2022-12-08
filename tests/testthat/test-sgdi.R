test_that("sgdi vs sgdi_qr", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  out1 = sgdi(x,y, model="qr")
  out2 = sgdi_qr(x,y)
  check = max(abs(out1$beta_hat - out2$beta_hat))
  expect_true(check==0)
})

test_that("sgdi vs sgdi_lm", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  out1 = sgdi(x,y, model="lm")
  out2 = sgdi_qr(x,y)
  check = max(abs(out1$beta_hat - out2$beta_hat))
  expect_true(check==0)
})

