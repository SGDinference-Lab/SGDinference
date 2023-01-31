test_that("sgdi_qr gamma_0 matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_qr(y~., data=my.dat, gamma_0=0.1)
  check = max(abs(out$coefficient - bt0))
  expect_true(check > 1)
})

test_that("sgdi_qr vs. sgd_qr", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_qr(y~., data=my.dat)
  out2 = sgdi_qr(y~., data=my.dat)
  check = max(abs(out1$coefficient - out2$coefficient))
  expect_true(check==0)
})

