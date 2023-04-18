test_that("sgdi vs sgdi_qr", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi(y~., data=my.dat, model="qr")
  out2 = sgdi_qr(y~., data=my.dat)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_true(check==0)
})

test_that("sgdi vs sgdi_lm", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi(y~., data=my.dat, model="lm")
  out2 = sgdi_lm(y~., data=my.dat)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_true(check==0)
})

