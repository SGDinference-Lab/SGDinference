########################################################
### There is an error message in quantile regression ###
########################################################
test_that("sgdi_qr should work with option rs", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_qr(y~., data=my.dat, inference="rs")
  expect_error(print(out), NA)
})
test_that("sgdi_qr should work with option rss", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_qr(y~., data=my.dat, inference="rss")
  expect_error(print(out), NA)
})
test_that("sgdi_qr should work with option rsd", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_qr(y~., data=my.dat, inference="rsd")
  expect_error(print(out), NA)
})
###############################################
### There is no problem for mean regression ###
###############################################
test_that("sgdi_lm should work with option rs", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_lm(y~., data=my.dat, inference="rs")
  expect_error(print(out), NA)
})
test_that("sgdi_lm should work with option rss", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_lm(y~., data=my.dat, inference="rss")
  expect_error(print(out), NA)
})
test_that("sgdi_lm should work with option rsd", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgdi_lm(y~., data=my.dat, inference="rsd")
  expect_error(print(out), NA)
})

