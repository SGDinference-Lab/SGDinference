###########################
### Quantile Regression ###
###########################
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
#############################################
### Quantile Regression without Inference ###
#############################################
## NB: Without inference, we do not choose any option on V.hat estimation.
test_that("sgd_qr should work", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgd_qr(y~., data=my.dat, path=F)
  expect_error(print(out), NA)
})

#######################
### Mean Regression ###
#######################
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
test_that("sgdi_lm: bt_start matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, bt_start=NULL)
  out2 = sgdi_lm(y~., data=my.dat, bt_start=bt0)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

#########################################
### Mean Regression without inference ###
#########################################
test_that("sgd_lm should work", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out = sgd_lm(y~., data=my.dat)
  expect_error(print(out), NA)
})
