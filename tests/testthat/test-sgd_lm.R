test_that("sgd_lm studentize matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_lm(y~., data=my.dat, studentize=TRUE)
  out2 = sgd_lm(y~., data=my.dat, studentize=FALSE)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgd_lm no_studentize matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_lm(y~., data=my.dat)
  out2 = sgd_lm(y~., data=my.dat, no_studentize=1e06)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgd_lm intercept matters 1", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_lm(y~., data=my.dat)
  out2 = sgd_lm(y~., data=my.dat, intercept=FALSE)
  check = max(abs(out1$coefficients[-1] - out2$coefficients))
  expect_false(check==0)
})

test_that("sgd_lm intercept matters 2", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_lm(y~., data=my.dat, studentize=FALSE)
  out2 = sgd_lm(y~., data=my.dat, studentize=FALSE, intercept=FALSE)
  check = max(abs(out1$coefficientscoefficients[-1] - out2$coefficients))
  expect_false(check==0)
})

test_that("sgd_lm bt_start matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_lm(y~., data=my.dat)
  out2 = sgd_lm(y~., data=my.dat, bt_start=bt0)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})
