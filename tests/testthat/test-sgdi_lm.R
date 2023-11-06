test_that("sgdi_lm gamma_0 matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, gamma_0=0.1)
  out2 = sgdi_lm(y~., data=my.dat, gamma_0=10)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_lm vs. sgd_lm", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  bt_start = bt0 * rnorm(5, mean=1, sd=0.25)
  out1 = sgd_lm(y~., data=my.dat, bt_start=bt_start)
  out2 = sgdi_lm(y~., data=my.dat, bt_start=bt_start)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_true(check==0)
})

test_that("Stop when rss_idx includes 0 with inference=rss", 
  {
    n = 1e05
    p = 5
    bt0 = rep(5,p)
    x = matrix(rnorm(n*(p-1)), n, (p-1))
    y = cbind(1,x) %*% bt0 + rnorm(n)
    my.dat = data.frame(y=y, x=x)
    expect_error(sgdi_lm(y~., data=my.dat, inference="rss", rss_idx=c(0,1)))
  }
)

test_that("Run the code with an option studentize=F", 
          {
            n = 1e05
            p = 5
            bt0 = rep(5,p)
            x = matrix(rnorm(n*(p-1)), n, (p-1))
            y = cbind(1,x) %*% bt0 + rnorm(n)
            my.dat = data.frame(y=y, x=x)
            out1 = sgdi_lm(y~., data=my.dat, studentize=F)
            out2 = sgdi_lm(y~., data=my.dat, studentize=T)
            check = max(abs(out1$coefficients - out2$coefficients))
            expect_true(check<1e-2)
          }
)

test_that("sgdi_lm studentize matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, studentize=TRUE)
  out2 = sgdi_lm(y~., data=my.dat, studentize=FALSE)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

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

test_that("sgdi_lm no_studentize matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat)
  out2 = sgdi_lm(y~., data=my.dat, no_studentize=1e06)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_lm with scalar X", {
  n = 1e05
  p = 2
  bt0 = rep(5,2)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, studentize=TRUE)
  out2 = sgdi_lm(y~., data=my.dat, studentize=FALSE)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_lm with scalar X with/out intercept", {
  n = 1e05
  p = 2
  bt0 = c(0,5)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, intercept=TRUE)
  out2 = sgdi_lm(y~., data=my.dat, intercept=FALSE)
  check = max(abs(out1$coefficients[-1] - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_lm burn matters 1", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat)
  out2 = sgdi_lm(y~., data=my.dat, burn=100)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_lm burn matters 2", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, studentize=FALSE)
  out2 = sgdi_lm(y~., data=my.dat, burn=100, studentize=FALSE)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})
