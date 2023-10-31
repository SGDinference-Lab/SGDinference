test_that("sgdi_qr gamma_0 matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat, gamma_0=0.01)
  check = max(abs(out1$coefficients - bt0))
  expect_true(check > 1)
})

test_that("sgdi_qr vs. sgd_qr", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_qr(y~., data=my.dat, studentize = F)
  out2 = sgdi_qr(y~., data=my.dat, studentize = F)
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
            expect_error(sgdi_qr(y~., data=my.dat, inference="rss", rss_idx=c(0,1)))
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
            out1 = sgdi_qr(y~., data=my.dat, studentize=F)
            out2 = sgdi_qr(y~., data=my.dat, studentize=T)
            check = max(abs(out1$coefficients - out2$coefficients))
            expect_true(check<1e-1)
          }
)

test_that("Different significance levels", 
          {
            n = 1e05
            p = 5
            bt0 = rep(5,p)
            x = matrix(rnorm(n*(p-1)), n, (p-1))
            y = cbind(1,x) %*% bt0 + rnorm(n)
            my.dat = data.frame(y=y, x=x)
            out1 = sgdi_qr(y~., data=my.dat, level=0.95)
            out2 = sgdi_qr(y~., data=my.dat, level=0.90)
            check = out1$ci.upper[1]- - out2$ci.upper[1]
            expect_true(check>0)
          }
)

test_that("sgdi_qr no_studentize matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat)
  out2 = sgdi_qr(y~., data=my.dat, no_studentize=1e06)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_qr with scalar X", {
  n = 1e05
  p = 2
  bt0 = rep(5,2)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat, studentize=TRUE)
  out2 = sgdi_qr(y~., data=my.dat, studentize=FALSE)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_qr with scalar X with/out intercept", {
  n = 1e05
  p = 2
  bt0 = c(0,5)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat, intercept=TRUE)
  out2 = sgdi_qr(y~., data=my.dat, intercept=FALSE)
  check = max(abs(out1$coefficients[-1] - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_qr burn matters 1", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat)
  out2 = sgdi_qr(y~., data=my.dat, burn=100)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_qr burn matters 2", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = 2*matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat, studentize=FALSE)
  out2 = sgdi_qr(y~., data=my.dat, burn=100, studentize=FALSE)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})

test_that("sgdi_qr path_index matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat, path=TRUE, path_index=1)
  out2 = sgdi_qr(y~., data=my.dat, path=TRUE, path_index=2)
  check = max(abs(out1$beta_hat_path[100] - out2$beta_hat_path[100]))
  expect_false(check==0)
})

test_that("sgdi_qr bt_start matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_qr(y~., data=my.dat)
  out2 = sgdi_qr(y~., data=my.dat, bt_start=bt0)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_false(check==0)
})