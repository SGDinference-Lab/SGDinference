#test_that("sgdi_lm rss", {
#  n = 1e05
#  p = 5
#  bt0 = rep(5,p)
#  x = matrix(rnorm(n*(p-1)), n, (p-1))
#  y = cbind(1,x) %*% bt0 + rnorm(n)
#  my.dat = data.frame(y=y, x=x)
#  sgdi.out = sgdi_lm(y~., data=my.dat)
#  out = sgdi_lm(y~., data=my.dat,inference="rss")
#  expect_identical(length(out$V_hat_sub), as.integer(1))
#})

test_that("sgdi_lm gamma_0 matters", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat, gamma_0=0.1)
  out2 = sgdi_lm(y~., data=my.dat, gamma_0=10)
  check = max(abs(out1$coefficient - out2$coefficient))
  expect_false(check==0)
})

test_that("sgdi_lm vs. sgd_lm", {
  n = 1e05
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgd_lm(y~., data=my.dat)
  out2 = sgdi_lm(y~., data=my.dat)
  check = max(abs(out1$coefficient - out2$coefficient))
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