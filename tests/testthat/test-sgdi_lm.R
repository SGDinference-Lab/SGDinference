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
  out1 = sgd_lm(y~., data=my.dat)
  out2 = sgdi_lm(y~., data=my.dat)
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
