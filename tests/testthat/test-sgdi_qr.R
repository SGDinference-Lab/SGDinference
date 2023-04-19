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
