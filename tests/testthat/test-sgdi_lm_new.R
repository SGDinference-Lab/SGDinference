test_that("sgdi_lm and sgdi_lm_new gives the same answer", {
  n = 1e06
  p = 5
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  out1 = sgdi_lm(y~., data=my.dat)
  out2 = sgdi_lm_new(y~., data=my.dat, no_studentize=n)
  check = max(abs(out1$coefficients - out2$coefficients))
  expect_true(check==0)
})


test_that("Is sgdi_lm_new faster than sgdi_lm", {
  library(microbenchmark)
  n = 1e06
  p = 100
  bt0 = rep(5,p)
  x = matrix(rnorm(n*(p-1)), n, (p-1))
  y = cbind(1,x) %*% bt0 + rnorm(n)
  my.dat = data.frame(y=y, x=x)
  mbm = microbenchmark(
  lm = lm(y~., data=my.dat),
  out = sgdi_lm(y~., data=my.dat, inference="rsd"),
  out.new = sgdi_lm_new(y~., data=my.dat, inference="rsd"),
  out.no.st = sgdi_lm(y~., data=my.dat, studentize = F, inference="rsd"),
  times = 5L
  )
)
  expect_true(check==0)
})
