test_that("expect no error using sgdi_lm", {
  y = Census2000$ln_hrwage 
  edu = Census2000$edyrs
  exp = Census2000$exp
  exp2 = exp^2
  x = cbind(edu, exp, exp2)
  out = sgdi_lm(y ~ x)
  expect_no_error(out)
})

test_that("expect difference between the 10th and 90th quantile regression estimates", {
  y = Census2000$ln_hrwage 
  edu = Census2000$edyrs
  exp = Census2000$exp
  exp2 = exp^2
  x = cbind(edu, exp, exp2)
  out1 = sgdi_qr(y ~ x, qt=0.1)
  out2 = sgdi_qr(y ~ x, qt=0.2)
  check = max(abs(out1$coefficients-out2$coefficients))
  expect_false(check==0)
})
