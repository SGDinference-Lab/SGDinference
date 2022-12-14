test_that("sgdi_lm using Census2000", {
    y = Census2000$ln_hrwage 
  edu = Census2000$edyrs
  exp = Census2000$exp
 exp2 = exp^2
    x = cbind(edu, exp, exp2)
 out1 = lm(y ~ x)
  bh1 = out1$coefficients[2]
 out2 = sgdi_lm(x,y,gamma_0=1)
  bh2 = out2$beta_hat[2]
  check = abs(bh1 - bh2)
  expect_true(check < 0.01)
})