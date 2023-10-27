test_that("expect no error using sgdi_lm", {
  y = Census2000$ln_hrwage 
  edu = Census2000$edyrs
  exp = Census2000$exp
  exp2 = exp^2
  x = cbind(edu, exp, exp2)
  out = sgdi_lm(y ~ x)
  expect_no_error(out)
})

test_that("expect no error using sgdi", {
  y = Census2000$ln_hrwage 
  edu = Census2000$edyrs
  exp = Census2000$exp
  exp2 = exp^2
  x = cbind(edu, exp, exp2)
  out = sgdi(y ~ x)
#  my.dat = data.frame(y,x)
#  out = sgdi(y ~ x, data=my.dat)
  expect_no_error(out)
# Currently, we have the following error message:
# Error in `eval(mf, parent.frame())`: argument "data" is missing, with no default  
# Note: sgdi_lm works but sgdi does not  
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

test_that("expect no error using sgdi", {
  y = Census2000$ln_hrwage 
  edu = Census2000$edyrs
  exp = Census2000$exp
  exp2 = exp^2
  x = cbind(edu, exp, exp2)
  out = sgdi_boot_qr(y ~ edu + exp + exp2, n_boot=100)
  expect_no_error(out)
  # Error in `dim(data) <- dim`:
  # ! invalid first argument, must be vector (list or atomic)
  #Backtrace:
  #  1. SGDinference::sgdi_boot_qr(y ~ edu + exp + exp2, n_boot = 100)
  #  3. base::as.matrix.default(x)
  #  4. base::array(...)
  #Execution halted
})
