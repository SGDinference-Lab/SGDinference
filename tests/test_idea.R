library(microbenchmark)
library(ggplot2)
library(SGDinference)

n = 1e05
p = 5
bt0 = rep(5,p)
set.seed(22023)

x = matrix(rnorm(n*(p-1)), n, (p-1))
y = cbind(1,x) %*% bt0 + rnorm(n)
my.dat = data.frame(y=y, x=x)

my.test = microbenchmark(
  sgdi.out.rs = sgdi_lm(y~., data=my.dat, inference="rs", studentize=F),
  sgdi.out.rs.new = sgdi_lm_new(y~., data=my.dat, inference="rs", studentize=F),
  times = 100L
)
print(my.test)
autoplot(my.test)
