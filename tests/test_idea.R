library(SGDinference)

n = 1e05
p = 5
bt0 = rep(5,p)
set.seed(22023)

x = matrix(rnorm(n*(p-1)), n, (p-1))
y = cbind(1,x) %*% bt0 + rnorm(n)
my.dat = data.frame(y=y, x=x)

sgdi.out.rs = sgdi_lm(y~., data=my.dat, inference="rs", studentize=T)
sgdi.out.rs
sgdi.out.rsd = sgdi_lm(y~., data=my.dat, inference="rsd", studentize=T)

sgdi.out.rss = sgdi_lm(y~., data=my.dat, inference="rss", rss_idx = c(1,3), studentize=T)
sgdi.out.rss

print(my.test)
autoplot(my.test)
