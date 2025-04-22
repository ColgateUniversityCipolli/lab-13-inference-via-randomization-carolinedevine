################################################################################
# Lab 13: Final Lab
# Caroline Devine
# Due: 4/29
################################################################################
library(tidyverse)
library(xtable)
library(VGAM)
library(e1071)
library(pwr)
library(effectsize)
################################################################################
# Question 1
################################################################################

### Part A ### -> using further observations
dat.finches <- read.csv("zebrafinches.csv")

further <- dat.finches$further
mu0 <- 0
t.further <- t.test(further,
                    mu = mu0,
                    alternative = "less")
t.stat.further <- t.further$statistic
n <- length(further)
skewness <- skewness(further)
fz <- dnorm(t.stat.further)
Fz <- pnorm(t.stat.further)

edgeworth.approx.error <- (skewness/sqrt(n)) * (((2*(t.stat.further)^2 + 1)/6)*(fz))
# Potential Error in the computation of the p-value is -1.226006e-13.
probability <- Fz + edgeworth.approx.error

### Part B ###
t.values <- seq(-10,10, length = 1000)
fz2 <- dnorm(t.values)

error.vals <- (skewness/sqrt(n)) * (((2*(t.values)^2 + 1)/6)*(fz2))

error.tvals <- tibble(
  t = t.values,
  error = error.vals
)

ggplot()+
  geom_line(data = error.tvals,
            aes(x = t, y = error),
            color = "lightblue")+
  theme_bw()+
  labs(title = "Edgeworth Approximation Error",
       x = "T Values (-10,10)",
       y = "Error")

### Part C ### - error in the rejection region
skewness <- skewness(further)
tval <- qnorm(0.05)
fz3 <- dnorm(tval)

min.n <- ((skewness/(6*(0.10*0.05))) * (2*tval^2 + 1) *fz3)^2

################################################################################
# Question 2: Boostrapping
################################################################################
library(boot)
### Part A ###

# Closer Resamples
R <- 10000
s.closer <- sd(dat.finches$closer)
n.closer <- length(dat.finches$closer)
resamples.closer <- tibble(tstat=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = dat.finches$closer,
                          size = n.closer,
                          replace = T)
  resamples.closer$tstat[i] <- (mean(curr.resample)-0)/(s.closer/sqrt(n.closer))
}
# shift so H0 is true

delta.closer <- mean(resamples.closer$tstat) - 0 # null mu0 = 0

resamples.null.closer <- resamples.closer |>
  mutate(tstat.shifted = tstat - delta.closer)

# CHECK: mean(resamples.null.closer$tstat.shifted)

# Further Resamples
R <- 10000
s.further <- sd(dat.finches$further)
n.further <- length(dat.finches$further)
resamples.further <- tibble(tstat=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = dat.finches$further,
                          size = n.further,
                          replace = T)
  resamples.further$tstat[i] <- (mean(curr.resample)-0)/(s.further/sqrt(n.further))
}
# shift so H0 is true
delta.further <- mean(resamples.further$tstat) - 0 # null mu0 = 0
resamples.null.further <- resamples.further |>
  mutate(tstat.shifted = tstat - delta.further)

# CHECK: mean(resamples.null.further$tstat.shifted)

# Difference Resamples
R <- 10000
s.diff <- sd(dat.finches$diff)
n.diff <- length(dat.finches$diff)
resamples.diff <- tibble(tstat=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = dat.finches$diff,
                          size = n.diff,
                          replace = T)
  resamples.diff$tstat[i] <- (mean(curr.resample)-0)/(s.diff/sqrt(n.diff))
}
# shift so H0 is true
delta.diff <- mean(resamples.diff$tstat) - 0 # null mu0 = 0
resamples.null.diff <- resamples.diff |>
  mutate(tstat.shifted = tstat - delta.diff)

# CHECK: mean(resamples.null.diff$tstat.shifted)

################# Part B ###################

## Closer ##

# Bootstrap P-Value
p.boot.closer <- mean(resamples.null.closer$tstat.shifted >= delta.closer)
# T-Test P-Value
p.t.closer <- t.test(dat.finches$closer,
                     mu = 0,
                     alternative = "greater")
p.t.closer <- p.t.closer$p.value

## Further ##
# Bootstrap P-Value
p.boot.further <- mean(resamples.null.further$tstat.shifted <= delta.further)
# T-Test P-Value
p.t.further <- t.test(dat.finches$further,
                     mu = 0,
                     alternative = "less")
p.t.further <- p.t.further$p.value

## Difference ##
# Bootstrap P-Value
low <- 0 - delta.diff
high <- 0 + delta.diff
p.low <- mean(resamples.null.diff$tstat.shifted <= low)
p.high <- mean(resamples.null.diff$tstat.shifted >= high)
p.boot.diff <- p.low + p.high
# T-Test P-Value
p.t.diff <- t.test(dat.finches$diff,
                    mu = 0,
                    alternative = "two.sided")
p.t.diff <- p.t.diff$p.value

comparison.pvals <- tibble(
  Data = c("Closer", "Further", "Difference"),
  `Bootstrap P-value` = c(p.boot.closer,
                          p.boot.further,
                          p.boot.diff),
  `T-Test P-value` = c(p.t.closer,
                       p.t.further,
                       p.t.diff)
)
view(comparison.pvals)

################# Part C ###################

#Closer
percentile.boot.closer <- quantile(resamples.null.closer$tstat.shifted, 0.05)
percentile.t.closer <- qt(0.05, df = n.closer -1)

#Further
percentile.boot.further <- quantile(resamples.null.further$tstat.shifted, 0.05)
percentile.t.further <- qt(0.05, df = n.further -1)

#Difference
percentile.boot.diff <- quantile(resamples.null.diff$tstat.shifted, 0.05)
percentile.t.diff <- qt(0.05, df = n.diff -1)

comparison.percentile <- tibble(
  Data = c("Closer", "Further", "Difference"),
  `Bootstrap Percentile` = c(percentile.boot.closer,
                             percentile.boot.further,
                             percentile.boot.diff),
  `T-Test Percentile` = c(percentile.t.closer,
                          percentile.t.further,
                          percentile.t.diff)
)
view(comparison.percentile)

################# Part D ###################
# Confidence Interval
### use resamples for boot -- need this for x bar (floruorences)
# Closer
quantile(resamples.null.closer$tstat.shifted, c(0.025, 0.975))

# Further
quantile(resamples.null.further$tstat.shifted, c(0.025, 0.975))

# Diff
quantile(resamples.null.diff$tstat.shifted, c(0.025, 0.975))
