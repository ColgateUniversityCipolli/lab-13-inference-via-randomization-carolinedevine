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
resamples.closer <-  tibble(tstat=rep(NA, R),
                            xbar=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = dat.finches$closer,
                          size = n.closer,
                          replace = T)
  resamples.closer$tstat[i] <- (mean(curr.resample)-0)/(s.closer/sqrt(n.closer))
  resamples.closer$xbar[i] <- mean(curr.resample)
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
resamples.further <- tibble(tstat=rep(NA, R),
                            xbar=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = dat.finches$further,
                          size = n.further,
                          replace = T)
  resamples.further$tstat[i] <- (mean(curr.resample)-0)/(s.further/sqrt(n.further))
  resamples.further$xbar[i] <- mean(curr.resample)
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
resamples.diff <-  tibble(tstat=rep(NA, R),
                          xbar=rep(NA, R))
for(i in 1:R){
  curr.resample <- sample(x = dat.finches$diff,
                          size = n.diff,
                          replace = T)
  resamples.diff$tstat[i] <- (mean(curr.resample)-0)/(s.diff/sqrt(n.diff))
  resamples.diff$xbar[i] <- mean(curr.resample)
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
# Want to use shifted t stat for the p-value
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
# want to use resample x bars (not shifted) for the confidence interval
# Closer
CI.boot.closer <- quantile(resamples.null.closer$xbar, c(0.025, 0.975))
CI.t.closer <- t.test(x=dat.finches$closer, mu = 0, conf.level = 0.95, alternative = "two.sided")
CI.t.closer <- CI.t.closer$conf.int

# Further
CI.boot.further <- quantile(resamples.null.further$xbar, c(0.025, 0.975))
CI.t.further <- t.test(x=dat.finches$further, mu = 0, conf.level = 0.95, alternative = "two.sided")
CI.t.further <- CI.t.further$conf.int

# Diff
CI.boot.diff <- quantile(resamples.null.diff$xbar, c(0.025, 0.975))
CI.t.diff <- t.test(x=dat.finches$diff, mu = 0, conf.level = 0.95, alternative = "two.sided")
CI.t.diff <- CI.t.diff$conf.int


################################################################################
# Question 3: Randomization Test
################################################################################

################################################################################
# Closer (#3)
################################################################################# Since we do not know the generating distribution for closer, furhter, and difference data
# We need to preform randomization procedure.
mu0 <- 0
R <- 10000
rand.closer <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift.closer <- dat.finches$closer - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift.closer *
    sample(x = c(-1, 1),
           size = length(x.shift.closer),
           replace = T)
  
  rand.closer$xbars[i] <- mean(curr.rand)
}
rand.closer <- rand.closer |>
  mutate(xbars = xbars + mu0) # shifting back

# p-value 
obs.mean.closer <- mean(dat.finches$closer)
p.rand.closer <- mean(rand.closer$xbars >= obs.mean.closer)

## Confidence Interval ##
R <- 1000
mu0.iterate <- 0.01
starting.point.closer <- mean(dat.finches$closer)

mu.lower.closer <- starting.point.closer
repeat{
  rand.closer <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.closer <- dat.finches$closer - mu.lower.closer
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.closer *
      sample(x = c(-1, 1),
             size = length(x.shift.closer),
             replace = T)
    
    rand.closer$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.closer <- rand.closer |>
    mutate(xbars = xbars + mu.lower.closer) # shifting back
  
  # p-value  (one-sided)
  obs.mean.closer <- mean(dat.finches$closer)
  p.val.closer <- mean(rand.closer$xbars >= obs.mean.closer)

if(p.val.closer < 0.05){
  break
}else{
  mu.lower.closer <- mu.lower.closer - mu0.iterate
}
}

mu.upper.closer <- starting.point.closer
repeat{
  rand.closer <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.closer <- dat.finches$closer - mu.upper.closer
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.closer *
      sample(x = c(-1, 1),
             size = length(x.shift.closer),
             replace = T)
    
    rand.closer$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.closer <- rand.closer |>
    mutate(xbars = xbars + mu.upper.closer) # shifting back
  
  # p-value  (one-sided)
  obs.mean.closer <- mean(dat.finches$closer)
  p.val.closer <- mean(rand.closer$xbars <= obs.mean.closer)
  
  if(p.val.closer < 0.05){
    break
  }else{
    mu.upper.closer <- mu.upper.closer + mu0.iterate
  }
}

closer.rand.CI <- c(mu.lower.closer, mu.upper.closer)


################################################################################
# Further (#3)
################################################################################
# We need to preform randomization procedure.
mu0 <- 0
R <- 10000
rand.further <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift.further <- dat.finches$further - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift.further *
    sample(x = c(-1, 1),
           size = length(x.shift.further),
           replace = T)
  
  rand.further$xbars[i] <- mean(curr.rand)
}
rand.further <- rand.further |>
  mutate(xbars = xbars + mu0) # shifting back

# p-value 
obs.mean.further <- mean(dat.finches$further)
p.rand.further <- mean(rand.further$xbars <= obs.mean.further)

## Confidence Interval ##
R <- 1000
mu0.iterate <- 0.01
starting.point.further <- mean(dat.finches$further)

mu.lower.further <- starting.point.further
repeat{
  rand.further <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.further <- dat.finches$further - mu.lower.further
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.further *
      sample(x = c(-1, 1),
             size = length(x.shift.further),
             replace = T)
    
    rand.further$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.further <- rand.further |>
    mutate(xbars = xbars + mu.lower.further) # shifting back
  
  # p-value  (one-sided)
  obs.mean.further <- mean(dat.finches$further)
  p.val.further <- mean(rand.further$xbars >= obs.mean.further)
  
  if(p.val.further < 0.05){
    break
  }else{
    mu.lower.further <- mu.lower.further - mu0.iterate
  }
}

mu.upper.further <- starting.point.further
repeat{
  rand.further <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.further <- dat.finches$further - mu.upper.further
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.further *
      sample(x = c(-1, 1),
             size = length(x.shift.further),
             replace = T)
    
    rand.further$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.further <- rand.further |>
    mutate(xbars = xbars + mu.upper.further) # shifting back
  
  # p-value  (one-sided)
  obs.mean.further <- mean(dat.finches$further)
  p.val.further <- mean(rand.further$xbars <= obs.mean.further)
  
  if(p.val.further < 0.05){
    break
  }else{
    mu.upper.further <- mu.upper.further + mu0.iterate
  }
}

further.rand.CI <- c(mu.lower.further, mu.upper.further)



################################################################################
# Diff (#3)
################################################################################

# part a
# We need to preform randomization procedure.
mu0 <- 0
R <- 10000
rand.diff <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift.diff <- dat.finches$diff - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift.diff *
    sample(x = c(-1, 1),
           size = length(x.shift.diff),
           replace = T)
  
  rand.diff$xbars[i] <- mean(curr.rand)
}
rand.diff <- rand.diff |>
  mutate(xbars = xbars + mu0) # shifting back

# part b
# p-value 
(delta.diff <- abs(mean(dat.finches$diff) - mu0))
(low.diff <- mu0 - delta.diff) # mirror
(high.diff <- mu0 + delta.diff)   # xbar

p.rand.diff <- mean(rand.diff$xbars <= low.diff) + 
  mean(rand.diff$xbars >= high.diff)

# part b
## Confidence Interval ##
R <- 1000
mu0.iterate <- 0.01
starting.point.diff <- mean(dat.finches$diff)

mu.lower.diff <- starting.point.diff
repeat{
  rand.diff <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.diff <- dat.finches$diff - mu.lower.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.diff *
      sample(x = c(-1, 1),
             size = length(x.shift.diff),
             replace = T)
    
    rand.diff$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.diff <- rand.diff |>
    mutate(xbars = xbars + mu.lower.diff) # shifting back
  
  # p-value  (one-sided)
  (delta.diff <- abs(mean(dat.finches$diff) - mu.lower.diff))
  (low.diff <- mu.lower.diff - delta.diff) # mirror
  (high.diff <- mu.lower.diff + delta.diff)   # xbar
  (p.val.diff <- mean(rand.diff$xbars <= low.diff) +
      mean(rand.diff$xbars >= high.diff))
  
  if(p.val.diff < 0.05){
    break
  }else{
    mu.lower.diff <- mu.lower.diff - mu0.iterate
  }
}

mu.upper.diff <- starting.point.diff
repeat{
  rand.diff <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift.diff <- dat.finches$diff - mu.upper.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift.diff *
      sample(x = c(-1, 1),
             size = length(x.shift.diff),
             replace = T)
    
    rand.diff$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand.diff <- rand.diff |>
    mutate(xbars = xbars + mu.upper.diff) # shifting back
  
  # p-value 
  (delta.diff <- abs(mean(dat.finches$closer) - mu.upper.diff))
  (low.diff <- mu.upper.diff - delta.diff) # mirror
  (high.diff <- mu.upper.diff + delta.diff)   # xbar
  (p.val.diff <- mean(rand.closer$xbars <= low.diff) +
      mean(rand.diff$xbars >= high.diff))
  
  if(p.val.diff < 0.05){
    break
  }else{
    mu.upper.diff <- mu.upper.diff + mu0.iterate
  }
}

c(mu.lower.diff, mu.upper.diff)


