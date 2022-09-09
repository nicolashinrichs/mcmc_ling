# load data and copy
# load csv file w first as names

library(MCMCpack)

tail(stats_clean)
nrow(stats_clean)
d <- stats_clean[1:900,]

# standardize variables
library(rethinking)
d$R <- standardize( d$RESPONSE)
d$S <- standardize( d$SUBJECT)
d$G <- standardize( d$GROUP)
d$P <- standardize( d$PACKAGING)
d$V <- standardize( d$VERB)
d$SS <- standardize( d$Segments)
d$ET <- standardize( d$`EVENT TYPE`)
d$ES <- standardize( d$`EVENT STRUCTURE`)
d$FT <- standardize( d$`FRAME TYPE`)

# regresion de packaging en mcmc

data.mcmcpack.add <- MCMCregress(P ~ R + S + G + ET + ES + FT, data = d)

# check through plot

plot(data.mcmcpack.add)

# compare raw vs model

mcmc = as.matrix(data.mcmcpack.add)
# generate a model matrix
Xmat = model.matrix(~ R + S + G + ET + ES + FT, data = d)
## get median parameter estimates
coefs = mcmc[, 1:7]
fit = coefs %*% t(Xmat)
## replace y with dep (P) and data with d
yRep = sapply(1:nrow(mcmc), function(i) rnorm(nrow(d), fit[i,
], sqrt(mcmc[i, "sigma2"])))
ggplot() + geom_density(data = NULL, aes(x = as.vector(yRep),
               fill = "Model"), alpha = 0.5) + geom_density(data = d,
                                        aes(x = P, fill = "Obs"), alpha = 0.5)

# posteriors

library(bayesplot)
mcmc_intervals(as.matrix(data.mcmcpack.add), regex_pars = "Intercept|R|S|G|ET|ES|FT|sigma")
mcmc_areas(as.matrix(data.mcmcpack.add), regex_pars = "Intercept|R|S|G|ET|ES|FT|sigma")

# posterior parameter summaries

library(broom)
tidyMCMC(data.mcmcpack.add, conf.int = TRUE, conf.method = "HPDinterval")


mcmc = data.mcmcpack.add
## Calculate the fitted values

mean.x1 = mean(d$R)
mean.x2 = mean(d$S)
mean.x3 = mean(d$G)
mean.x4 = mean(d$ET)
mean.x5 = mean(d$ES)
mean.x6 = mean(d$FT)

newdata = rbind(
     data.frame(R = seq(min(d$R, na.rm = TRUE), max(d$R,na.rm = TRUE), len = 100), S = 0, G = 0, ET = 0, ES = 0, FT = 0, Pred = 1),
     data.frame(R = 0, S = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), G = 0, ET = 0, ES = 0, FT = 0, Pred = 2),
     data.frame(R = 0, S = 0, G = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), ET = 0, ES = 0, FT = 0, Pred = 3),
     data.frame(R = 0, S = 0, G = 0, ET = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), ES = 0, FT = 0, Pred = 4),
     data.frame(R = 0, S = 0, G = 0, ET = 0, ES = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), FT = 0, Pred = 5),
     data.frame(R = 0, S = 0, G = 0, ET = 0, ES = 0, FT = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), Pred = 6)
     )
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)
coefs = mcmc[, 1:7]
fit = coefs %*% t(Xmat)
newdata = newdata %>% mutate(x1 = R + mean.x1, x2 = S + mean.x2, x3 = G + mean.x3, x4 = ET + mean.x4, x5 = ES + mean.x5, x6 = FT + mean.x6) %>%
  cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "HPDinterval")) %>%
  mutate(x = dplyr:::recode(Pred, x1, x2, x3, x4, x5, x6))

ggplot(newdata, aes(y = estimate, x = x)) + geom_line() + geom_ribbon(aes(ymin = conf.low,
                                                                          ymax = conf.high), fill = "blue", alpha = 0.3) + scale_y_continuous("Y") +
  scale_x_continuous("X") + theme_classic() + facet_wrap(~Pred)

## Calculate partial residuals fitted values
fdata = rdata = rbind(data.frame(R = d$R, S = 0, G = 0, ET = 0, ES = 0, FT = 0, Pred = 1), 
                      data.frame(R = 0, S = d$S, G = 0, ET = 0, ES = 0, FT = 0, Pred = 2),
                      data.frame(R = 0, S = 0, G = d$G, ET = 0, ES = 0, FT = 0, Pred = 3),
                      data.frame(R = 0, S = 0, G = 0, ET = d$ET, ES = 0, FT = 0, Pred = 4),
                      data.frame(R = 0, S = 0, G = 0, ET = 0, ES = d$ES, FT = 0, Pred = 5),
                      data.frame(R = 0, S = 0, G = 0, ET = 0, ES = 0, FT = d$FT, Pred = 6))
fMat = rMat = model.matrix(~R + S + G + ET + ES + FT, fdata)
fit = as.vector(apply(coefs, 2, median) %*% t(fMat))
resid = as.vector(d$P - apply(coefs, 2, median) %*% t(rMat))
rdata = rdata %>% mutate(partial.resid = resid + fit) %>% mutate(x1 = R +
mean.x1, x2 = S + mean.x2, x3 = G + mean.x3, x4 = ET + mean.x4, x5 = ES + mean.x5, x6 = FT + mean.x6) %>% mutate(x = dplyr:::recode(Pred, x1,
                          x2, x3, x4, x5, x6))

ggplot(newdata, aes(y = estimate, x = x)) + geom_point(data = rdata, aes(y = partial.resid),
                                           color = "gray") + geom_line() + geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                                                                     fill = "blue", alpha = 0.3) + scale_y_continuous("Y") + theme_classic() +
  facet_wrap(~Pred, strip.position = "bottom", labeller = label_bquote("x" *
                                                                         .(Pred))) + theme(axis.title.x = element_blank(), strip.background = element_blank(),
                                                                                           strip.placement = "outside")
# effect sizes

mcmc = data.mcmcpack.add
newdata = expand.grid(R = c(min(d$R), max(d$R)), S = c(min(d$S), max(d$S)), G = c(min(d$G), max(d$G)), ET = c(min(d$ET), max(d$ET)), ES = c(min(d$ES), max(d$ES)), FT = c(min(d$FT), max(d$FT))) 
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)
coefs = mcmc[, c("(Intercept)", "R", "S", "G", "ET", "ES", "FT")]
fit = coefs %*% t(Xmat)
s1 = seq(1, 9, b = 2)
s2 = seq(2, 10, b = 2)

## Raw effect size
(RES = tidyMCMC(as.mcmc(fit[, s2] - fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval"))

## Cohen's D
cohenD = (fit[, s2] - fit[, s1])/sqrt(mcmc[, "sigma2"])
(cohenDES = tidyMCMC(as.mcmc(cohenD), conf.int = TRUE, conf.method = "HPDinterval"))

# Percentage change (relative to Group A)
ESp = 100 * (fit[, s2] - fit[, s1])/fit[, s1]
(PES = tidyMCMC(as.mcmc(ESp), conf.int = TRUE, conf.method = "HPDinterval"))

# Probability that the effect is greater than 50% (an increase of >50%)
(p50 = apply(ESp, 2, function(x) sum(x > 50)/length(x)))

## fractional change
(FES = tidyMCMC(as.mcmc(fit[, s2]/fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval"))


## SDs

mcmc = data.mcmcpack.add
Xmat = model.matrix(~R + S + G + ET + ES + FT, data = d)
sd.x1 = abs(mcmc[, "R"]) * sd(Xmat[, "R"])
sd.x2 = abs(mcmc[, "S"]) * sd(Xmat[, "S"])
sd.x3 = abs(mcmc[, "G"]) * sd(Xmat[, "G"])
sd.x4 = abs(mcmc[, "ET"]) * sd(Xmat[, "ES"])
sd.x5 = abs(mcmc[, "ES"]) * sd(Xmat[, "ES"])
sd.x6 = abs(mcmc[, "FT"]) * sd(Xmat[, "FT"])
sd.x = sd.x1 + sd.x2 + sd.x3 + sd.x4 + sd.x5 + sd.x6

# generate a model matrix
newdata = d
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)

## get median parameter estimates
coefs = mcmc[, c("(Intercept)", "R", "S", "G", "ET", "ES", "FT")]
fit = coefs %*% t(Xmat)
resid = sweep(fit, 2, d$P, "-")
sd.resid = apply(resid, 1, sd)

sd.all = cbind(sd.x1, sd.x2, sd.x3, sd.x4, sd.x5, sd.x6, sd.resid)
(fpsd = tidyMCMC(sd.all, conf.int = TRUE, conf.method = "HPDinterval"))

# OR expressed as a percentage
(fpsd.p = tidyMCMC(100 * sd.all/rowSums(sd.all), estimate.method = "median",
                   conf.int = TRUE, conf.method = "HPDinterval"))

## we can even plot this as a Bayesian ANOVA table
ggplot(fpsd, aes(y = estimate, x = term)) + geom_pointrange(aes(ymin = conf.low,
         ymax = conf.high)) + geom_text(aes(label = sprintf("%.2f%%", fpsd.p$estimate),
          vjust = -1)) + scale_y_continuous("Finite population standard deviation") +
  scale_x_discrete() + coord_flip() + theme_classic()

# define beta of VERB

# regresion de packaging en mcmc 

data.mcmcpack.add <- MCMCregress(SS ~ R + S + G + ET + ES + FT, data = d)

# check through plot

plot(data.mcmcpack.add)

# compare raw vs model

mcmc = as.matrix(data.mcmcpack.add)
# generate a model matrix
Xmat = model.matrix(~ R + S + G + ET + ES + FT, data = d)
## get median parameter estimates
coefs = mcmc[, 1:7]
fit = coefs %*% t(Xmat)
## replace y with dep (P) and data with d
yRep = sapply(1:nrow(mcmc), function(i) rnorm(nrow(d), fit[i,
], sqrt(mcmc[i, "sigma2"])))
ggplot() + geom_density(data = NULL, aes(x = as.vector(yRep),
                                         fill = "Model"), alpha = 0.5) + geom_density(data = d,
                                                                                      aes(x = SS, fill = "Obs"), alpha = 0.5)

# posteriors

library(bayesplot)
mcmc_intervals(as.matrix(data.mcmcpack.add), regex_pars = "Intercept|R|S|G|ET|ES|FT|sigma")
mcmc_areas(as.matrix(data.mcmcpack.add), regex_pars = "Intercept|R|S|G|ET|ES|FT|sigma")

# posterior parameter summaries

library(broom)
tidyMCMC(data.mcmcpack.add, conf.int = TRUE, conf.method = "HPDinterval")


mcmc = data.mcmcpack.add
## Calculate the fitted values

mean.x1 = mean(d$R)
mean.x2 = mean(d$S)
mean.x3 = mean(d$G)
mean.x4 = mean(d$ET)
mean.x5 = mean(d$ES)
mean.x6 = mean(d$FT)

newdata = rbind(
  data.frame(R = seq(min(d$R, na.rm = TRUE), max(d$R,na.rm = TRUE), len = 100), S = 0, G = 0, ET = 0, ES = 0, FT = 0, Pred = 1),
  data.frame(R = 0, S = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), G = 0, ET = 0, ES = 0, FT = 0, Pred = 2),
  data.frame(R = 0, S = 0, G = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), ET = 0, ES = 0, FT = 0, Pred = 3),
  data.frame(R = 0, S = 0, G = 0, ET = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), ES = 0, FT = 0, Pred = 4),
  data.frame(R = 0, S = 0, G = 0, ET = 0, ES = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), FT = 0, Pred = 5),
  data.frame(R = 0, S = 0, G = 0, ET = 0, ES = 0, FT = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), Pred = 6)
)
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)
coefs = mcmc[, 1:7]
fit = coefs %*% t(Xmat)
newdata = newdata %>% mutate(x1 = R + mean.x1, x2 = S + mean.x2, x3 = G + mean.x3, x4 = ET + mean.x4, x5 = ES + mean.x5, x6 = FT + mean.x6) %>%
  cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "HPDinterval")) %>%
  mutate(x = dplyr:::recode(Pred, x1, x2, x3, x4, x5, x6))

ggplot(newdata, aes(y = estimate, x = x)) + geom_line() + geom_ribbon(aes(ymin = conf.low,
                                                                          ymax = conf.high), fill = "blue", alpha = 0.3) + scale_y_continuous("Y") +
  scale_x_continuous("X") + theme_classic() + facet_wrap(~Pred)

## Calculate partial residuals fitted values
fdata = rdata = rbind(data.frame(R = d$R, S = 0, G = 0, ET = 0, ES = 0, FT = 0, Pred = 1), 
                      data.frame(R = 0, S = d$S, G = 0, ET = 0, ES = 0, FT = 0, Pred = 2),
                      data.frame(R = 0, S = 0, G = d$G, ET = 0, ES = 0, FT = 0, Pred = 3),
                      data.frame(R = 0, S = 0, G = 0, ET = d$ET, ES = 0, FT = 0, Pred = 4),
                      data.frame(R = 0, S = 0, G = 0, ET = 0, ES = d$ES, FT = 0, Pred = 5),
                      data.frame(R = 0, S = 0, G = 0, ET = 0, ES = 0, FT = d$FT, Pred = 6))
fMat = rMat = model.matrix(~R + S + G + ET + ES + FT, fdata)
fit = as.vector(apply(coefs, 2, median) %*% t(fMat))
resid = as.vector(d$SS - apply(coefs, 2, median) %*% t(rMat))
rdata = rdata %>% mutate(partial.resid = resid + fit) %>% mutate(x1 = R +
                                                                   mean.x1, x2 = S + mean.x2, x3 = G + mean.x3, x4 = ET + mean.x4, x5 = ES + mean.x5, x6 = FT + mean.x6) %>% mutate(x = dplyr:::recode(Pred, x1,
                                                                                                                                                                                                       x2, x3, x4, x5, x6))

ggplot(newdata, aes(y = estimate, x = x)) + geom_point(data = rdata, aes(y = partial.resid),
                                                       color = "gray") + geom_line() + geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                                                                                                   fill = "blue", alpha = 0.3) + scale_y_continuous("Y") + theme_classic() +
  facet_wrap(~Pred, strip.position = "bottom", labeller = label_bquote("x" *
                                                                         .(Pred))) + theme(axis.title.x = element_blank(), strip.background = element_blank(),
                                                                                           strip.placement = "outside")
# effect sizes

mcmc = data.mcmcpack.add
newdata = expand.grid(R = c(min(d$R), max(d$R)), S = c(min(d$S), max(d$S)), G = c(min(d$G), max(d$G)), ET = c(min(d$ET), max(d$ET)), ES = c(min(d$ES), max(d$ES)), FT = c(min(d$FT), max(d$FT))) 
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)
coefs = mcmc[, c("(Intercept)", "R", "S", "G", "ET", "ES", "FT")]
fit = coefs %*% t(Xmat)
s1 = seq(1, 9, b = 2)
s2 = seq(2, 10, b = 2)

## Raw effect size
(RES = tidyMCMC(as.mcmc(fit[, s2] - fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval"))

## Cohen's D
cohenD = (fit[, s2] - fit[, s1])/sqrt(mcmc[, "sigma2"])
(cohenDES = tidyMCMC(as.mcmc(cohenD), conf.int = TRUE, conf.method = "HPDinterval"))

# Percentage change (relative to Group A)
ESp = 100 * (fit[, s2] - fit[, s1])/fit[, s1]
(PES = tidyMCMC(as.mcmc(ESp), conf.int = TRUE, conf.method = "HPDinterval"))

# Probability that the effect is greater than 50% (an increase of >50%)
(p50 = apply(ESp, 2, function(x) sum(x > 50)/length(x)))

## fractional change
(FES = tidyMCMC(as.mcmc(fit[, s2]/fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval"))


## SDs

mcmc = data.mcmcpack.add
Xmat = model.matrix(~R + S + G + ET + ES + FT, data = d)
sd.x1 = abs(mcmc[, "R"]) * sd(Xmat[, "R"])
sd.x2 = abs(mcmc[, "S"]) * sd(Xmat[, "S"])
sd.x3 = abs(mcmc[, "G"]) * sd(Xmat[, "G"])
sd.x4 = abs(mcmc[, "ET"]) * sd(Xmat[, "ES"])
sd.x5 = abs(mcmc[, "ES"]) * sd(Xmat[, "ES"])
sd.x6 = abs(mcmc[, "FT"]) * sd(Xmat[, "FT"])
sd.x = sd.x1 + sd.x2 + sd.x3 + sd.x4 + sd.x5 + sd.x6

# generate a model matrix
newdata = d
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)

## get median parameter estimates
coefs = mcmc[, c("(Intercept)", "R", "S", "G", "ET", "ES", "FT")]
fit = coefs %*% t(Xmat)
resid = sweep(fit, 2, d$SS, "-")
sd.resid = apply(resid, 1, sd)

sd.all = cbind(sd.x1, sd.x2, sd.x3, sd.x4, sd.x5, sd.x6, sd.resid)
(fpsd = tidyMCMC(sd.all, conf.int = TRUE, conf.method = "HPDinterval"))

# OR expressed as a percentage
(fpsd.p = tidyMCMC(100 * sd.all/rowSums(sd.all), estimate.method = "median",
                   conf.int = TRUE, conf.method = "HPDinterval"))

## we can even plot this as a Bayesian ANOVA table
ggplot(fpsd, aes(y = estimate, x = term)) + geom_pointrange(aes(ymin = conf.low,
                                                                ymax = conf.high)) + geom_text(aes(label = sprintf("%.2f%%", fpsd.p$estimate),
                                                                                                   vjust = -1)) + scale_y_continuous("Finite population standard deviation") +
  scale_x_discrete() + coord_flip() + theme_classic()


# define beta of Segments

# regresion de packaging en mcmc 

data.mcmcpack.add <- MCMCregress(SS ~ R + S + G + ET + ES + FT, data = d)

# check through plot

plot(data.mcmcpack.add)

# compare raw vs model

mcmc = as.matrix(data.mcmcpack.add)
# generate a model matrix
Xmat = model.matrix(~ R + S + G + ET + ES + FT, data = d)
## get median parameter estimates
coefs = mcmc[, 1:7]
fit = coefs %*% t(Xmat)
## replace y with dep (P) and data with d
yRep = sapply(1:nrow(mcmc), function(i) rnorm(nrow(d), fit[i,
], sqrt(mcmc[i, "sigma2"])))
ggplot() + geom_density(data = NULL, aes(x = as.vector(yRep),
                                         fill = "Model"), alpha = 0.5) + geom_density(data = d,
                                                                                      aes(x = SS, fill = "Obs"), alpha = 0.5)

# posteriors

library(bayesplot)
mcmc_intervals(as.matrix(data.mcmcpack.add), regex_pars = "Intercept|R|S|G|ET|ES|FT|sigma")
mcmc_areas(as.matrix(data.mcmcpack.add), regex_pars = "Intercept|R|S|G|ET|ES|FT|sigma")

# posterior parameter summaries

library(broom)
tidyMCMC(data.mcmcpack.add, conf.int = TRUE, conf.method = "HPDinterval")


mcmc = data.mcmcpack.add
## Calculate the fitted values

mean.x1 = mean(d$R)
mean.x2 = mean(d$S)
mean.x3 = mean(d$G)
mean.x4 = mean(d$ET)
mean.x5 = mean(d$ES)
mean.x6 = mean(d$FT)

newdata = rbind(
  data.frame(R = seq(min(d$R, na.rm = TRUE), max(d$R,na.rm = TRUE), len = 100), S = 0, G = 0, ET = 0, ES = 0, FT = 0, Pred = 1),
  data.frame(R = 0, S = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), G = 0, ET = 0, ES = 0, FT = 0, Pred = 2),
  data.frame(R = 0, S = 0, G = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), ET = 0, ES = 0, FT = 0, Pred = 3),
  data.frame(R = 0, S = 0, G = 0, ET = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), ES = 0, FT = 0, Pred = 4),
  data.frame(R = 0, S = 0, G = 0, ET = 0, ES = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), FT = 0, Pred = 5),
  data.frame(R = 0, S = 0, G = 0, ET = 0, ES = 0, FT = seq(min(d$S, na.rm = TRUE), max(d$S, na.rm = TRUE), len = 100), Pred = 6)
)
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)
coefs = mcmc[, 1:7]
fit = coefs %*% t(Xmat)
newdata = newdata %>% mutate(x1 = R + mean.x1, x2 = S + mean.x2, x3 = G + mean.x3, x4 = ET + mean.x4, x5 = ES + mean.x5, x6 = FT + mean.x6) %>%
  cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "HPDinterval")) %>%
  mutate(x = dplyr:::recode(Pred, x1, x2, x3, x4, x5, x6))

ggplot(newdata, aes(y = estimate, x = x)) + geom_line() + geom_ribbon(aes(ymin = conf.low,
                                                                          ymax = conf.high), fill = "blue", alpha = 0.3) + scale_y_continuous("Y") +
  scale_x_continuous("X") + theme_classic() + facet_wrap(~Pred)

## Calculate partial residuals fitted values
fdata = rdata = rbind(data.frame(R = d$R, S = 0, G = 0, ET = 0, ES = 0, FT = 0, Pred = 1), 
                      data.frame(R = 0, S = d$S, G = 0, ET = 0, ES = 0, FT = 0, Pred = 2),
                      data.frame(R = 0, S = 0, G = d$G, ET = 0, ES = 0, FT = 0, Pred = 3),
                      data.frame(R = 0, S = 0, G = 0, ET = d$ET, ES = 0, FT = 0, Pred = 4),
                      data.frame(R = 0, S = 0, G = 0, ET = 0, ES = d$ES, FT = 0, Pred = 5),
                      data.frame(R = 0, S = 0, G = 0, ET = 0, ES = 0, FT = d$FT, Pred = 6))
fMat = rMat = model.matrix(~R + S + G + ET + ES + FT, fdata)
fit = as.vector(apply(coefs, 2, median) %*% t(fMat))
resid = as.vector(d$SS - apply(coefs, 2, median) %*% t(rMat))
rdata = rdata %>% mutate(partial.resid = resid + fit) %>% mutate(x1 = R +
                                                                   mean.x1, x2 = S + mean.x2, x3 = G + mean.x3, x4 = ET + mean.x4, x5 = ES + mean.x5, x6 = FT + mean.x6) %>% mutate(x = dplyr:::recode(Pred, x1,
                                                                                                                                                                                                       x2, x3, x4, x5, x6))

ggplot(newdata, aes(y = estimate, x = x)) + geom_point(data = rdata, aes(y = partial.resid),
                                                       color = "gray") + geom_line() + geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
                                                                                                   fill = "blue", alpha = 0.3) + scale_y_continuous("Y") + theme_classic() +
  facet_wrap(~Pred, strip.position = "bottom", labeller = label_bquote("x" *
                                                                         .(Pred))) + theme(axis.title.x = element_blank(), strip.background = element_blank(),
                                                                                           strip.placement = "outside")
# effect sizes

mcmc = data.mcmcpack.add
newdata = expand.grid(R = c(min(d$R), max(d$R)), S = c(min(d$S), max(d$S)), G = c(min(d$G), max(d$G)), ET = c(min(d$ET), max(d$ET)), ES = c(min(d$ES), max(d$ES)), FT = c(min(d$FT), max(d$FT))) 
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)
coefs = mcmc[, c("(Intercept)", "R", "S", "G", "ET", "ES", "FT")]
fit = coefs %*% t(Xmat)
s1 = seq(1, 9, b = 2)
s2 = seq(2, 10, b = 2)

## Raw effect size
(RES = tidyMCMC(as.mcmc(fit[, s2] - fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval"))

## Cohen's D
cohenD = (fit[, s2] - fit[, s1])/sqrt(mcmc[, "sigma2"])
(cohenDES = tidyMCMC(as.mcmc(cohenD), conf.int = TRUE, conf.method = "HPDinterval"))

# Percentage change (relative to Group A)
ESp = 100 * (fit[, s2] - fit[, s1])/fit[, s1]
(PES = tidyMCMC(as.mcmc(ESp), conf.int = TRUE, conf.method = "HPDinterval"))

# Probability that the effect is greater than 50% (an increase of >50%)
(p50 = apply(ESp, 2, function(x) sum(x > 50)/length(x)))

## fractional change
(FES = tidyMCMC(as.mcmc(fit[, s2]/fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval"))


## SDs

mcmc = data.mcmcpack.add
Xmat = model.matrix(~R + S + G + ET + ES + FT, data = d)
sd.x1 = abs(mcmc[, "R"]) * sd(Xmat[, "R"])
sd.x2 = abs(mcmc[, "S"]) * sd(Xmat[, "S"])
sd.x3 = abs(mcmc[, "G"]) * sd(Xmat[, "G"])
sd.x4 = abs(mcmc[, "ET"]) * sd(Xmat[, "ES"])
sd.x5 = abs(mcmc[, "ES"]) * sd(Xmat[, "ES"])
sd.x6 = abs(mcmc[, "FT"]) * sd(Xmat[, "FT"])
sd.x = sd.x1 + sd.x2 + sd.x3 + sd.x4 + sd.x5 + sd.x6

# generate a model matrix
newdata = d
Xmat = model.matrix(~R + S + G + ET + ES + FT, newdata)

## get median parameter estimates
coefs = mcmc[, c("(Intercept)", "R", "S", "G", "ET", "ES", "FT")]
fit = coefs %*% t(Xmat)
resid = sweep(fit, 2, d$SS, "-")
sd.resid = apply(resid, 1, sd)

sd.all = cbind(sd.x1, sd.x2, sd.x3, sd.x4, sd.x5, sd.x6, sd.resid)
(fpsd = tidyMCMC(sd.all, conf.int = TRUE, conf.method = "HPDinterval"))

# OR expressed as a percentage
(fpsd.p = tidyMCMC(100 * sd.all/rowSums(sd.all), estimate.method = "median",
                   conf.int = TRUE, conf.method = "HPDinterval"))

## we can even plot this as a Bayesian ANOVA table
ggplot(fpsd, aes(y = estimate, x = term)) + geom_pointrange(aes(ymin = conf.low,
                                                                ymax = conf.high)) + geom_text(aes(label = sprintf("%.2f%%", fpsd.p$estimate),
                                                                                                   vjust = -1)) + scale_y_continuous("Finite population standard deviation") +
  scale_x_discrete() + coord_flip() + theme_classic()
