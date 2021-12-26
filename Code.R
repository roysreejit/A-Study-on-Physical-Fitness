rm(list = ls())

setwd("/Users/roysreejit/OneDrive/Academics/ISI/Sem 1/Regression Techniques/Project") # nolint

fitdata <- read.csv("fitness1.csv")
fitdata <- fitdata[, c("oxy", "age", "weight", "runtime",
                       "restpulse", "runpulse", "maxpulse")
]
colnames(fitdata) <- c("Y", "X1", "X2", "X3", "X4", "X5", "X6")
fitdata1 <- fitdata
attach(fitdata)

library(MASS)
library(lattice)
library(car)
library(lmtest)
library(olsrr)
library(randtests)
library(ridge)
library(knitr)
library(leaps)

#2
kable(cor(fitdata1))
pairs(fitdata1)

#3
fm <- lm(Y~ ., fitdata1)
summary(fm)

fm_pred <- predict(fm)
xyplot(fitdata1[, "Y"] + fm_pred ~ 1:31,
       type = c("g", "b"),
       ylab = "Response and Predicted Values",
       xlab = "Index",
       main = list("Full Regression Model"),
       auto.key = list(columns = 2,
                       text = c("Response", "Predicted")
       )
)

#4.1
ols_plot_resid_fit(fm)

#4.1.1
bptest(fm)

#4.2
ols_plot_resid_qq(fm)

#4.2.1
shapiro.test(residuals(fm))

#4.2.2
ols_test_normality(fm)[1]

#4.3
acf(residuals(fm))

#4.3.1
dwtest(fm)

#4.3.2
runs.test(residuals(fm), plot = TRUE)

#4.4.1
vif(fm)

#4.4.2
kappa(fitdata1[, -c(1)])


#5

#5.1
ols_plot_cooksd_bar(fm)

#5.2
ols_plot_dffits(fm)

#5.3
ols_plot_dfbetas(fm)

#5.4
avPlots(fm)

#5.5
ols_plot_resid_stud_fit(fm)
outlierTest(fm)

#5.6
plot(hatvalues(fm), col = 4)
abline(h = 14 / 31, col = 2)
se2 <- seq(1, 31)
idc2 <- (hatvalues(fm) > 14 / 31)
text(se2[idc2], hatvalues(fm)[idc2],
labels = rownames(fitdata1)[idc2], cex = 0.6, pos = 4)

#5.7
plot(fm, which = 5)

#Remedy
fitdata2 <- fitdata[-c(15, 17, 20), ]
fitdata$Y[c(15, 17, 20)] <- predict(lm(Y~ ., fitdata2),
fitdata1[c(15, 17, 20), ])

fm2 <- lm(Y~ ., fitdata)
summary(fm2)

bptest(fm2)
shapiro.test(residuals(fm2))
dwtest(fm2)
vif(fm2)

#6.1
rm <- linearRidge(fitdata, alpha = 0)
summary(rm)
rm_pred <- predict(rm)
xyplot(fitdata1[, "Y"] + rm_pred ~ 1:31,
       type = c("g", "b"),
       ylab = "Response and Predicted Values",
       xlab = "Index",
       main = list("Ridge Regression Model"),
       auto.key = list(columns = 2,
                       text = c("Response", "Predicted")
       )
)

#6.2.1
m1 <- lm(Y ~ X1 + X2 + X3 + X4 + X6, fitdata[, -c(6)])
summary(m1)

ols_plot_resid_fit(m1)
bptest(m1)

ols_plot_resid_qq(m1)
shapiro.test(residuals(m1))
ols_test_normality(m1)[1]

acf(residuals(m1))
dwtest(m1)
runs.test(residuals(m1), plot = TRUE)

vif(m1)
kappa(fitdata[, -c(1, 6)])

#6.2.2
m2 <- lm(Y ~ X1 + X2 + X3 + X4 + X5, fitdata[, -c(7)])
summary(m2)

ols_plot_resid_fit(m2)
bptest(m2)

ols_plot_resid_qq(m2)
shapiro.test(residuals(m2))
ols_test_normality(m2)[1]

acf(residuals(m2))
dwtest(m2)
runs.test(residuals(m2), plot = TRUE)

vif(m2)
kappa(fitdata[, -c(1, 7)])

#7.1
avPlots(m1)

regfit1 <- regsubsets(Y ~ X1 + X2 + X3 + X4 + X6, data = fitdata)
regsumm1 <- summary(regfit1)
regsumm1

mat1 <- cbind(regsumm1$cp, regsumm1$adjr2, regsumm1$bic,
ols_step_all_possible(m1)$aic[c(1, 6, 16, 26, 31)])

rownames(mat1) <- c("X3", "X1, X3", "X1, X2, X3",
"X1, X2, X3, X6", "X1, X2, X3, X4, X6")
colnames(mat1) <- c("cp", "adjr2", "bic", "aic")
kable(mat1)
#stepwise
ols_step_both_aic(m1, details = TRUE)

#7.2
avPlots(m2)

regfit2 <- regsubsets(Y ~ X1 + X2 + X3 + X4 + X5, data = fitdata)
regsumm2 <- summary(regfit2)
regsumm2

mat2 <- cbind(regsumm2$cp, regsumm2$adjr2, regsumm2$bic,
ols_step_all_possible(m2)$aic[c(1, 6, 16, 26, 31)])

rownames(mat2) <- c("X3", "X1, X3", "X1, X3, X5",
"X1, X2, X3, X5", "X1, X2, X3, X4, X5")
colnames(mat2) <- c("cp", "adjr2", "bic", "aic")
kable(mat2)
#stepwise
ols_step_both_aic(m2, details = TRUE)

#7.3
ma <- lm(Y ~ X1 + X2 + X3 + X6, fitdata)
mb <- lm(Y ~ X1 + X2 + X3 + X5, fitdata)

matf <- matrix(c(regsumm1$adjr2[4], AIC(ma), BIC(ma), regsumm1$cp[4],
regsumm2$adjr2[4], AIC(mb), BIC(mb), regsumm2$cp[4]), byrow = TRUE, nrow = 2)
rownames(matf) <- c("Model A", "Model B")
colnames(matf) <- c("adjr2", "aic", "bic", "cp")
kable(matf)

#8
summary(mb)
ols_plot_cooksd_bar(mb)
ols_plot_dffits(mb)
plot(mb, which = 5)

fitdata3 <- fitdata[-c(10), ]
m <- lm(Y ~ X1 + X2 + X3 + X5, fitdata3)
summary(m)

bptest(m)
shapiro.test(residuals(m))
ols_test_normality(m)[1]
dwtest(m)
vif(m)

m_pred <- predict(m)
fitdata4 <- fitdata1[-c(10), ]
xyplot(fitdata4[, "Y"] + m_pred ~ 1:30,
       type = c("g", "b"),
       ylab = "Response and Predicted Values",
       xlab = "Index",
       main = list("Final Regression Model"),
       auto.key = list(columns = 2,
                       text = c("Response", "Predicted")
       )
)
