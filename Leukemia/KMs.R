library(spBayesSurv)
library(ggplot2)
library(survminer)
library(survival)

data("LeukSurv")

################################################################################
# Comparison of Kaplan-Meier curves by district
################################################################################

fit <- survfit(Surv(time, cens) ~ district,
               data = LeukSurv)



ggsurvplot(fit, data = LeukSurv, risk.table = FALSE, linetype = 1, palette = "grey")

library(RColorBrewer)

colorscheme <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'RdYlGn'))(24), 1)

mort <- vector()

for(i in 1:24) mort[i] <- mean(LeukSurv$cens[which(LeukSurv$district == i)])



ggsurvplot(fit, data = LeukSurv, risk.table = FALSE, linetype = 1, palette = colorscheme[order(mort)])

################################################################################
# G-rho test to compare the survival curves by district
################################################################################

surv_diff <- survdiff(Surv(time, cens) ~ district, data = LeukSurv)
surv_diff
