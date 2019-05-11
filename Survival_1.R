## Sorce : http://www.sthda.com/english/wiki/survival-analysis-basics#log-rank-test-comparing-survival-curves-survdiff

install.packages(c("survival", "survminer"))
library("survival")
library("survminer")

## Cox modle function : coxph(formula, data, method)
data("lung")
head(lung)

########### kaplan-Meier survival estimate, survfit() ########
fit <- survfit(Surv(Month..5.Year., vital_status..5.year.) ~ GOLGA4_, data = CRC)
print(fit)

## We want to compute the survival probability by sex
# Summary of survival curves
summary(fit)

# Access to the sort summary table
summary(fit)$table

d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)

# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

ggsurvplot(fit,             # survfit object with calculated statistics.
  pval = TRUE,              # show p-value of log-rank test.
  conf.int = TRUE,          # show confidence intervals for point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",    # customize X axis label.
  break.time.by = 200,      # break X axis in time intervals by 200.
  ggtheme = theme_light(),  # customize plot and risk table with a theme.
  risk.table = "abs_pct",   # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c("Male", "Female"), # change legend labels.
  palette = c("#E7B800", "#2E9FDF")  # custom color palettes.
)

## The survival curves can be shorten using the argument xlim as follow:
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           xlim = c(0, 600))

##  cumulative events
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           fun = "event")

## cummulative hazard
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           fun = "cumhaz")


## Kaplan-Meier life table: summary of survival curves
summary(fit)
res.sum <- surv_summary(fit)
head(res.sum)
attr(res.sum, "table")

######### Log-Rank test comparing survival curves: survdiff()######
surv_diff <- survdiff(Surv(time, status) ~ sex, data = lung)
surv_diff

## Fit complex survival curves
require("survival")
fit2 <- survfit( Surv(time, status) ~ sex + rx + adhere,
                 data = colon )
fit2

# Plot survival curves by sex and facet by rx and adhere
ggsurv <- ggsurvplot(fit2, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw())

ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(rx ~ adhere)

## Sorce: http://www.sthda.com/english/wiki/cox-proportional-hazards-model
## We'll fit the Cox regression using the following covariates: age, sex, ph.ecog and wt.loss.
## We start by computing univariate Cox analyses for all these variables; 
## then we'll fit multivariate cox analyses using two variables to describe how the factors jointly impact on survival.
## Univariate Cox regression
res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
res.cox

covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

## Multivariate Cox regression analysis
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.cox)

## Visualizing the estimated distribution of survival times
# Plot the baseline survival function
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())

# Create the new data  
sex_df <- with(lung,
               data.frame(sex = c(1, 2), 
                          age = rep(mean(age, na.rm = TRUE), 2),
                          ph.ecog = c(1, 1)
               )
)
sex_df

# Survival curves
fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
           ggtheme = theme_minimal())

