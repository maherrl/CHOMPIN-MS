###########################################################
# This script is the relative abundace analysis for the ##
# CHOMPIN project
# Here I use the abundance of the top five most abundant #
# OTUs and run models
# Created June 21, 2018
# Edited 06/25/18
#Created by Becca Maher
##########################################################

# Libraries------------------------
library(lme4)

# Data-----------------------
# Data frame with meta data and the relative abundances of the top
# five most abundant OTUs by sample
endo <- read.csv(file = "~/Box Sync/CHOMPIN/Manuscript_mal16/R-info/endomodels.csv")
summary(endo)

# Make relative abundances
endo$Endo_rel <- endo$Endozoicimonaceae/endo$sequencing_depth
endo$Desulf_rel <- endo$Desulfovibrionaceae/endo$sequencing_depth
endo$Entero_rel <- endo$Enterobacteriaceae/endo$sequencing_depth
endo$Amoeb_rel <- endo$Amoebophilaceae/endo$sequencing_depth
endo$Morax_rel <- endo$Moraxellaceae/endo$sequencing_depth


# Starting with Endozoicomonaceae, then repeat for each of the five:
# Desulfovibrionaceae, Enterobacteraceae, Ameobaphilaceae, Moraxellaceae
# Check for normality in the data
hist(endo$Endo_rel) # not normally distribute
hist(logitTransform(endo$Endo_rel)) # better
shapiro.test(logitTransform(endo$Endo_rel)) #but sig different from null (normal)

# Add transformed column to data
endo$asin.endo <- asinTransform(endo$Endo_rel)
endo$log.desulf <- log(endo$Desulfovibrionaceae)

# Checking for heteroscedasticity, differences in group spread
dotchart(endo$log.desulf, groups = endo$nutrient)
dotchart(endo$log.desulf, groups = endo$temp)
dotchart(endo$log.desulf, groups = endo$corallivory)
# really no differences

# Mixed effects model selection for inclusion of tank and colony as random effects
# models
full.model <- lmer(asin.endo ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = endo, REML = T)
tank.model <- lmer(asin.endo ~ nutrient*temp*corallivory + (1|tank), data = endo, REML = T)
colony.model <- lmer(asin.endo ~ nutrient*temp*corallivory + (1|colony), data = endo, REML = T)

# Shows the variance in the response variable across random effects
dotplot(ranef(full.model, condVar = TRUE))

# need to check assumptions for random effects:
# 1) residuals are independently distributed
# 2) residuals are from a normal distrubution
qqnorm(ranef(full.model)$colony[,1])
qqline(ranef(full.model)$colony[,1])
qqnorm(ranef(full.model)$tank[,1])
qqline(ranef(full.model)$tank[,1])
# not sure about this

# Determine the random effect structure
# determining if the random effects of tank and colony are significant in the model
anova(tank.model, full.model) # random effect of colony not significant (p=1)
anova(colony.model, full.model) # random effect of tank not significant (p=1)

AIC(tank.model, colony.model, full.model) # full.model highest AIC

# Can remove both random effects for all three measures of alpha diversity
# So now I need to compare a model with both removed to the full model

full.model <- lmer(Endozoicimonaceae ~ temp*nutrient*corallivory + (1|tank) + (1|colony), data = endo, REML = F)
fixed.model <- glm(Endozoicimonaceae ~ temp*nutrient*corallivory, family = quasipoisson, data = endo)

anova(full.model, fixed.model) # p=1
AIC(full.model, fixed.model) # fixed is less
summary(fixed.model)

coeftest(fm_pois, vcov. = sandwich)
summary(fm_pois <- glmer(Endozoicimonaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory + (1|tank) + (1|colony), data = endo, family = "poisson"))
summary(fm_qpois <- glm(Endozoicimonaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory, data = endo, family = "quasipoisson"))
summary(fm_nbin <- glm.nb(Endozoicimonaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory, data = endo))
summary(fm_hurdle0 <- hurdle(Endozoicimonaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory, data = endo, dist = "negbin"))
summary(fm_zinb0 <- zeroinfl(Endozoicimonaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory, data = endo, dist = "negbin"))
fm <- list("ML-Pois" = fm_pois, "Quasi-Pois" = fm_qpois, "NB" = fm_nbin)
sapply(fm, function(x) coef(x)[1:8])

# Test for inclusion of random effects
full.model <- glmer(Moraxellaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory + (1|tank) + (1|colony), data = endo, family = "poisson")
tank.model <- glmer(Moraxellaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory + (1|tank), data = endo, family = "poisson")
colony.model <- glmer(Moraxellaceae ~ offset(log(sequencing_depth)) + temp*nutrient*corallivory + (1|colony), data = endo, family = "poisson")

anova(tank.model, full.model)
anova(colony.model, full.model)


with(endo, tapply(Endozoicimonaceae, corallivory, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
          }))

# Is a mixed model right for your needs?
str(endo)
head(endo)

library(car)
library(MASS)

qqp(endo$Endozoicimonaceae, "norm")
# lognormal
qqp(endo$Endozoicimonaceae, "lnorm")
# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- fitdistr(endo$Desulfovibrionaceae, "Negative Binomial")
qqp(endo$Desulfovibrionaceae, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

poisson <- fitdistr(endo$Desulfovibrionaceae, "Poisson")
qqp(endo$Desulfovibrionaceae, "pois", poisson$estimate)

gamma <- fitdistr(endo$Desulfovibrionaceae, "gamma")
qqp(endo$Endozoicimonaceae, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
