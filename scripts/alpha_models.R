## Mixed effects models for alpha diversity

#load alphadiv dataset

# first looking at distribution of response variables
# Richness
hist(alphadiv$richness)
qqnorm(alphadiv$richness)
abline(0,1)
shapiro.test(alphadiv$richness) # null-hypothesis that the population is normally distributed
# failed the test, log-transformation?
alphadiv$logrich <- log(alphadiv$richness)
hist(alphadiv$logrich)
qqnorm(alphadiv$logrich)
abline(0,1) # this line still doesn't look right
shapiro.test(alphadiv$logrich) # passes the test!
#Evenness
# Tried log transformation (made it worse) and sqrt made it better but still not normal
# Arcsin of the square root because it is proportion data!
hist(alphadiv$evenness)
shapiro.test(alphadiv$evenness)
alphadiv$logeven <- log(alphadiv$evenness)
alphadiv$asineven <- asin(sqrt(alphadiv$evenness))
hist(alphadiv$asineven)
shapiro.test(alphadiv$asineven)
# Faith's PD
hist(alphadiv$faithPD)
qqnorm(alphadiv$faithPD)
alphadiv$logfaithPD <- log(alphadiv$faithPD)
shapiro.test(alphadiv$logfaithPD)
hist(alphadiv$logfaithPD)

# richness
# models
full.model <- lmer(logrich ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = alphadiv, REML = T)
tank.model <- lmer(logrich ~ nutrient*temp*corallivory + (1|tank), data = alphadiv, REML = T)
colony.model <- lmer(logrich ~ nutrient*temp*corallivory + (1|colony), data = alphadiv, REML = T)

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
anova(tank.model, full.model) # random effect of colony not significant (p=0.4)
anova(colony.model, full.model) # random effect of tank not significant (p=1)

# Evenness
# models
full.model.e <- lmer(asineven ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = alphadiv, REML = T)
tank.model.e <- lmer(asineven ~ nutrient*temp*corallivory + (1|tank), data = alphadiv, REML = T)
colony.model.e <- lmer(asineven ~ nutrient*temp*corallivory + (1|colony), data = alphadiv, REML = T)

# Shows the variance in the response variable across random effects
dotplot(ranef(full.model.e, condVar = TRUE))

# need to check assumptions for random effects:
# 1) residuals are independently distributed
# 2) residuals are from a normal distrubution
qqnorm(ranef(full.model.e)$colony[,1])
qqline(ranef(full.model.e)$colony[,1])
qqnorm(ranef(full.model.e)$tank[,1])
qqline(ranef(full.model.e)$tank[,1])
# not sure about this, what does a flat line mean?

# Determine the random effect structure
# determining if the random effects of tank and colony are significant in the model
anova(tank.model.e, full.model.e) # random effect of colony not significant (p=1)
anova(colony.model.e, full.model.e) # random effect of tank not significant (p=1)

# Faith's PD
# models
full.model.f <- lmer(logfaithPD ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = alphadiv, REML = T)
tank.model.f <- lmer(logfaithPD ~ nutrient*temp*corallivory + (1|tank), data = alphadiv, REML = T)
colony.model.f <- lmer(logfaithPD ~ nutrient*temp*corallivory + (1|colony), data = alphadiv, REML = T)

# Shows the variance in the response variable across random effects
dotplot(ranef(full.model.f, condVar = TRUE))

# need to check assumptions for random effects:
# 1) residuals are independently distributed
# 2) residuals are from a normal distrubution
qqnorm(ranef(full.model.f)$colony[,1])
qqline(ranef(full.model.f)$colony[,1])
qqnorm(ranef(full.model.f)$tank[,1])
qqline(ranef(full.model.f)$tank[,1])
# not sure about this, all flat lines

# Determine the random effect structure
# determining if the random effects of tank and colony are significant in the model
anova(tank.model.f, full.model.f) # random effect of colony not significant (p=1)
anova(colony.model.f, full.model.f) # random effect of tank not significant (p=1)

AIC(tank.model, colony.model, full.model) # full.model highest AIC

# Can remove both random effects for all three measures of alpha diversity
# So now I need to compare a model with both removed to the full model

full.model <- lmer(logrich ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = alphadiv, REML = F)
fixed.model <- lm(logrich ~ nutrient*temp*corallivory, data = alphadiv)

anova(full.model, fixed.model) # p=0.7
AIC(full.model, fixed.model) # fixed is less

full.model.ef <- lmer(asineven ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = alphadiv, REML = F)
fixed.model.ef <- lm(asineven ~ nutrient*temp*corallivory, data = alphadiv)

anova(full.model.ef, fixed.model.ef) # p=1
AIC(full.model.ef, fixed.model.ef) # fixed is less

full.model.ff <- lmer(logfaithPD ~ nutrient*temp*corallivory + (1|tank) + (1|colony), data = alphadiv, REML = F)
fixed.model.ff <- lm(logfaithPD ~ nutrient*temp*corallivory, data = alphadiv)

anova(full.model.ff, fixed.model.ff) # p=1
AIC(full.model.ff, fixed.model.ff) # fixed is less

summary(fixed.model)
summary(fixed.model.ef)
summary(fixed.model.ff)

lm.full.maind.r <- lm(logfaithPD ~ nutrient +temp +corallivory, data = alphadiv)
lm.main.r.n <- lm(logfaithPD ~ nutrient, data = alphadiv)
lm.main.r.nt <- lm(logfaithPD ~ nutrient +temp, data = alphadiv)
lm.main.r.nc <- lm(logfaithPD ~ nutrient +corallivory, data = alphadiv)
lm.main.r.tc <- lm(logfaithPD ~ temp + corallivory, data = alphadiv)
anova(lm.full.maind.r, lm.main.r.tc)
anova(lm.main.r.n,lm.main.r.nc)

summary(fixed.model.ef)
anova(fixed.model.ef)
