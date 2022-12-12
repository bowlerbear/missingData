
# https://github.com/bowlerbear/urbanBias/tree/main/simulations

#### weighing #######

require(ipw)
require(survey)

require(WeightIt)

#https://www.andrewheiss.com/blog/2020/12/01/ipw-binary-continuous/#ipw-manually-binary-treatment


temp <- ipwpoint(exposure = sampled, 
                 family = "binomial", link = "logit",
                 numerator = ~ 1, denominator = ~ y, data = grid_df)
df$Weights2 <- temp$ipw.weights

weights_weightit <- weightit(net ~ income + health,  # Model net use with confounders
                             data = net_data, 
                             estimand = "ATE",  # Find the ATE
                             method = "ps") 

sglm1 <- svyglm(z ~ 1,family=binomial,design = svydesign(~ 1,  weights = ~ weighta, data = df))

### imputation #####