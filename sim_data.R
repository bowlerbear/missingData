library(spatstat)
library(raster)
library(tidyverse)
library(tmap)

### Params #####


### State ######

# Study Region
bbox <- owin(c(0,1),c(0,1))

# Linear predictor of intensity - affected by y
intercept = 4; beta.x = 0; beta.y = 2
lp <-  function(x,y){ exp(intercept + beta.x * x +  beta.y * y)}

# Point process intensity - points per unit area
lambda <- as.im(lp, bbox)
plot(lambda)

# Simulate point pattern - total number of points is lambda x area of window
pp <- rpoispp(lambda)
plot(pp)

# Grid size
n_grids = 20

# Calc points per grid 
grid_raster <- raster(pixellate(pp, dimyx = c(n_grids, n_grids)))
names(grid_raster) <- "state"
grid_df <- as.data.frame(grid_raster, xy=T)
t0 <- tm_shape(grid_raster) + tm_raster(style = "cont")

#### summarise #######

#modelled relationships 
lm1 <- glm(state ~ x + y, data = grid_df, family="poisson")
coef(lm1)

#get total N
state_N = sum(predict(lm1, type="response"))


### MCAR #############

# Probability of grid cell being sampled - assume only affected by x
intercept = -2; beta.x = 4; beta.y = 0 
grid_df$prob_sampled <- boot::inv.logit(intercept + beta.x * grid_df$x + beta.y * grid_df$y) 
grid_df$sampled <- sapply(grid_df$prob_sampled, FUN = function(x){rbinom(1,1,x)})

# Insert missing data at unsampled cells
grid_df <- grid_df %>% mutate(count_1 = ifelse(sampled==0, NA, state))

# Make new raster with missing values
obs_raster <- rasterFromXYZ(grid_df[,c("x", "y", "count_1")])
t1 <- tm_shape(obs_raster) + tm_raster(style = "cont", title="MCAR Count", colorNA= "grey")

#### summarise #######

#modelled relationships 
lm1 <- glm(count_1 ~ x + y, data = grid_df, family="poisson")
sample1_coefs <- coef(lm1)

#get total N
sample1_N = sum(predict(lm1, type="response", newdata = grid_df))


### MAR #############

# Probability of grid cell being sampled - also affected by y
intercept = -2; beta.x = 0; beta.y = 4 
grid_df$prob_sampled <- boot::inv.logit(intercept + beta.x * grid_df$x + beta.y * grid_df$y) 
grid_df$sampled <- sapply(grid_df$prob_sampled, FUN = function(x){rbinom(1,1,x)})

# Insert missing data at unsampled cells
grid_df <- grid_df %>% mutate(count_2 = ifelse(sampled==0, NA, state))

# Make new raster with missing values
obs_raster <- rasterFromXYZ(grid_df[,c("x", "y", "count_2")])
t2 <- tm_shape(obs_raster) + tm_raster(style = "cont", title="MAR Count", colorNA= "grey")

#### summarise #######

#modelled relationships 
lm1 <- glm(count_2 ~ x + y, data = grid_df, family="poisson")
sample2_coefs <- coef(lm1)

#get total N
sample2_N = sum(predict(lm1, type="response", newdata = grid_df))

### MNAR #############

# Probability of grid cell being sampled - also affected by y 
intercept = -2; beta.x = 0; beta.y = 4 
grid_df$prob_sampled <- boot::inv.logit(intercept + beta.x * grid_df$x + beta.y * grid_df$y) 
grid_df$sampled <- sapply(grid_df$prob_sampled, FUN = function(x){rbinom(1,1,x)})

# Insert missing data at unsampled cells
grid_df <- grid_df %>% mutate(count_3 = ifelse(sampled==0, NA, state))

# Make new raster with missing values
obs_raster <- rasterFromXYZ(grid_df[,c("x", "y", "count_3")])
t3 <- tm_shape(obs_raster) + tm_raster(style = "cont", title="MNAR Count", colorNA= "grey")

#### summarise #######

#modelled relationships 
lm1 <- glm(count_1 ~ x, data = grid_df, family="poisson")
sample2_coefs <- coef(lm1)


#get total N
sample3_N = sum(predict(lm1, type="response", newdata = grid_df))

### Pref sampling #############

# Probability of grid cell being sampled - also affected by y 
intercept = -2; beta.x = 0; beta.y = 4 
grid_df$state_scale = (grid_df$state - min(grid_df$state))/(max(grid_df$state) - min(grid_df$state))
grid_df$prob_sampled <- boot::inv.logit(intercept + beta.y * grid_df$state_scale) 
grid_df$sampled <- sapply(grid_df$prob_sampled, FUN = function(x){rbinom(1,1,x)})

# Insert missing data at unsampled cells
grid_df <- grid_df %>% mutate(count_3 = ifelse(sampled==0, NA, state))

# Make new raster with missing values
obs_raster <- rasterFromXYZ(grid_df[,c("x", "y", "count_3")])
t3 <- tm_shape(obs_raster) + tm_raster(style = "cont", title="MNAR Count", colorNA= "grey")

#### summarise #######

#modelled relationships 
lm1 <- glm(count_1 ~ x, data = grid_df, family="poisson")
sample2_coefs <- coef(lm1)


#get total N
sample3_N = sum(predict(lm1, type="response", newdata = grid_df))

### plot ####

tmap_arrange(t0, t1, t2, t3, nrow=1)

### package up ####


### end #######
