library(tidyverse)

#generate an x and z - our environmental covariates

x <- runif(50, 0, 1)
z <- runif(50, 0, 1)

#generate species abundance - affected by x

linpred <- 2 + 2 * x
count <- sapply(linpred, function(x){ rpois(1,exp(x))})
qplot(x, count) + scale_y_log10()

#make full dataset
df_full <- data.frame(x, z, count)

#MCAR - sampling affected by z

sampled <- ifelse(z < 0.5, 0, 1)
df_mcar <- df_full[sampled==1,]

#MAR/MNAR - sampling affected by x

sampled <- ifelse(x > 0.5, 0, 1)
df_mar <- df_full[sampled==1,]

#MAR/MNAR - intermediate - reduced sampled by x

linpred <- 4 + -8 * x
prob <- boot::inv.logit(linpred)
#qplot(x,prob)
sampled <- sapply(linpred, function(x){ rbinom(1,1,boot::inv.logit(x))})
df_int <- df_full[sampled==1,]


#MAR/MNAR - non-linear relationship

linpred <- 2 + (2 * x) + (-1 * x^2)
count <- sapply(linpred, function(x){ rpois(1,exp(x))})
qplot(x, count) + scale_y_log10()
df_full_nonlinear <- data.frame(x, z, count)
sampled <- ifelse(x > 0.5, 0, 1)
df_mar_nonlinear <- df_full_nonlinear[sampled==1,]

# plot 1

df_full <- df_full %>% add_column(Missing = "No")
df_mcar <- df_mcar %>% add_column(Missing = "Yes")

all_df <- bind_rows(df_full, df_mcar)

scatter_mcar <- ggplot(data = all_df,
                        aes(x = x, y = count, colour = Missing, fill = Missing)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_brewer("Missing data", type = "qual") +
  scale_fill_brewer("Missing data", type = "qual") +
  theme_classic() +
  ylab("ln Abundance") + xlab("Covariate")

dens_mcar <- ggplot(all_df, aes(x = count, fill = Missing))+
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(type = "qual") +
  ylab("Abundance") + xlab("Density")


#plot 2

df_mar <- df_mar %>% add_column(Missing = "Yes")
all_df <- bind_rows(df_full, df_mar)

scatter_mar <- ggplot(data = all_df,
                       aes(x = x, y = count, colour = Missing, fill = Missing)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_brewer("Missing data", type = "qual") +
  scale_fill_brewer("Missing data", type = "qual") +
  theme_classic() +
  ylab("ln Abundance") + xlab("Covariate")

dens_mar <- ggplot(all_df, aes(x = count, fill = Missing))+
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(type = "qual") +
  ylab("Abundance") + xlab("Density")

#plot 3

df_int <- df_int %>% add_column(Missing = "Yes")
all_df <- bind_rows(df_full, df_int)

scatter_int <- ggplot(data = all_df,
                      aes(x = x, y = count, colour = Missing, fill = Missing)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_y_log10() +
  scale_color_brewer("Missing data", type = "qual") +
  scale_fill_brewer("Missing data", type = "qual") +
  theme_classic() +
  ylab("ln Abundance") + xlab("Covariate")

dens_int <- ggplot(all_df, aes(x = count, fill = Missing))+
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(type = "qual") +
  ylab("Abundance") + xlab("Density")


#plot 4

df_full_nonlinear <- df_full_nonlinear %>% add_column(Missing = "No")
df_mar_nonlinear <- df_mar_nonlinear %>% add_column(Missing = "Yes")
all_df <- bind_rows(df_full_nonlinear, df_mar_nonlinear)

scatter_nonlinear <- ggplot(data = all_df,
                      aes(x = x, y = count, colour = Missing, fill = Missing)) +
  geom_point() +
  stat_smooth(method = "gam") +
  scale_y_log10() +
  scale_color_brewer("Missing data", type = "qual") +
  scale_fill_brewer("Missing data", type = "qual") +
  theme_classic() +
  ylab("ln Abundance") + xlab("Covariate")

dens_nonlinear <- ggplot(all_df, aes(x = count, fill = Missing))+
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(type = "qual") +
  ylab("Abundance") + xlab("Density")


#package all up

library(cowplot)

plot_grid(scatter_mcar,dens_mcar,
          scatter_mar, dens_mar,
          scatter_int, dens_int,
          scatter_nonlinear, dens_nonlinear,
          ncol = 2,
          labels=c("a - MCAR (covariate)","b - MCAR (sample)",
                   "c - MAR (covariate)","d - MAR (sample)",
                   "e - MAR intermediate (covariate)","f - MAR intermediate (sample)",
                   "g - MAR nonlinear (covariate)", "h - MAR nonlinear (sample)"),
          label_x = .1, hjust = 0, vjust = 0.95, label_size = 10)
          
          
          
