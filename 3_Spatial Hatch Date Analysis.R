#R version 4.5.1, code was last run on X/X/2026
#Final code for Climate extremes reshape reproductive phenology: Regional variability of hatch timing in the subarctic North Pacific

#LMM for the spatial analysis of hatch dates for 2018 to 2023

# Paired with:

#Data from all samples, including predicted ages - created in "1_Age-Length Key.R" file
# "2_All_fish_ages_applied_agelengthkey.csv"

# Marginal means from "2_Temporal Hatch Date Analysis.R" file - This file is created in this R script
# "3_Eastern Kodiak Marginal Means 07-23.csv"




#------------------Load Libraries-----------------------------------------------

library(rstudioapi)
library(car)
library(lme4)
library(performance)
library(MuMIn)
library(emmeans)
library(lubridate)
library(patchwork)
library(tidyverse)


#-----------------Set Working Directory-----------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




#-----------------Import csv files----------------------------------------------

all.fish.18.23 <- read.csv("2_All_fish_ages_applied_agelengthkey.csv")

# This is created in "2_Temporal Hatch Date Analysis.R"
eastern.kodiak.07.23 <- read.csv("3_Eastern Kodiak Marginal Means 07-23.csv")




#-----------------Formatting Data-----------------------------------------------

# Define heatwave categories
hw.cats <- tibble(
  year = 2018:2023,
  HW   = c("Adjacent", "During", "Adjacent", "After", "After", "After")
)


# Apply to all years & combine
all.fish.18.23 <- all.fish.18.23 %>%
  left_join(hw.cats, by = "year") %>%
  mutate(
    hatch_DATE       = as.Date(hdate, origin = paste0(year - 1, "-12-31")),
    hatch_DATE_noYear = format(hatch_DATE, "%b-%d"),
    year   = factor(year),
    region = factor(region),
    HW     = factor(HW)
  )


#Creating a "fake" date for plots
all.fish.18.23$hatch_dummyDate <- as.Date(paste0("2020-", all.fish.18.23$hatch_DATE_noYear), format = "%Y-%b-%d")



# Ordering bays/regions and heatwave categories
custom.order <- c("Anton Larsen Bay", "Cook Bay", "Kaiugnak Bay", "Agripina Bay", "Falmouth Bay", "Sand Point", "Baralof Bay")

all.fish.18.23$bay <- fct_relevel(all.fish.18.23$bay, custom.order)

custom.order2 <- c("During", "Adjacent", "After")

all.fish.18.23$HW <- fct_relevel(all.fish.18.23$HW , custom.order2)

custom.order3 <- c("Eastern Kodiak", "Kaiugnak Bay", "Agripina Bay","Shumagin Islands")

all.fish.18.23$region <- fct_relevel(all.fish.18.23$region , custom.order3)




#--------------------Defining Colors for Plots----------------------------------

bay.colors <- c("Anton Larsen Bay" = "blue", "Cook Bay" = "blue", "Kaiugnak Bay" = "forestgreen", "Agripina Bay" = "gold2", 
                 "Falmouth Bay" = "turquoise1", "Sand Point" = "turquoise1", "Baralof Bay" = "turquoise1")

HW.colors <- c("Before" = "darkblue", "During" = "darkred", "Adjacent" = "lightpink", "After" = "steelblue1")

region.colors <- c("Eastern Kodiak" = "blue", "Kaiugnak Bay" = "forestgreen", "Agripina Bay" = "gold2", 
                   "Shumagin Islands" = "turquoise1")




#------------------Plots--------------------------------------------------------

# Hatch date density plot by regions
HD.density <- ggplot(all.fish.18.23, aes(as.Date(hatch_dummyDate, format = "%b-%d"), fill = region)) + 
  geom_density(alpha = 0.5, position = "identity") +
  xlab("Hatch Date") + 
  ylab("Density") +
  labs(fill = NULL) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b-%d", limits = as.Date(c("2020-01-15", "2020-06-15"))) +
  scale_fill_manual(values = region.colors) +
  geom_vline(xintercept = as.Date("2020-04-09"),
             color = "black",
             linewidth = 0.75,
             linetype = "dashed") +
  theme_classic() +
  theme(axis.title = element_text(size = 36, face = "bold"),
        axis.text = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "top",
        strip.text = element_text(size = 36)) +
  guides(fill = guide_legend(nrow = 1))


#Hatch date box plots by bay
HD.bays <- ggplot(all.fish.18.23, aes(bay, as.Date(hatch_dummyDate, format = "%b-%d"), fill = bay)) + 
  geom_boxplot(size = 0.5, alpha = 0.6, color = "gray") +
  ylab("Hatch Date") +
  labs(fill = NULL,
       x = "Bay") +
  scale_y_date(date_breaks = "1 month", date_labels = "%b", limits = as.Date(c("2020-01-15", "2020-06-15"))) +
  theme_classic() +
  guides(fill = "none") +
  scale_fill_manual(values = bay.colors) +
  theme(axis.title = element_text(size = 36, face = "bold"),
        axis.text = element_text(size = 30),
        axis.text.x = element_text(angle = 0),
        legend.text =  element_text(size = 30),
        legend.position = "none",
        strip.text = element_text(size = 36)) 


# Combined plots

combined.hdate <- (HD.density / HD.bays) +
  plot_layout(ncol = 1, nrow = 2) & 
  theme(legend.position = "top")

combined.hdate




#--------------Spatial Models---------------------------------------------------


#Determining best random effects - bay is rank deficient (not all bays sampled in every MHW category) so using region 

HD.model1 <- lmer(hdate ~ region * HW + (1|bay) + (1|year), data = all.fish.18.23, REML = FALSE) 
HD.model2 <- lmer(hdate ~ region * HW + (1|year), data = all.fish.18.23, REML = FALSE)
HD.model3 <- lmer(hdate ~ region * HW + (1|bay), data = all.fish.18.23, REML = FALSE)
HD.model4 <- lm(hdate ~ region * HW, data = all.fish.18.23)

AIC(HD.model1, HD.model2, HD.model3, HD.model4)

#Model 1 has the best random effects

# R squared values and weights of models 1 to 4
r.squaredGLMM(HD.model1)
r.squaredGLMM(HD.model2)
r.squaredGLMM(HD.model3)
summary(HD.model4)

Weights(AIC(HD.model1, HD.model2, HD.model3, HD.model4))


#Determining best fixed effects
HD.model5 <- lmer(hdate ~ region + HW + (1|bay) + (1|year), data = all.fish.18.23, REML = FALSE)
HD.model6 <- lmer(hdate ~ region + (1|bay) + (1|year), data = all.fish.18.23, REML = FALSE)
HD.model7 <- lmer(hdate ~ HW + (1|bay) + (1|year), data = all.fish.18.23, REML = FALSE)

AIC(HD.model1, HD.model5, HD.model6, HD.model7)

#Model 1 has the best fixed effects

#R squared of models 5 to 6
r.squaredGLMM(HD.model5)
r.squaredGLMM(HD.model6)
r.squaredGLMM(HD.model7)

#Weights of all the models

Weights(AIC(HD.model1, HD.model2, HD.model3, HD.model4, HD.model5, HD.model6, HD.model7))


#Function to determine K for each model
calc_K <- function(model) {
  if(inherits(model, "lm")) {
    return(length(coef(model)) + 1)
  } else if(inherits(model, "lmerMod")) {
    fixed <- length(fixef(model))
    random <- length(VarCorr(model))  
    return(fixed + random + 1)
  }
}

# Insert each model here
calc_K(HD.model7) 



#Final model
HD.model.final <- lmer(hdate ~ region * HW + (1|bay) + (1|year), data = all.fish.18.23, REML = TRUE) 

summary(HD.model.final)

check_model(HD.model.final)

hist(residuals(HD.model.final), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)

r.squaredGLMM(HD.model.final)

Anova(HD.model.final, type = "III")




#--------------------Marginal Means---------------------------------------------

mmwestgulf <- emmeans(HD.model.final, ~ HW * region, pbkrtest.limit = 5839, lmerTest.limit = 5839)
mmwestgulf 

mmwestgulf.HW <- emmeans(HD.model.final, ~ HW, pbkrtest.limit = 5839, lmerTest.limit = 5839)
mmwestgulf.HW

mmwestgulf.region <- emmeans(HD.model.final, ~ region, pbkrtest.limit = 5839, lmerTest.limit = 5839)
mmwestgulf.region




#--------------Marginal Means Plot----------------------------------------------

mmwestgulf.df <- as.data.frame(mmwestgulf)

mmwestgulf.df <- mmwestgulf.df %>%
  mutate(dummy.hdate = as.Date(mmwestgulf.df$emmean - 1, origin = "2020-01-01"),
         dummy.lwr.Cl = as.Date(mmwestgulf.df$lower.CL - 1, origin = "2020-01-01"),
         dummy.upr.Cl = as.Date(mmwestgulf.df$upper.CL - 1, origin = "2020-01-01"),
         dataset = "West Gulf")

eastern.kodiak.07.23 <- eastern.kodiak.07.23 %>%
  mutate(dummy.hdate = as.Date(eastern.kodiak.07.23$emmean - 1, origin = "2020-01-01"),
         dummy.lwr.Cl = as.Date(eastern.kodiak.07.23$lower.CL - 1, origin = "2020-01-01"),
         dummy.upr.Cl = as.Date(eastern.kodiak.07.23$upper.CL - 1, origin = "2020-01-01"),
         region = "Eastern Kodiak",
         dataset = "Eastern Kodiak")


# Combine the datasets
combined.mm <- bind_rows(eastern.kodiak.07.23, mmwestgulf.df)

# Create an interaction variable for positioning
combined.mm <- combined.mm %>%
  mutate(bay.dataset = interaction(region, dataset, sep = "_"),
         HW = factor(HW, levels = c("Before", "During", "Adjacent","After")),
         bay.dataset = factor(bay.dataset, levels = c("Eastern Kodiak_Eastern Kodiak", "Eastern Kodiak_West Gulf", 
                              "Kaiugnak Bay_West Gulf", "Agripina Bay_West Gulf", "Shumagin Islands_West Gulf")))



#Plot of marginal means for both the temporal model and the spatial model

ggplot(combined.mm, aes(bay.dataset, dummy.hdate, color = HW)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = dummy.lwr.Cl, ymax = dummy.upr.Cl), linewidth = 1, width = 0.3, 
                position = position_dodge(width = 0.5)) +
  scale_y_date(date_breaks = "1 month", date_labels = "%b-%d", limits = as.Date(c("2020-02-15", "2020-05-15"))) +
  scale_x_discrete(labels = function(x) gsub("_.*", "", x)) +  # Remove dataset suffix from labels
  theme_classic() +
  ylab("Hatch Date") +
  xlab("Region") +
  labs(color = NULL) +
  scale_color_manual(values = HW.colors) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 24),
        #axis.text.x = element_text(hjust = 1, angle = 45),
        legend.text = element_text(size = 24),
        legend.position = "top",
        strip.text = element_text(size = 24)) +
  guides(color = guide_legend(nrow = 1))
  



