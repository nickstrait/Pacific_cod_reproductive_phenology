#R version 4.5.3, code was last run on 5/19/2026
#Final code for Strait et al. 2026
# "Climate extremes reshape reproductive phenology: Regional variability of hatch timing in the subarctic North Pacific"

#Creates annual age-length keys with a linear model

# Paired with:

#Data from aged individuals including increments in "wide" format
# "1_All_aged_Pcod.csv"

#Collection and size information from all Pcod juveniles captured (not including aged ones)
# "1_All_unaged_Pcod.csv"


#------------------------Load libraries-----------------------------------------

library(rstudioapi)
library(car)
library(performance)
library(emmeans)
library(tidyverse)




#------------------------Set working directory----------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




#------------------------Import csv files---------------------------------------

# Contains all aged juvenile Pacific cod from 2018 to 2023
Aged.Pcod <- read.csv("1_All_aged_Pcod.csv")

# Contains all unaged juvenile Pacific cod from 2018 to 2023
All.unaged.pcod <- read.csv("1_All_unaged_Pcod.csv")




#------------------------Formatting data----------------------------------------

Aged.Pcod$year <- as.factor(Aged.Pcod$year)
Aged.Pcod$bay <- as.factor(Aged.Pcod$bay)
Aged.Pcod$region <- as.factor(Aged.Pcod$region)

All.unaged.pcod$year <- as.factor(All.unaged.pcod$year)




#------------------------Exploring Data-----------------------------------------

# Age and fish length correlation
cor(Aged.Pcod$TL_mm, Aged.Pcod$age) 

ggplot(Aged.Pcod, aes(TL_mm, age)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()


# Age x Length relationship per year
ggplot(Aged.Pcod, aes(TL_mm, age, color = as.factor(year))) + 
  geom_point() +
  geom_smooth(method = lm) +
  xlab("TL (mm)") + 
  ylab("Age (days)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_blank())


# Age x Length relationship per region
ggplot(Aged.Pcod, aes(TL_mm, age, color = as.factor(region))) + 
  geom_point() +
  geom_smooth(method = lm) +
  xlab("TL (mm)") + 
  ylab("Age (days)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12),
        legend.text = element_text(size = 10), legend.title = element_blank())


##------------Regional Differences----------------------------------------------

# Used to validate our region groupings and determine whether there were differences in size-at-age among the bays within Eastern Kodiak and Shumagin Islands

# Comparing Eastern Kodiak Bays (Anton Larsen and Cook Bay)
eastkodiak <- Aged.Pcod %>%
  filter(region == "Eastern Kodiak")

e.kodiak.dif <- lm(age ~  TL_mm + year + doc + bay, data = eastkodiak)
Anova(e.kodiak.dif) # no differences between bays in Eastern Kodiak


# Comparing the Shumagin Islands (Baralof Bay, Sand Point, Falmouth Bay)

shumagin.islands <- Aged.Pcod %>%
  filter(region == "Shumagin Islands",
         year %in% c("2018", "2019", "2020", "2021")) # only Falmouth Bay available in 2022 and 2023

shum.isl.dif <- lm(age ~  TL_mm + year + doc + bay, data = shumagin.islands)
Anova(shum.isl.dif) # significant differences between bays in the Shumagin Islands


ggplot(shumagin.islands, aes(TL_mm, age, color = bay)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~year)

# Calculating marginal means of the Shumagin Islands bays

emm <- emmeans(shum.isl.dif, ~ bay)

# Pairwise comparison of each bay

pairwise_results <- pairs(emm, adjust = "tukey")
summary(pairwise_results) 

# Differences arises between Falmouth and Sand Point 
# Checked if difference matters for model in code below




#-----------------------Age length key models-----------------------------------

lm.TL.year.region <- lm(age ~ TL_mm * year + region, data = Aged.Pcod)
lm.TL.year.doc.region.2 <- lm(age ~ TL_mm * year * doc + region, data = Aged.Pcod)
lm.TL.year.doc.region.3 <- lm(age ~ TL_mm + year * doc + region, data = Aged.Pcod)
lm.TL.year.doc.region.4 <- lm(age ~ TL_mm * year + doc + region, data = Aged.Pcod)

check_model(lm.TL.year.region)
check_model(lm.TL.year.doc.region.2)
check_model(lm.TL.year.doc.region.3)
check_model(lm.TL.year.doc.region.4)

#all have high collinearity due to interactions, so removing and comparing models

model.1 <- lm(age ~ TL_mm, data = Aged.Pcod)
model.2 <- lm(age ~ TL_mm + year, data = Aged.Pcod)
model.3 <- lm(age ~ TL_mm + region, data = Aged.Pcod)
model.4 <- lm(age ~ TL_mm + doc, data = Aged.Pcod)
model.5 <- lm(age ~ TL_mm + year + region, data = Aged.Pcod)
model.6 <- lm(age ~ TL_mm + year + doc, data = Aged.Pcod)
model.7 <- lm(age ~ TL_mm + year + doc + region, data = Aged.Pcod)
model.8 <- lm(age ~ TL_mm + year + doc + latitude + longitude, data = Aged.Pcod)
model.9 <- lm(age ~ TL_mm + year + latitude + longitude, data = Aged.Pcod)
model.10 <- lm(age ~ TL_mm + year + doc + bay, data = Aged.Pcod)

AIC(model.1, model.2, model.3, model.4, model.5, model.6, model.7, model.8, model.9, model.10)

# Model 10 is best, but due to unbalanced sampling of bays, going with the region grouping - model 7

summary(model.7) # 0.8585

check_model(model.7)

hist(residuals(model.7, type = 'pearson'), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)

Anova(model.7, type = "II")




#--------------------Verifying Region Grouping----------------------------------
# Due to differences in size-at-age in the Shumagin Islands - checking impact on model predictions

# Model with only Sand Point samples
model.11 <- lm(age ~ TL_mm + year + doc , data = Aged.Pcod %>% filter(bay == "Sand Point"))
check_model(model.11)

# Model with only Baralof Bay samples
model.12 <- lm(age ~ TL_mm + year + doc , data = Aged.Pcod %>% filter(bay == "Baralof Bay"))
check_model(model.12)

# Model with only Falmouth Bay samples in years that overlap with Sand Point and Baralof Bay
model.13 <- lm(age ~ TL_mm + year + doc , data = Aged.Pcod %>% filter(bay == "Falmouth Bay", year %in% c("2018", "2019", "2020", "2021")))
check_model(model.13) # high collinearity is due to most years only having one doc

summary(model.11)
summary(model.12)
summary(model.13)


# Creating a test dataframe with unique IDs for comparison
test.data <- All.unaged.pcod %>%
  mutate(unique.ID = row_number())

sandpoint <- test.data %>%
  filter(bay == "Sand Point",
         !year %in% c(2022, 2023))

baralof <- test.data %>%
  filter(bay == "Baralof Bay",
         !year %in% c(2022, 2023))

falmouth <- test.data %>%
  filter(bay == "Falmouth Bay",
         !year %in% c(2022, 2023))


#Predicting the ages using four different models
test.data$pred.age.region <- predict(model.7, newdata = test.data) # model with region grouping

test.data$pred.age.bay <- predict(model.10, newdata = test.data) # model with bay grouping

sandpoint$pred.age.SP <- predict(model.11, newdata = sandpoint) # model using only Sand Point 

baralof$pred.age.B <- predict(model.12, newdata = baralof) # model using only Baralof Bay

falmouth$pred.age.F <- predict(model.13, newdata = falmouth) # model using only Falmouth Bay


# Joining predicted ages
sp.baralof.fal <- test.data %>%
  filter(bay %in% c("Sand Point", "Baralof Bay", "Falmouth Bay")) %>%
  select(bay, year, doc, TL_mm, unique.ID, pred.age.region, pred.age.bay) %>%
  full_join(falmouth %>% select(unique.ID, pred.age.F), by = "unique.ID") %>%
  left_join(sandpoint %>% select(unique.ID, pred.age.SP), by = "unique.ID") %>%
  left_join(baralof %>% select(unique.ID, pred.age.B), by = "unique.ID")

# calculating percent difference between the predicted ages
sandpoint.perdiff <- sp.baralof.fal %>%
  filter(bay == "Sand Point") %>%
  mutate(pd.mod7.mod10 = (pred.age.bay - pred.age.region) / pred.age.region *100,
         pd.mod10.mod11 = (pred.age.SP - pred.age.bay) / pred.age.bay * 100,
         pd.mod7.mod11 = (pred.age.SP - pred.age.region) / pred.age.bay * 100)

baralof.perdiff <- sp.baralof.fal %>%
  filter(bay == "Baralof Bay") %>%
  mutate(pd.mod7.mod10 = (pred.age.bay - pred.age.region) / pred.age.region *100,
         pd.mod10.mod12 = (pred.age.B - pred.age.bay) / pred.age.bay * 100,
         pd.mod7.mod12 = (pred.age.B - pred.age.region) / pred.age.bay * 100)

falmouth.perdiff <- sp.baralof.fal %>%
  filter(bay == "Falmouth Bay") %>%
  mutate(pd.mod7.mod10 = (pred.age.bay - pred.age.region) / pred.age.region *100,
         pd.mod10.mod13 = (pred.age.F - pred.age.bay) / pred.age.bay * 100,
         pd.mod7.mod13 = (pred.age.F - pred.age.region) / pred.age.bay * 100)


#average % difference between the models
mean(sandpoint.perdiff$pd.mod7.mod10)
mean(sandpoint.perdiff$pd.mod7.mod11, na.rm = TRUE)
mean(sandpoint.perdiff$pd.mod10.mod11, na.rm = TRUE)
mean(baralof.perdiff$pd.mod7.mod10, na.rm = TRUE)
mean(baralof.perdiff$pd.mod7.mod12, na.rm = TRUE)
mean(baralof.perdiff$pd.mod10.mod12, na.rm = TRUE)
mean(falmouth.perdiff$pd.mod7.mod10, na.rm = TRUE)
mean(falmouth.perdiff$pd.mod7.mod13, na.rm = TRUE)
mean(falmouth.perdiff$pd.mod10.mod13, na.rm = TRUE)
# All average percent differences are less than 5% 



# Quantifying how different the model outputs are
# Calculating the average difference in age predictions from the models

# Falmouth Bay sample predictions
falmouth.perdiff <- falmouth.perdiff %>%
  mutate(age.diff.1 = pred.age.region - pred.age.bay, # region model prediction - bay model prediction
         age.diff.2 = pred.age.region - pred.age.F,   # region model prediction - falmouth model prediction
         age.diff.3 = pred.age.bay - pred.age.F)      # bay model prediction - falmouth model prediction

mean(falmouth.perdiff$age.diff.1)
mean(falmouth.perdiff$age.diff.2, na.rm = TRUE)
mean(falmouth.perdiff$age.diff.3, na.rm = TRUE)
# There is a 2 to 4 day difference on average between all model predicted ages for Falmouth Bay



# Baralof Bay sample predictions
baralof.perdiff <- baralof.perdiff %>%
  mutate(age.diff.1 = pred.age.region - pred.age.bay, # region model prediction - bay model prediction
         age.diff.2 = pred.age.region - pred.age.B,   # region model prediction - baralof model prediction
         age.diff.3 = pred.age.bay - pred.age.B)      # bay model prediction - baralof model prediction

mean(baralof.perdiff$age.diff.1)
mean(baralof.perdiff$age.diff.2, na.rm = TRUE)
mean(baralof.perdiff$age.diff.3, na.rm = TRUE)
# There is less than a day difference on average between all model predicted ages for Baralof Bay



#Sand Point sample predictions
sandpoint.perdiff <- sandpoint.perdiff %>%
  mutate(age.diff.1 = pred.age.region - pred.age.bay,
         age.diff.2 = pred.age.region - pred.age.SP,
         age.diff.3 = pred.age.bay - pred.age.SP)

mean(sandpoint.perdiff$age.diff.1)
mean(sandpoint.perdiff$age.diff.2, na.rm = TRUE)
mean(sandpoint.perdiff$age.diff.3, na.rm = TRUE)
# There is <1 to 3 day difference on average between all model predicted ages for Sand Point


#Based on the average % difference and average difference in ages, we will use the model with region grouping (model 7)




#-----------------------Predicting age for unaged fish--------------------------

predicted.ages <- All.unaged.pcod


# This predicts using the model, but uses confidence intervals to randomly select values between them

conf.dist <- predict(model.7, newdata = predicted.ages, interval = "confidence", level = .95)

head(conf.dist)

predicted.ages[c("fit","lwr.conf", "upr.conf")] <- conf.dist

set.seed(25)
predicted.ages <- predicted.ages %>%
  mutate(predict.age = round(runif(n(), min = lwr.conf, max = upr.conf ))) %>%
  select(-fit, -upr.conf, -lwr.conf)

ggplot(predicted.ages, aes(predict.age, TL_mm)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()



# Bind aged and predicted age dataframes and create HD based off of both aged samples and predicted age
all.predicted.and.aged <- Aged.Pcod %>%
  select(ID, year, month, day, doc, bay, region, latitude, longitude, TL_mm, age) %>%
  mutate(predict.age = NA)

predicted.ages <- predicted.ages %>%
  mutate(age = NA,
         ID = NA)

all.predicted.and.aged <- all.predicted.and.aged %>%
  rbind(predicted.ages) %>%
  mutate(final.age = coalesce(age, predict.age),    # final age combines predicted ages with fish that were aged
         hdate = doc - final.age)

#write_csv(all.predicted.and.aged, "2_All_fish_ages_applied_agelengthkey.csv")




#---------------Length at Age Plot----------------------------------------------

region.colors <- c("Eastern Kodiak" = "blue", "Kaiugnak Bay" = "forestgreen", "Agripina Bay" = "gold2", 

                                      "Shumagin Islands" = "turquoise1")
# Setting order of regions
custom.order <- c("Eastern Kodiak", "Kaiugnak Bay", "Agripina Bay","Shumagin Islands")

Aged.Pcod$region <- fct_relevel(Aged.Pcod$region , custom.order)

# This plot only uses samples we aged (no predicted ages)

ggplot(Aged.Pcod, aes(TL_mm, age, color = region)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  facet_wrap(~year) +
  theme_classic() +
  labs(x = "TL (mm)",
       y = "Age (days)",
       color = NULL) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.position = "top",
        strip.text = element_text(size = 30)) +
  scale_color_manual(values = region.colors)

#ggsave("Total Length x Age - Aged Samples.jpg", width = 20, height = 10, dpi = 600)





