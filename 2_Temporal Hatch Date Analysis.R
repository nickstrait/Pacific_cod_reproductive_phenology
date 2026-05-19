#R version 4.5.3, code was last run on 5/19/2026
#Final code for Strait et al. 2026
# "Climate extremes reshape reproductive phenology: Regional variability of hatch timing in the subarctic North Pacific"


#LMM to examine hatch dates at Eastern Kodiak before, during, adjacent to, and after MHWs

# Paired with:

#Data from all samples, including predicted ages - created in "1_Age-Length Key.R" file
# "2_All_fish_ages_applied_agelengthkey.csv"

#Samples from Almeida et al. 2024 that were included in this analysis
# "2_Almeidaetal_data.csv"



#------------------Load Libraries-----------------------------------------------

library(rstudioapi)
library(car)
library(lme4)
library(performance)
library(patchwork)
library(grid)
library(gridExtra)
library(emmeans)
library(lubridate)
library(ggh4x)
library(MuMIn)
library(tidyverse)




#-----------------Set Working Directory-----------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




#-----------------Import csv files----------------------------------------------

all.fish.18.23 <- read.csv("2_All_fish_ages_applied_agelengthkey.csv")
all.fish.07.17 <- read.csv("2_Almeidaetal_data.csv")




#----------------Data Formatting------------------------------------------------

# Define heatwave categories for 2018 to 2023 data
hw.cats <- tibble(
  year = 2018:2023,
  HW   = c("Adjacent", "During", "Adjacent", "After", "After", "After")
)

# Apply across 2018-2023 years
all.fish.18.23.kodiak <- all.fish.18.23 %>%
  left_join(hw.cats, by = "year") %>%
  mutate(
    hatch_DATE        = as.Date(hdate, origin = paste0(year - 1, "-12-31")),
    hatch_DATE_noYear = format(hatch_DATE, "%b-%d")
  ) %>%
  filter(region == "Eastern Kodiak") %>%
  transmute(
    ID = as.character(ID),
    year, bay, hdate,
    hatch_DATE = ymd(hatch_DATE),
    hatch_DATE_noYear = format(hatch_DATE, "%b-%d"),
    final.age, TL_mm,
    HW = factor(HW)
  )


# Define heatwave categories for 2007 to 2017 data
hw.cat.07.17 <- tibble(
  year   = c(2007, 2009, 2010, 2012, 2013, 2014, 2015, 2016, 2017),
  origin = as.Date(paste0(c(2006, 2008, 2009, 2011, 2012, 2013, 2014, 2015, 2016), "-12-31")),
  HW     = c(rep("Before", 6), rep("During", 3))
)

# Apply across 2007–2017 years
all.fish.07.17 <- all.fish.07.17 %>%
  inner_join(hw.cat.07.17, by = "year") %>%
  mutate(
    hdate            = doc - final.age,
    hatch_DATE       = as.Date(hdate, origin = origin),
    hatch_DATE_noYear = format(hatch_DATE, "%b-%d")
  ) %>%
  select(ID, year, bay, hdate, hatch_DATE, hatch_DATE_noYear, final.age, TL_mm, HW) %>%
  mutate(ID = as.character(ID))


  
#Binding datasets into one
kodiak.07.23 <- bind_rows(all.fish.18.23.kodiak, all.fish.07.17)  

# Creating a "fake" date to help create figures
kodiak.07.23 <- kodiak.07.23 %>%
  mutate(year = as.factor(year),
         bay = as.factor(bay),
         dummy_date = as.Date(paste0("2000-", kodiak.07.23$hatch_DATE_noYear), format = "%Y-%b-%d")) 


# Specifying the order of bays and HWs
custom.order1 <- c("Anton Larsen Bay", "Cook Bay")

kodiak.07.23$bay <- fct_relevel(kodiak.07.23$bay, custom.order1)

custom.order2 <- c("Before", "During", "Adjacent", "After")

kodiak.07.23$HW <- fct_relevel(kodiak.07.23$HW , custom.order2)




#----------------Defining Colors for plots--------------------------------------

bay.colors <- c("Anton Larsen Bay" = "blue", "Cook Bay" = "seagreen2")
HW.colors <- c("Before" = "darkblue", "During" = "darkred", "Adjacent" = "lightpink", "After" = "steelblue1")
year.colors <- c("2007" = "blue", "2009" = "red", "2010" = "forestgreen", "2012" = "gold2", "2013" = "turquoise1",
                  "2014" = "lightpink", "2015" = "blue", "2016" = "red", "2017" = "forestgreen", "2019" = "gold2",
                  "2018" = "blue", "2020" = "red", "2021" = "blue", 
                  "2022" = "red", "2023" = "forestgreen")
year.colors2 <- c("2007" = "darkblue", "2009" = "darkblue", "2010" = "darkblue", "2012" = "darkblue", "2013" = "darkblue",
                  "2014" = "darkblue", "2015" = "darkred", "2016" = "darkred", "2017" = "darkred", "2019" = "darkred",
                  "2018" = "lightpink", "2020" = "lightpink", "2021" = "steelblue1", 
                  "2022" = "steelblue1", "2023" = "steelblue1")




#-------------------Temporal Models---------------------------------------------

# Determining Random Effects
kodiak.model.1 <- lmer(hdate ~ HW * bay + (1 | year), data = kodiak.07.23, REML = FALSE) 
kodiak.model.2 <- lm(hdate ~ HW * bay, data = kodiak.07.23)

AIC(kodiak.model.1, kodiak.model.2)

# model.1 has best random effects

Weights(AIC(kodiak.model.1, kodiak.model.2))


# Determining Fixed Effects
kodiak.model.3 <- lmer(hdate ~ HW + bay + (1 | year), data = kodiak.07.23, REML = FALSE)
kodiak.model.4 <- lmer(hdate ~ HW + (1 | year), data = kodiak.07.23, REML = FALSE)

AIC(kodiak.model.1, kodiak.model.3, kodiak.model.4)

# model.1 has best fixed effects
summary(kodiak.model.1)


# r squared and weights of the models
r.squaredGLMM(kodiak.model.1)
r.squaredGLMM(kodiak.model.2)
r.squaredGLMM(kodiak.model.3)
r.squaredGLMM(kodiak.model.4)

Weights(AIC(kodiak.model.1, kodiak.model.2, kodiak.model.3, kodiak.model.4))


# Function to determine K for each model
calc_K <- function(model) {
  if(inherits(model, "lm")) {
    return(length(coef(model)) + 1)
  } else if(inherits(model, "lmerMod")) {
    fixed <- length(fixef(model))
    random <- length(VarCorr(model))  
    return(fixed + random + 1)
  }
}

# Insert each model here to determine
calc_K(kodiak.model.1) 


# Checking model 1

check_model(kodiak.model.1) #High collinearity issues 
vif(kodiak.model.1)

#this is due to how the model handles two factors with an interaction. 
#It creates dummies for the factors and the interaction column is a direct product of that so it is structurally collinear, 
#but not isn't an issue as long as the design is balanced and SE aren't large



#Final Model
kodiak.model.final <- lmer(hdate ~ HW * bay + (1 | year), data = kodiak.07.23, REML = TRUE)

summary(kodiak.model.final)

check_model(kodiak.model.final)

vif(kodiak.model.final)

hist(residuals(kodiak.model.final), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)

r.squaredGLMM(kodiak.model.final)

Anova(kodiak.model.final, type = "III")

emmeans(kodiak.model.final, pairwise ~ bay | HW, adjust = "tukey") # bays statistically differ only in Adjacent and After




#-----------------Marginal Means------------------------------------------------

mmkodiak.hw <- emmeans(kodiak.model.final, ~ HW)
mmkodiak.hw

mmkodiak.bay <- emmeans(kodiak.model.final, ~ bay)
mmkodiak.bay

mmkodiak.region.df <- as.data.frame(mmkodiak.hw)

mmkodiak <- emmeans(kodiak.model.final, ~ bay * HW)
mmkodiak

# This file will be used in "3_Spatial Hatch Date Analysis.R"

#write_csv(mmkodiak.region.df, "3_Eastern Kodiak Marginal Means 07-23.csv")




#-----------------Plots---------------------------------------------------------

#Creating step up for coloring facet wrap headers

strip <- strip_themed(
  text_x = elem_list_text(colour = year.colors2[as.character(sort(unique(kodiak.07.23$year)))])
)



#Hatch Date Densities by each bay and year
ggplot(kodiak.07.23, aes(as.Date(dummy_date, format = "%b-%d"), fill = bay)) + 
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap2(~year, strip = strip) +
  xlab("Hatch Date") + 
  ylab("Density") +
  geom_vline(xintercept = as.Date("2000-04-09"),
             color = "black",
             linewidth = 0.75,
             linetype = "dashed") +
  labs(fill = NULL) +
  scale_x_date(breaks = as.Date(c("2000-02-01", "2000-04-01", "2000-06-01")), date_labels = "%b-%d", limits = as.Date(c("2000-01-15", "2000-06-15"))) +
  scale_fill_manual(values = bay.colors) +
  theme_classic() +
  theme(axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.position = "top",
        strip.text = element_text(size = 30)) +
  guides(fill = guide_legend(nrow = 1))

#The warning message is due to only one sample from Anton Larsen Bay in 2015














