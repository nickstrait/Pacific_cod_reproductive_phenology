#R version 4.5.1, code was last run on X/X/2026
#Final code for Climate extremes reshape reproductive phenology: Regional variability of hatch timing in the subarctic North Pacific

#Code to calculate hatch success based of deep water temperatures

# Paired with:

#Data from all samples, including predicted ages - created in "1_Age-Length Key.R" file
# "2_All_fish_ages_applied_agelengthkey.csv"

#Samples from Almeida et al. 2024 that were included in this analysis
# "2_Almeidaetal_data.csv"

#Temperature Data from Shumagin Islands mooring - temps were predicted using GAK1
# "4_Shumagin Predicted Temps.csv"

#Average prior 22-day temperature in Before years for each GAK1 depth set
# "4_before_HW_incubation_temps.csv"

#Average prior 22-day temperature in Adjacent years for each GAK1 depth set
# "4_adjacent_HW_incubation_temps.csv"

#Average prior 22-day temperature in During years for each GAK1 depth set
# "4_during_HW_incubation_temps.csv"

#Average prior 22-day temperature in After years for each GAK1 depth set
# "4_after_HW_incubation_temps.csv"

#Average daily temperature from GAK1 50 to 115 m depths
# "4_GAK Mooring 50-115m averages.csv"

https://doi.org/10.5281/zenodo.18623697


#------------------Load Libraries-----------------------------------------------

library(rstudioapi)
library(grid)
library(gridExtra)
library(patchwork)
library(slider)
library(grid)
library(tidyverse)




#-----------------Set Working Directory-----------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))




#-----------------Import csv files----------------------------------------------

Shumagin <- read.csv("4_Shumagin Predicted Temps.csv")
before.incub.temps <- read.csv("4_before_HW_incubation_temps.csv")
adjacent.incub.temps <- read.csv("4_adjacent_HW_incubation_temps.csv")
during.incub.temps <- read.csv("4_during_HW_incubation_temps.csv")
after.incub.temps <- read.csv("4_after_HW_incubation_temps.csv")
GAK.50.115.allyrs <- read.csv("4_GAK Mooring 50-115m averages.csv")
all.fish.18.23 <- read.csv("2_All_fish_ages_applied_agelengthkey.csv")
all.fish.07.17 <- read.csv("2_Almeidaetal_data.csv")




#----------------Data Formatting------------------------------------------------

Shumagin <- Shumagin %>%
  mutate(date = ymd(date),
         hatch_DATE = date,
         doy = as.numeric(format(date, "%j")),
         dummy_date = as.Date(doy - 1, origin = "2000-01-01")) %>%
  select(-temp, -shum.predict) %>%
  rename(temp = temp.final)


#Formatting Temperature Data
GAK.50.115.allyrs <- GAK.50.115.allyrs %>%
  mutate(date = ymd(date),
         doy = yday(date),
         HW = case_when(
           year(date) < 2015 ~ "Before",
           year(date) %in% c(2015, 2016, 2017, 2019) ~ "During",
           year(date) %in% c(2018, 2020) ~ "Adjacent",
           between(year(date), 2021, 2023) ~ "After"),
         month.day = format(date, "%m-%d")) %>%
  filter(!(month(date) == 2 & day(date) == 29))


# Define heatwave categories for 2018 to 2023 data
hw.cats <- tibble(
  year = 2018:2023,
  HW   = c("Adjacent", "During", "Adjacent", "After", "After", "After")
)

# Apply across 2018-2023 years
all.fish.18.23 <- all.fish.18.23 %>%
  left_join(hw.cats, by = "year")

all.fish.18.23.kodiak <- all.fish.18.23 %>%
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


#Bind two dataframes
kodiak.07.23 <- bind_rows(all.fish.18.23.kodiak, all.fish.07.17)  

kodiak.07.23 <- kodiak.07.23 %>%
  mutate(year = as.factor(year),
         bay = as.factor(bay),
         dummy_date = as.Date(hdate - 1, origin = "2000-01-01")) # this creates a "fake" date to help create figures



#Creating separate data frames for each HW category for plotting later
kodiak.07.23.before <- kodiak.07.23 %>%
  filter(HW == "Before")

kodiak.07.23.adjacent <- kodiak.07.23 %>%
  filter(HW == "Adjacent")

kodiak.07.23.during <- kodiak.07.23 %>%
  filter(HW == "During")

kodiak.07.23.after <- kodiak.07.23 %>%
  filter(HW == "After")




#----------Calculating the 22 day average temperatures and daily hatch success-------------------------


# Using GAK1 for Eastern Kodiak
GAK.50.115.allyrs.hs <- GAK.50.115.allyrs %>%
  arrange(date) %>%
  mutate(avg.temp.22d = slide_index_dbl(
    .x = temp,
    .i = date,
    .f = ~ mean(.x, na.rm = TRUE),
    .before = 23,     # Looks back 23 days (including current date)
    .after = -1       # Excludes the current date
  ),
  year = as.factor(year(date)),
  hatch.success = (1 / (1 + ((avg.temp.22d - 4.192) / 2.125)^2)),
  dummy_date = as.Date(doy - 1, origin = "2000-01-01"))


# Creating different dataframes for GAK1 for each HW category
GAK.50.115.before.hs <- GAK.50.115.allyrs.hs %>%
  filter(HW == "Before",
         year %in% c(2007, 2009, 2010, 2012, 2013, 2014))

GAK.50.115.adjacent.hs <- GAK.50.115.allyrs.hs %>%
  filter(HW == "Adjacent")

GAK.50.115.during.hs <- GAK.50.115.allyrs.hs %>%
  filter(HW == "During")

GAK.50.115.after.hs <- GAK.50.115.allyrs.hs %>%
  filter(HW == "After")


# Using Shumagin Buoy for Shumagin Islands
Shumagin <- Shumagin %>%
  arrange(date) %>%
  mutate(avg.temp.22d = slide_index_dbl(
    .x = temp,
    .i = date,
    .f = ~ mean(.x, na.rm = TRUE),
    .before = 23,     # Looks back 23 days (including current date)
    .after = -1       # Excludes the current date
  )) %>%
  mutate(year = year(date),
         HW = case_when(year %in% c(2015, 2016, 2017, 2019) ~ "During",
                        year %in% c(2007, 2009, 2010, 2012, 2013, 2014) ~ "Before", 
                        year %in% c(2018, 2020) ~ "Adjacent",
                        year %in% c(2021, 2022, 2023) ~ "After"),
         year = as.factor(year),
         doy = yday(date),
         hatch.success = (1 / (1 + ((avg.temp.22d - 4.192) / 2.125)^2)))


# Creating different dataframes for Shumagin Buoy for each HW category
Shumagin.hs.before <- Shumagin %>%
  filter(HW == "Before")

Shumagin.hs.during <- Shumagin %>%
  filter(HW == "During")

Shumagin.hs.adjacent <- Shumagin %>%
  filter(HW == "Adjacent")

Shumagin.hs.after <- Shumagin %>%
  filter(HW == "After")




#------------Calculating the daily hatch success for each HW category-----------

#These data will be used to create hatch success by depth plot

# Before 

before.incub.temps.HS <- before.incub.temps %>%
  mutate(hatch.success = (1 / (1 + ((avg.22.temp - 4.192) / 2.125)^2)),
         dummy_date = mdy(dummy_date)) 

# During

during.incub.temps.HS <- during.incub.temps %>%
  mutate(hatch.success = (1 / (1 + ((avg.22.temp - 4.192) / 2.125)^2)),
         dummy_date = mdy(dummy_date)) 

# Adjacent

adjacent.incub.temps.HS <- adjacent.incub.temps %>%
  mutate(hatch.success = (1 / (1 + ((avg.22.temp - 4.192) / 2.125)^2)),
         dummy_date = mdy(dummy_date))

# After

after.incub.temps.HS <- after.incub.temps %>%
  mutate(hatch.success = (1 / (1 + ((avg.22.temp - 4.192) / 2.125)^2)),
         dummy_date = mdy(dummy_date))


# Specify the order of the GAK depths
custom.order <- c("GAK1 50 to 115 m", "GAK1 120 to 185 m", "GAK1 190 to 250 m")

# Convert the variable to a factor with custom order using forcats
before.incub.temps.HS$depth <- fct_relevel(before.incub.temps.HS$depth, custom.order)
during.incub.temps.HS$depth <- fct_relevel(during.incub.temps.HS$depth, custom.order)
adjacent.incub.temps.HS$depth <- fct_relevel(adjacent.incub.temps.HS$depth, custom.order)
after.incub.temps.HS$depth <- fct_relevel(after.incub.temps.HS$depth, custom.order)




#--------------Defining colors for plots----------------------------------------

year.colors <- c("2007" = "blue", "2009" = "red", "2010" = "forestgreen", "2012" = "gold2", "2013" = "turquoise1",
                  "2014" = "lightpink", "2015" = "blue", "2016" = "red", "2017" = "forestgreen", "2019" = "gold2",
                  "2018" = "blue", "2020" = "red", "2021" = "blue", 
                  "2022" = "red", "2023" = "forestgreen")

year.colors2 <- c("2018" = "blue", "2019" = "red", "2020" = "red", "2021" = "blue", "2022" = "red", "2023" = "forestgreen")

depth.colors <- c("GAK1 50 to 115 m" = "lightskyblue1", "GAK1 120 to 185 m" = "dodgerblue2", "GAK1 190 to 250 m" = "darkblue")




#-----------------------Eastern Kodiak Hatch Success Plots----------------------

# The following sections create each hatch success for Eastern Kodiak and the Shumagin Islands

##---------------Before---------------------------------------------------------

doy.all.before <- as.numeric(format(kodiak.07.23.before$dummy_date, "%j"))
density.all.before <- density(doy.all.before)
max.density.all.before <- max(density.all.before$y)

density.by.year.df.before <- kodiak.07.23.before %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    d <- density(.$doy)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.before <- max(density.by.year.df.before$y)
max.density.before <- max(max.density.all.before, max.density.by.year.before)

# 2. Scale all density curves so their max is, 1
density.scale.before <- 1 / max.density.before

density.all.df.before <- data.frame(
  x = as.Date(density.all.before$x, origin = "2000-01-01"),
  y = density.all.before$y * density.scale.before
)
density.by.year.df.before$y <- density.by.year.df.before$y * density.scale.before

# 3. Plot
before.kodiak.hs.plot <- ggplot(GAK.50.115.before.hs, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.before, aes(x = x, y = y, fill = year), alpha = 0.3, position = "identity", inherit.aes = FALSE) +
  geom_area(data = density.all.df.before, aes(x = x, y = y), fill = "black",alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * max.density.before)
  ) +
  labs(x = NULL, color = "Before", fill = "Before", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)


# Function to find the narrowest interval around the peak density containing 95% proportion of the density
find.centered.range <- function(df, target.prop) {
  y <- df$y
  dates <- df$x
  total <- sum(y)
  
  peak <- which.max(y)
  left <- peak
  right <- peak
  current_sum <- y[peak]
  n <- length(y)
  
  while (current_sum / total < target.prop) {
    expand_left <- left > 1
    expand_right <- right < n
    
    if (expand_left & expand_right) {
      if (y[left - 1] > y[right + 1]) {
        left <- left - 1
        current_sum <- current_sum + y[left]
      } else {
        right <- right + 1
        current_sum <- current_sum + y[right]
      }
    } else if (expand_left) {
      left <- left - 1
      current_sum <- current_sum + y[left]
    } else if (expand_right) {
      right <- right + 1
      current_sum <- current_sum + y[right]
    } else {
      break
    }
  }
  
  return(c(dates[left], dates[right]))
}

###------ Date range where 95% of data falls------------------------------------
kodiak.range_95.before <- find.centered.range(density.all.df.before, 0.95)
kodiak.range_95.before

# March 12 to May 5th - Falls within area of highest hatch success


###---- Determining how many days on average were above 75% before MHW (Using the geom_smooth line)----------------

# Fit the model 
fit.kod.before <- loess(hatch.success ~ as.numeric(dummy_date), data = GAK.50.115.before.hs, span = 0.75)

# Generate a date sequence for the full range
daily.dates.kod.before <- tibble(
  date = seq(min(GAK.50.115.before.hs$dummy_date), max(GAK.50.115.before.hs$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.kod.before <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.kod.before$predicted.hatch.success <- predict(fit.kod.before, newdata = daily.dates.kod.before)

# Count days with hatch success above 75%
sum(daily.dates.kod.before$predicted.hatch.success > 0.75, na.rm = TRUE)

# 151 days above 75%




##------------Adjacent----------------------------------------------------------

doy.all.adjacent <- as.numeric(format(kodiak.07.23.adjacent$dummy_date, "%j"))
density.all.adjacent <- density(doy.all.adjacent)
max.density.all.adjacent <- max(density.all.adjacent$y)

density.by.year.df.adjacent <- kodiak.07.23.adjacent %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    d <- density(.$doy)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.adjacent <- max(density.by.year.df.adjacent$y)
max.density.adjacent <- max(max.density.all.adjacent, max.density.by.year.adjacent)

# 2. Scale all density curves so their max is, 1
density.scale.adjacent <- 1 / max.density.adjacent

density.all.df.adjacent <- data.frame(
  x = as.Date(density.all.adjacent$x, origin = "2000-01-01"),
  y = density.all.adjacent$y * density.scale.adjacent
)
density.by.year.df.adjacent$y <- density.by.year.df.adjacent$y * density.scale.adjacent

# 3. Plot
adjacent.kodiak.hs.plot <- ggplot(GAK.50.115.adjacent.hs, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.adjacent, aes(x = x, y = y, fill = year), alpha = 0.3, position = "identity", inherit.aes = FALSE) +
  geom_area(data = density.all.df.adjacent, aes(x = x, y = y), fill = "black", alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * max.density.adjacent)
  ) +
  labs(x = NULL, color = "Adjacent", fill = "Adjacent", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32),
    axis.text.y.left = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y.left = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)




###------ Date range where 95% of data falls-----------------
kodiak.range_95.adj <- find.centered.range(density.all.df.adjacent, 0.95)
kodiak.range_95.adj 

# Feb 19 to April 7th - Falls within area with highest hatch success


###---- Determining how many days on average were above 75% adjacent to MHW (Using the geom_smooth line)---------------
# Fit the model 
fit.kod.adj <- loess(hatch.success ~ as.numeric(dummy_date), data = GAK.50.115.adjacent.hs, span = 0.75)

# Generate a date sequence for the full range
daily.dates.kod.adj <- tibble(
  date = seq(min(GAK.50.115.adjacent.hs$dummy_date), max(GAK.50.115.adjacent.hs$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.kod.adj <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.kod.adj$predicted.hatch.success <- predict(fit.kod.adj, newdata = daily.dates.kod.adj)

# Count days with hatch success above 75%
sum(daily.dates.kod.adj$predicted.hatch.success > 0.75, na.rm = TRUE)

# 106 days above 75%




##-------------During-----------------------------------------------------------

doy.all.during <- as.numeric(format(kodiak.07.23.during$dummy_date, "%j"))
density.all.during <- density(doy.all.during)
max.density.all.during <- max(density.all.during$y)

density.by.year.df.during <- kodiak.07.23.during %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    d <- density(.$doy)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.during <- max(density.by.year.df.during$y)
max.density.during <- max(max.density.all.during, max.density.by.year.during)

# 2. Scale all density curves so their max is, 1
density.scale.during <- 1 / max.density.during

density.all.df.during <- data.frame(
  x = as.Date(density.all.during$x, origin = "2000-01-01"),
  y = density.all.during$y * density.scale.during
)
density.by.year.df.during$y <- density.by.year.df.during$y * density.scale.during

# 3. Plot
during.kodiak.hs.plot <- ggplot(GAK.50.115.during.hs, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.during, aes(x = x, y = y, fill = year), alpha = 0.3, position = "identity", inherit.aes = FALSE) +
  geom_area(data = density.all.df.during, aes(x = x, y = y), fill = "black", alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * max.density.during)
  ) +
  labs(x = NULL, color = "During", fill = "During", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32)
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)




###------ Date range where 95% of data falls------------------------------------
kodiak.range_95.during <- find.centered.range(density.all.df.during, 0.95)
kodiak.range_95.during

# Feb 27 to April 16 - Falls within area of highest hatch success


###---- Determining how many days on average were above 75% during to MHW (Using the geom_smooth line)---------------
# Fit the model 
fit.kod.during <- loess(hatch.success ~ as.numeric(dummy_date), data = GAK.50.115.during.hs, span = 0.75)

# Generate a date sequence for the full range
daily.dates.kod.during <- tibble(
  date = seq(min(GAK.50.115.during.hs$dummy_date), max(GAK.50.115.during.hs$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.kod.during <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.kod.during$predicted.hatch.success <- predict(fit.kod.during, newdata = daily.dates.kod.during)

# Count days with hatch success above 75%
sum(daily.dates.kod.during$predicted.hatch.success > 0.75, na.rm = TRUE)

# 0 days above 75%




##-----------------After--------------------------------------------------------

doy.all.after <- as.numeric(format(kodiak.07.23.after$dummy_date, "%j"))
density.all.after <- density(doy.all.after)
max.density.all.after <- max(density.all.after$y)

density.by.year.df.after <- kodiak.07.23.after %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    d <- density(.$doy)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.after <- max(density.by.year.df.after$y)
max.density.after <- max(max.density.all.after, max.density.by.year.after)

# 2. Scale all density curves so their max is, 1
density.scale.after <- 1 / max.density.after

density.all.df.after <- data.frame(
  x = as.Date(density.all.after$x, origin = "2000-01-01"),
  y = density.all.after$y * density.scale.after
)
density.by.year.df.after$y <- density.by.year.df.after$y * density.scale.after

# 3. Plot
after.kodiak.hs.plot <- ggplot(GAK.50.115.after.hs, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.after, aes(x = x, y = y, fill = year), alpha = 0.3, position = "identity", inherit.aes = FALSE) +
  geom_area(data = density.all.df.after, aes(x = x, y = y), fill = "black", alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * max.density.after)
  ) +
  labs(x = NULL, color = "After", fill = "After", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)




###------ Date range where 95% of data falls------------------------------------
kodiak.range_95.after <- find.centered.range(density.all.df.after, 0.95)
kodiak.range_95.after

# Feb 17 to May 2nd - Falls within area of highest hatch success


###---- Determining how many days on average were above 75% after to MHW (Using the geom_smooth line)---------------
# Fit the model 
fit.kod.after <- loess(hatch.success ~ as.numeric(dummy_date), data = GAK.50.115.after.hs, span = 0.75)

# Generate a date sequence for the full range
daily.dates.kod.after <- tibble(
  date = seq(min(GAK.50.115.after.hs$dummy_date), max(GAK.50.115.after.hs$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.kod.after <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.kod.after$predicted.hatch.success <- predict(fit.kod.after, newdata = daily.dates.kod.adj)

# Count days with hatch success above 75%
sum(daily.dates.kod.after$predicted.hatch.success > 0.75, na.rm = TRUE)

# 125 days above 75%


##------------------Merging plots-----------------------------------------------

# This creates some margins around the plots so when they are spaced out more
before.kodiak.hs.plot <- before.kodiak.hs.plot + theme(plot.margin = margin(b = 15, r = 15))
adjacent.kodiak.hs.plot <- adjacent.kodiak.hs.plot + theme(plot.margin = margin(l = 15, b= 15))
during.kodiak.hs.plot <- during.kodiak.hs.plot + theme(plot.margin = margin(t = 15, r = 15))
after.kodiak.hs.plot <- after.kodiak.hs.plot + theme(plot.margin = margin(t = 15, l = 15))


# This combines the four above plots into one plot
combined.hs.kodiak <- (before.kodiak.hs.plot + adjacent.kodiak.hs.plot) / (during.kodiak.hs.plot + after.kodiak.hs.plot) +
  plot_layout(axes = "collect_x")

#Adding axis labels to the plot
combined.hs.kodiak <- grid.arrange(
  patchworkGrob(combined.hs.kodiak),
  left = textGrob("Hatch Success", rot = 90, gp = gpar(fontsize = 32, fontface = "bold")),
  bottom = textGrob("Day of Year", gp = gpar(fontsize = 32, fontface = "bold"), vjust = 0.3),
  right = textGrob("Hatch Date Density", rot = -90, gp = gpar(fontsize = 32, fontface = "bold"))
)




#------------------- Shumagin Island Hatch Success plot-------------------------

# Creating dataframes for each HW category (no samples from before at Shumagin Islands)
all.fish.18.23.during <- all.fish.18.23 %>%
  filter(HW == "During",
         region == "Shumagin Islands") %>%
  mutate(dummy_date = as.Date(hdate - 1, origin = "2000-01-01"),
         year = as.factor(year))

all.fish.18.23.adjacent <- all.fish.18.23 %>%
  filter(HW == "Adjacent", 
         region == "Shumagin Islands") %>%
  mutate(dummy_date = as.Date(hdate - 1, origin = "2000-01-01"),
         year = as.factor(year))


all.fish.18.23.after <- all.fish.18.23 %>%
  filter(HW == "After",
         region == "Shumagin Islands") %>%
  mutate(dummy_date = as.Date(hdate - 1, origin = "2000-01-01"),
         year = as.factor(year))




##------------ Before ----------------------------------------------------------

# 3. Plot
before.shum.hs.plot <- ggplot(Shumagin.hs.before, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) + 
  geom_col(aes(x = dummy_date, y = 0, fill = year), alpha = 0.2, show.legend = TRUE) + # using to create a better legend
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ .)  # Blank secondary axis
  ) +
  labs(x = NULL, fill = "Before", color = "Before", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32),
    plot.margin = margin(b = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank()
  ) +
  scale_fill_manual(values = year.colors) +
  scale_color_manual(values = year.colors) +
  guides(fill = guide_legend(nrow = 1))




###---- Determining how many days on average were above 75% before to MHW (Using the geom_smooth line)----
# Fit the model 
fit.shum.before <- loess(hatch.success ~ as.numeric(dummy_date), data = Shumagin.hs.before, span = 0.75)

# Generate a date sequence for the full range
daily.dates.shum.before <- tibble(
  date = seq(min(Shumagin.hs.before$dummy_date), max(Shumagin.hs.before$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.shum.before <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.shum.before$predicted.hatch.success <- predict(fit.shum.before, newdata = daily.dates.shum.before)

# Count days with hatch success above 75%
sum(daily.dates.shum.before$predicted.hatch.success > 0.75, na.rm = TRUE)

# 153 days above 75%




##------------ During ----------------------------------------------------------

doy.all.shum.during <- as.numeric(format(all.fish.18.23.during$dummy_date, "%j"))
doy.all.clean.shum.during <- na.omit(doy.all.shum.during)
density.all.shum.during <- density(doy.all.clean.shum.during)
max.density.all.shum.during <- max(density.all.shum.during$y)

density.by.year.df.shum.during <- all.fish.18.23.during %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    doy_no_na <- na.omit(.$doy)  # Remove NAs
    d <- density(doy_no_na)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.shum.during <- max(density.by.year.df.shum.during$y)
max.density.shum.during <- max(max.density.all.shum.during, max.density.by.year.shum.during)

# 2. Scale all density curves so their max is, 1
density.scale.shum.during <- 1 / max.density.shum.during

density.all.df.shum.during <- data.frame(
  x = as.Date(density.all.shum.during$x, origin = "2000-01-01"),
  y = density.all.shum.during$y * density.scale.shum.during
)
density.by.year.df.shum.during$y <- density.by.year.df.shum.during$y * density.scale.shum.during

# 3. Plot
during.shum.hs.plot <- ggplot(Shumagin.hs.during, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +
  geom_col(aes(x = dummy_date, y = 0, fill = year), alpha = 0.2, show.legend = TRUE) + # using to create a better legend
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.shum.during, aes(x = x, y = y, fill = year), alpha = 0.3, position = "identity", 
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_area(data = density.all.df.shum.during, aes(x = x, y = y), fill = "black", alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . / density.scale.shum.during)
  ) +
  labs(x = NULL, color = "During", fill = "During", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32),
    plot.margin = margin(b = 15),
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)




###------ Date range where 95% of data falls------------------------------------
shum.range_95.during <- find.centered.range(density.all.df.shum.during, 0.95)
shum.range_95.during

# March 21 to May 12th - Falls within area of highest hatch success


###---- Determining how many days on average were above 75% during to MHW (Using the geom_smooth line)-------------------
# Fit the model 
fit.shum.during <- loess(hatch.success ~ as.numeric(dummy_date), data = Shumagin.hs.during, span = 0.75)

# Generate a date sequence for the full range
daily.dates.shum.during <- tibble(
  date = seq(min(Shumagin.hs.during$dummy_date), max(Shumagin.hs.during$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.shum.during <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.shum.during$predicted.hatch.success <- predict(fit.shum.during, newdata = daily.dates.shum.during)

# Count days with hatch success above 75%
sum(daily.dates.shum.during$predicted.hatch.success > 0.75, na.rm = TRUE)

# 104 days above 75%




##------------ Adjacent---------------------------------------------------------

doy.all.shum.adjacent <- as.numeric(format(all.fish.18.23.adjacent$dummy_date, "%j"))
doy.all.clean.shum.adjacent <- na.omit(doy.all.shum.adjacent)
density.all.shum.adjacent <- density(doy.all.clean.shum.adjacent)
max.density.all.shum.adjacent <- max(density.all.shum.adjacent$y)

density.by.year.df.shum.adjacent <- all.fish.18.23.adjacent %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    doy_no_na <- na.omit(.$doy)  # Remove NAs
    d <- density(doy_no_na)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.shum.adjacent <- max(density.by.year.df.shum.adjacent$y)
max.density.shum.adjacent <- max(max.density.all.shum.adjacent, max.density.by.year.shum.adjacent)

# 2. Scale all density curves so their max is, 1
density.scale.shum.adjacent <- 1 / max.density.shum.adjacent

density.all.df.shum.adjacent <- data.frame(
  x = as.Date(density.all.shum.adjacent$x, origin = "2000-01-01"),
  y = density.all.shum.adjacent$y * density.scale.shum.adjacent
)
density.by.year.df.shum.adjacent$y <- density.by.year.df.shum.adjacent$y * density.scale.shum.adjacent

# 3. Plot
adjacent.shum.hs.plot <- ggplot(Shumagin.hs.adjacent, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +  # hide line key
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.shum.adjacent,
            aes(x = x, y = y, fill = year),
            colour = NA,      
            alpha = 0.3,
            position = "identity",
            inherit.aes = FALSE) +
  geom_area(data = density.all.df.shum.adjacent,
            aes(x = x, y = y),
            fill = "black", alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(date_labels = "%b-%d",
               date_breaks = "1 month",
               limits = as.Date(c("2000-01-15", "2000-06-15"))) +
  scale_y_continuous(limits = c(0, 1),
                     sec.axis = sec_axis(~ . / density.scale.shum.adjacent)) +
  labs(x = NULL, color = "Adjacent", fill = "Adjacent", y = NULL) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(size = 26, face = "bold"),
        legend.text  = element_text(size = 26),
        axis.text = element_text(size = 26),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank()) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)




###------ Date range where 95% of data falls------------------------------------
shum.range_95.adj <- find.centered.range(density.all.df.shum.adjacent, 0.95)
shum.range_95.adj

# March 20 to May 11th - Falls within area of highest hatch success


###---- Determining how many days on average were above 75% adjacent to MHW (Using the geom_smooth line)-------------------
# Fit the model 
fit.shum.adj <- loess(hatch.success ~ as.numeric(dummy_date), data = Shumagin.hs.adjacent, span = 0.75)

# Generate a date sequence for the full range
daily.dates.shum.adj <- tibble(
  date = seq(min(Shumagin.hs.adjacent$dummy_date), max(Shumagin.hs.adjacent$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.shum.adj <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.shum.adj$predicted.hatch.success <- predict(fit.shum.adj, newdata = daily.dates.shum.adj)

# Count days with hatch success above 75%
sum(daily.dates.shum.adj$predicted.hatch.success > 0.75, na.rm = TRUE)

# 151 days above 75%




##-----------------After--------------------------------------------------------

doy.all.shum.after <- as.numeric(format(all.fish.18.23.after$dummy_date, "%j"))
doy.all.clean.shum.after <- na.omit(doy.all.shum.after)
density.all.shum.after <- density(doy.all.clean.shum.after)
max.density.all.shum.after <- max(density.all.shum.after$y)

density.by.year.df.shum.after <- all.fish.18.23.after %>%
  mutate(doy = as.numeric(format(dummy_date, "%j"))) %>%
  group_by(year) %>%
  do({
    doy_no_na <- na.omit(.$doy)  # Remove NAs
    d <- density(doy_no_na)
    data.frame(
      x = as.Date(d$x, origin = "2000-01-01"),
      y = d$y,
      year = unique(.$year)
    )
  }) %>%
  ungroup()
max.density.by.year.shum.after <- max(density.by.year.df.shum.after$y)
max.density.shum.after <- max(max.density.all.shum.after, max.density.by.year.shum.after)

# 2. Scale all density curves so their max is, 1
density.scale.shum.after <- 1 / max.density.shum.after

density.all.df.shum.after <- data.frame(
  x = as.Date(density.all.shum.after$x, origin = "2000-01-01"),
  y = density.all.shum.after$y * density.scale.shum.after
)
density.by.year.df.shum.after$y <- density.by.year.df.shum.after$y * density.scale.shum.after

# 3. Plot
after.shum.hs.plot <- ggplot(Shumagin.hs.after, aes(x = dummy_date, y = hatch.success)) +
  geom_line(aes(color = year), alpha = 0.2, linewidth = 1.5, show.legend = FALSE) +
  geom_smooth(color = "black", linewidth = 2, se = FALSE) +
  geom_area(data = density.by.year.df.shum.after, aes(x = x, y = y, fill = year), alpha = 0.3, position = "identity", inherit.aes = FALSE) +
  geom_area(data = density.all.df.shum.after, aes(x = x, y = y), fill = "black", alpha = 0.6, position = "identity", inherit.aes = FALSE) +
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ . / density.scale.shum.after)
  ) +
  labs(x = NULL, color = "After", fill = "After", y = NULL) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 26),
    legend.text = element_text(size = 26),
    legend.title = element_text(size = 26, face = "bold"),
    legend.position = "top",
    strip.text = element_text(size = 32),
    plot.margin = margin(t = 15),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_manual(values = year.colors) +
  scale_fill_manual(values = year.colors)




###------ Date range where 95% of data falls------------------------------------
shum.range_95.after <- find.centered.range(density.all.df.shum.after, 0.95)
shum.range_95.after

# April 1 to May 14th - Falls within area of highest hatch success



###---- Determining how many days on average were above 75% after to MHW (Using the geom_smooth line)---------------------
# Fit the model 
fit.shum.after <- loess(hatch.success ~ as.numeric(dummy_date), data = Shumagin.hs.after, span = 0.75)

# Generate a date sequence for the full range
daily.dates.shum.after <- tibble(
  date = seq(min(Shumagin.hs.after$dummy_date), max(Shumagin.hs.after$dummy_date), by = "1 day")
)

# Convert date to numeric for prediction
daily.dates.shum.after <- tibble(
  date = seq(as.Date("2000-01-15"), as.Date("2000-06-15"), by = "1 day")
) %>%
  mutate(dummy_date = as.numeric(date))

daily.dates.shum.after$predicted.hatch.success <- predict(fit.shum.after, newdata = daily.dates.shum.after)

# Count days with hatch success above 75%
sum(daily.dates.shum.after$predicted.hatch.success > 0.75, na.rm = TRUE)

# 153 days above 75%




##------------------Merging plots-----------------------------------------------

# This creates some margins around the plots so when they are spaced out more
before.shum.hs.plot <- before.shum.hs.plot + theme(plot.margin = margin(b = 15, r = 40))
adjacent.shum.hs.plot <- adjacent.shum.hs.plot + theme(plot.margin = margin(l = 39, b= 15))
during.shum.hs.plot <- during.shum.hs.plot + theme(plot.margin = margin(t = 15, r = 15))
after.shum.hs.plot <- after.shum.hs.plot + theme(plot.margin = margin(t = 15, l = 15))


# This combines the four above plots into one plot
combined.hs.shumagin <- (before.shum.hs.plot + adjacent.shum.hs.plot) / (during.shum.hs.plot + after.shum.hs.plot)

#Adding axis labels to the plot
combined.hs.shumagin <- grid.arrange(
  patchworkGrob(combined.hs.shumagin),
  left = textGrob("Hatch Success", rot = 90, gp = gpar(fontsize = 32, fontface = "bold")),
  bottom = textGrob("Day of Year", gp = gpar(fontsize = 32, fontface = "bold"), vjust = 0.3),
  right = textGrob("Hatch Date Density", rot = -90, gp = gpar(fontsize = 32, fontface = "bold"))
)




#-----------------Hatch Success by Depth Plot-----------------------------------


##-----------Before-------------------------------------------------------------
before.hs.depth <- ggplot(before.incub.temps.HS, aes(dummy_date, hatch.success, color = depth)) +
  geom_line(linewidth = 1.5) +
  theme_classic() + 
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_color_manual(values = depth.colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL) +
  ylim(0, 1) +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 28),
    legend.position = "top",
    plot.margin = margin(b = 15, r = 15),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) + 
  annotate("text", x = as.Date("2000-01-15"), y = 0.1, label = "Before", size = 10, fontface = "bold", hjust = 0)




##-----------Adjacent-----------------------------------------------------------
adjacent.hs.depth <- ggplot(adjacent.incub.temps.HS, aes(dummy_date, hatch.success, color = depth)) +
  geom_line(linewidth = 1.5) +
  theme_classic() + 
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_color_manual(values = depth.colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL) +
  ylim(0, 1) +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_blank(),
    legend.text = element_text(size = 28),
    legend.position = "top",
    plot.margin = margin(b = 15, l = 15),
    axis.ticks = element_blank()
  ) + 
  annotate("text", x = as.Date("2000-01-15"), y = 0.1, label = "Adjacent", size = 10, fontface = "bold", hjust = 0)




##--------During----------------------------------------------------------------
during.hs.depth <- ggplot(during.incub.temps.HS, aes(dummy_date, hatch.success, color = depth)) +
  geom_line(linewidth = 1.5) +
  theme_classic() + 
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_color_manual(values = depth.colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL) +
  ylim(0, 1) +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 28),
    legend.position = "top",
    plot.margin = margin(t = 15, r = 15)
  ) + 
  annotate("text", x = as.Date("2000-01-15"), y = 0.1, label = "During", size = 10, fontface = "bold", hjust = 0)




##----------After---------------------------------------------------------------
after.hs.depth <- ggplot(after.incub.temps.HS, aes(dummy_date, hatch.success, color = depth)) +
  geom_line(linewidth = 1.5) +
  theme_classic() + 
  scale_x_date(
    date_labels = "%b-%d",
    date_breaks = "1 month",
    limits = as.Date(c("2000-01-15", "2000-06-15"))
  ) +
  scale_color_manual(values = depth.colors) +
  labs(x = NULL,
       y = NULL,
       color = NULL) +
  ylim(0, 1) +
  theme(
    axis.title = element_text(size = 32, face = "bold"),
    axis.text = element_text(size = 28),
    legend.text = element_text(size = 28),
    legend.position = "top",
    plot.margin = margin(t = 15, l = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) + 
  annotate("text", x = as.Date("2000-01-15"), y = 0.1, label = "After", size = 10, fontface = "bold", hjust = 0)


#Combining the above plots
combined.hs.depth <- (before.hs.depth + adjacent.hs.depth) / (during.hs.depth + after.hs.depth)
combined.hs.depth

#Adding the legend to the top
combined.hs.depth <- combined.hs.depth + plot_layout(guides = "collect") & 
  theme(legend.position = "top")

#Adding axis labels
combined.hs.depth <- grid.arrange(
  patchworkGrob(combined.hs.depth),
  left = textGrob("Hatch Success", rot = 90, gp = gpar(fontsize = 32, fontface = "bold")),
  bottom = textGrob("Day of Year", gp = gpar(fontsize = 32, fontface = "bold"), vjust = 0.3)
) 




#------------Hatch Success through the years------------------------------------ 

# The first 22 rows of data don't use the prior 22 days for the averages since there aren't 22 days before
# Cutting off the data at row 23

GAK.50.115.allyrs.hsplot <- GAK.50.115.allyrs.hs %>%
  slice(-1:-22)


# Create the grouping variable
GAK.50.115.allyrs.hsplot <- GAK.50.115.allyrs.hsplot %>%
  mutate(year_group = case_when(
    year %in% c(2007, 2009, 2010, 2012, 2013, 2014) ~ "Group1",
    year %in% c(2015, 2016, 2017, 2019) ~ "Group2",
    year %in% c(2018, 2020) ~ "Group3",
    year %in% c(2021, 2022, 2023) ~ "Group4",
    TRUE ~ "Group5"
  ))

# Create the base plot
base.year.hs.plot <- ggplot(GAK.50.115.allyrs.hsplot, aes(dummy_date, hatch.success, color = year_group)) + 
  geom_line(linewidth = 2) +
  facet_wrap(~year) +
  scale_color_manual(values = c("Group1" = "darkblue",
                                "Group2" = "darkred",
                                "Group3" = "lightpink",
                                "Group4" = "steelblue1",
                                "Group5" = "black")) +
  xlab("Date") + 
  ylab("E. Kodiak Hatch Success") +
  labs(color = NULL) +
  scale_x_date(
    breaks = as.Date(c("2000-02-01", "2000-04-01", "2000-06-01")), 
    date_labels = "%b-%d", 
    limits = as.Date(c("2000-01-15", "2000-06-15"))) +
  theme_classic() +
  theme(axis.title = element_text(size = 30, face = "bold"),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.position = "none",
        strip.text = element_text(size = 30)) +
  guides(fill = guide_legend(nrow = 1))


# Code to make the year text match the grouping colors

# Convert to grob and modify strip colors
base.year.hs.plot.final <- ggplotGrob(base.year.hs.plot)

# Function to get color for a year
get.year.color <- function(year) {
  if(length(year) == 0 || is.na(year)) return("black")
  if(year %in% c(2007, 2009, 2010, 2012, 2013, 2014)) return("darkblue")
  if(year %in% c(2015, 2016, 2017, 2019)) return("darkred")
  if(year %in% c(2018, 2020)) return("lightpink")
  if(year %in% c(2021, 2022, 2023)) return("steelblue1")
  return("black")
}

# Find strip indices
strip.indices <- which(grepl('strip', base.year.hs.plot.final$layout$name))

# Modify each strip
for(i in strip.indices) {
  tryCatch({
    # Navigate through the grob structure to find the text label
    strip_grob <- base.year.hs.plot.final$grobs[[i]]
    
    # Try to extract text from different possible locations
    strip_text <- NULL
    if(inherits(strip_grob$grobs[[1]]$children[[2]]$children[[1]], "text")) {
      strip_text <- strip_grob$grobs[[1]]$children[[2]]$children[[1]]$label
    }
    
    if(!is.null(strip_text) && length(strip_text) > 0) {
      year_value <- as.numeric(strip_text)
      
      if(!is.na(year_value)) {
        # Set the color based on the year
        base.year.hs.plot.final$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- get.year.color(year_value)
      }
    }
  }, error = function(e) {
    message(paste("Error processing strip", i, ":", e$message))
  })
}

# Draw the plot
grid.newpage()
grid.draw(base.year.hs.plot.final)




##-------------------- Number of days from 99 to 06 above 75% ------------------


# Count days with hatch success above 75%
GAK.50.115.allyrs.hsplot %>%
  mutate(md = format(date, "%m-%d")) %>%
  filter(md >= "01-15", md <= "06-15") %>%
  group_by(year) %>%
  count(hatch.success > 0.75)




