# Load packages
library(tidyverse)
library(janitor)
library(readxl)
library(lme4)
library(lmerTest)
library(ggpubr)
library(car)

# Load and clean data
response <- read.csv("mill_mesocosm_full_data4.csv") %>% clean_names()
meta <- read.csv("mesocosm_meta.csv") %>% clean_names()

# Join and prepare data
meso_join <- response %>% 
  full_join(meta, by = "seaweed_id") %>%
  group_by(seaweed_id, date_sampled, sp) %>%
  mutate(total_colony_area = sum(colony_area, na.rm = TRUE)) %>%
  ungroup()

meso_adult <- meso_join %>%
  filter(!is.na(turbulence), turbulence != "", sp %in% c("mm", "ep"),
         grepl("^A", seaweed_id), !(notes %in% c("mia", "MIA", "missing"))) %>%
  mutate(colony_area = ifelse(is.na(colony_area), 0, colony_area))


# Adult colony area plot
ggplot(meso_adult, aes(x = turbulence, y = colony_area, colour = sp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.75)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  facet_grid(date_sampled ~ nutrient_treatment) +
  theme_minimal(base_size = 13) +
  labs(title = "Adult:Colony Area", y = "Colony Area")


# Adult analysis (non-YO)
anova_data <- meso_join %>%
  filter(!is.na(colony_area),
         sp %in% c("ep", "mm"),
         !is.na(turbulence), !is.na(nutrient_treatment), !is.na(date_sampled),
         !grepl("^A", seaweed_id)) %>%
  group_by(seaweed_id, sp, turbulence, nutrient_treatment, date_sampled) %>%
  summarise(total_colony_area = sum(colony_area, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(c(turbulence, nutrient_treatment, date_sampled), as.factor))

# LMM for adults
adult_model <- lmer(
  total_colony_area ~ sp * turbulence * nutrient_treatment + date_sampled +
    (1 | main_experiment_tank),
  data = anova_data
)
summary(adult_model)
Anova(adult_model, type = 3)

# Adult plot
ggplot(anova_data, aes(x = turbulence, y = total_colony_area, colour = sp)) +
  geom_boxplot() +
  facet_grid(date_sampled ~ nutrient_treatment) +
  theme_minimal() +
  labs(title = "Adult: Total Colony Area", y = "Total Colony Area")






# Young analysis (YO only)
meso_young <- meso_join %>%
  filter(!is.na(turbulence), turbulence != "", sp %in% c("mm", "ep"),
         grepl("^YO", seaweed_id), !(notes %in% c("mia", "MIA", "missing"))) %>%
  mutate(colony_area = ifelse(is.na(colony_area), 0, colony_area)) %>%
  group_by(seaweed_id, blade_no, sp, turbulence, nutrient_treatment,
           date_sampled, main_experiment_tank) %>%
  summarise(
    blade_colony_avg = mean(colony_area, na.rm = TRUE),
    blade_colony_count = sum(colony_area > 0, na.rm = TRUE),
    blade_length = mean(blade_length, na.rm = TRUE),
    .groups = "drop"
  )

# LMM for young (avg area)
lmm_young_area <- lmer(
  blade_colony_avg ~ sp * turbulence + nutrient_treatment + date_sampled +
    (1 | main_experiment_tank) + (1 | seaweed_id),
  data = meso_young
)
summary(lmm_young_area)
Anova(lmm_young_area, type = 3)

# GLMM for colony count (Poisson)
glmm_young_count <- glmer(
  blade_colony_count ~ sp * turbulence + nutrient_treatment + date_sampled +
    (1 | main_experiment_tank) + (1 | seaweed_id),
  data = meso_young,
  family = poisson(link = "log")
)
summary(glmm_young_count)

# Young colony area plot
ggplot(meso_young, aes(x = turbulence, y = blade_colony_avg, colour = sp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.75)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  facet_grid(date_sampled ~ nutrient_treatment) +
  theme_minimal(base_size = 13) +
  labs(title = "Young: Mean Colony Area per Blade", y = "Avg Colony Area")

# Young colony count plot
ggplot(meso_young, aes(x = turbulence, y = blade_colony_count, colour = sp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.75)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  facet_grid(date_sampled ~ nutrient_treatment) +
  theme_minimal(base_size = 13) +
  labs(title = "Young: Colony Count per Blade", y = "Colony Count")

