# Libraries
library(dplyr)
library(ggplot2)
library(lubridate)
library(viridis)

# ==== 1) Create data frame ====
df <- read.csv("../seaweed_area_kc.csv") %>% 
  clean_names()

# ==== 2) Clean ====
df <- df %>%
  mutate(
    Date = dmy(date),
    Month = month(Date, label = TRUE, abbr = TRUE),
    for_standardizing = as.numeric(for_standardizing)
  ) %>%
  filter(!is.na(for_standardizing)) %>%
  filter(seaweed_species != "Digitata") %>%
  group_by(seaweed_species, Month) %>%
  summarise(mean_blade_area = mean(for_standardizing, na.rm = TRUE), .groups = "drop")

# ==== 3) Plot ====
ggplot(df, aes(x = Month, y = mean_blade_area, fill = seaweed_species)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d(option = "plasma", end = 0.85) +
  labs(
    x = NULL,
    y = expression("Mean blade area (cm"^2*")"),
    fill = "Seaweed species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )
