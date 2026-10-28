library(tidyverse)
library(glue)
library(ggpattern)
library(ggpubr)
library(readxl) 
library(janitor) 
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(scales)

rm(list = ls())

# Set working directory
setwd("/biotransformation_k_R/")

colors_plot <- c("up1" = "chartreuse2",
                 "up2" = "seagreen4",  
                 "dw1" = "navy",
                 "dw2" = "deepskyblue4", 
                 "dw3" = "mediumturquoise")

#################################################################################
#river data
data_air_in <- read_xlsx("data/HOBO/HOBO_air.xlsx") %>% 
  clean_names() %>%  
  mutate(datetime = {start <- as.POSIXct("04.20.22 10:00:00 AM",
                          format = "%m.%d.%y %I:%M:%S %p",
                          tz = "UTC")  
      seq.POSIXt(from = start, by = "10 min", length.out = n())},
      date_only = format(datetime, "%d/%m/%Y"),
      time_only = format(datetime, "%H:%M")) %>%
  select(!date) %>%
  filter(datetime <= as.POSIXct("2022-05-09 23:50:00", tz = "UTC"))

data_oli_in <- read_xlsx("data/HOBO/HOBO_oli.xlsx") %>% 
  clean_names() %>%  
  mutate(datetime = {start <- as.POSIXct("04.20.22 12:00:00 PM",
                                         format = "%m.%d.%y %I:%M:%S %p",
                                         tz = "UTC")  
  seq.POSIXt(from = start, by = "10 min", length.out = n())},
  date_only = format(datetime, "%d/%m/%Y"),
  time_only = format(datetime, "%H:%M")) %>%
  select(!date) %>%
  filter(datetime <= as.POSIXct("2022-05-09 23:50:00", tz = "UTC"))
  
data_bir_in <- read_xlsx("data/HOBO/HOBO_bir.xlsx") %>% 
  clean_names() %>%  	
  mutate(datetime = {start <- as.POSIXct("06.13.22 10:00:00 PM",
                                         format = "%m.%d.%y %I:%M:%S %p",
                                         tz = "UTC")  
  seq.POSIXt(from = start, by = "5 min", length.out = n())},
  date_only = format(datetime, "%d/%m/%Y"),
  time_only = format(datetime, "%H:%M")) %>%
  select(!date)
  
data_kol_in <- read_xlsx("data/HOBO/HOBO_kol.xlsx") %>% 
  clean_names() %>%  	
  mutate(datetime = {start <- as.POSIXct("06.14.22 10:00:00 AM",
                                         format = "%m.%d.%y %I:%M:%S %p",
                                         tz = "UTC")  
  seq.POSIXt(from = start, by = "5 min", length.out = n())},
  date_only = format(datetime, "%d/%m/%Y"),
  time_only = format(datetime, "%H:%M")) %>%
  select(!date) %>%
  filter(datetime >= as.POSIXct("2022-06-14 12:00:00", tz = "UTC"))

data_alt_in <- read_xlsx("data/HOBO/HOBO_alt.xlsx") %>% 
  clean_names() %>%  
  mutate(datetime = {start <- as.POSIXct("09.02.22 06:00:00 AM",
                                         format = "%m.%d.%y %I:%M:%S %p",
                                         tz = "UTC")  
  seq.POSIXt(from = start, by = "5 min", length.out = n())},
  date_only = format(datetime, "%d/%m/%Y"),
  time_only = format(datetime, "%H:%M")) %>%
  select(!date)

data_ehr_in <- read_xlsx("data/HOBO/HOBO_ehr.xlsx") %>% 
  clean_names() %>%  	
  mutate(datetime = {start <- as.POSIXct("09.02.22 06:00:00 AM",
                                         format = "%m.%d.%y %I:%M:%S %p",
                                         tz = "UTC")  
  seq.POSIXt(from = start, by = "5 min", length.out = n())},
  date_only = format(datetime, "%d/%m/%Y"),
  time_only = format(datetime, "%H:%M")) %>%
  select(!date)
#############################################################################
#summarize
df_air_sum <- data_air_in %>%
  select(ends_with("temp")) %>%
  pivot_longer(cols = everything(),names_to = "site",values_to = "temp") %>%
  mutate(site = sub("_temp$", "", site)) %>%
  group_by(site) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), .groups = "drop")

df_oli_sum <- data_oli_in %>%
  select(ends_with("temp")) %>%
  pivot_longer(cols = everything(),names_to = "site",values_to = "temp") %>%
  mutate(site = sub("_temp$", "", site)) %>%
  group_by(site) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), .groups = "drop") %>%
  bind_rows(tibble(site = "oli_dw3", mean_temp = NA_real_, sd_temp = NA_real_))
 
df_bir_sum <- data_bir_in %>%
  select(ends_with("temp")) %>%
  pivot_longer(cols = everything(),names_to = "site",values_to = "temp") %>%
  mutate(site = sub("_temp$", "", site)) %>%
  group_by(site) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), .groups = "drop")

df_kol_sum <- data_kol_in %>%
  select(ends_with("temp")) %>%
  pivot_longer(cols = everything(),names_to = "site",values_to = "temp") %>%
  mutate(site = sub("_temp$", "", site)) %>%
  group_by(site) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), .groups = "drop") %>%
  bind_rows(tibble(site = "kol_up2", mean_temp = NA_real_, sd_temp = NA_real_))

df_alt_sum <- data_alt_in %>%
  select(ends_with("temp")) %>%
  pivot_longer(cols = everything(),names_to = "site",values_to = "temp") %>%
  mutate(site = sub("_temp$", "", site)) %>%
  group_by(site) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), .groups = "drop")

df_ehr_sum <- data_ehr_in %>%
  select(ends_with("temp")) %>%
  pivot_longer(cols = everything(),names_to = "site",values_to = "temp") %>%
  mutate(site = sub("_temp$", "", site)) %>%
  group_by(site) %>%
  summarise(mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), .groups = "drop")

tot_temp_sum <- bind_rows(df_air_sum, df_alt_sum, df_bir_sum, df_ehr_sum, df_kol_sum, df_oli_sum) %>%
  mutate(across(c(mean_temp, sd_temp), ~ round(.x, 2)))

write.table(tot_temp_sum, "data/riverT_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#############################################################################
#temperature
air_temp_long <- data_air_in %>%
  select(!c(ends_with("_lux"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^air_") %>% str_remove("_temp$"), starts_with("air")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

air_temp <- ggplot(air_temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(air_temp_long$datetime), to = max(air_temp_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title    = element_text(size = 16, hjust = 0),
        legend.position = c(0.9, 0.8),                      
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA)) +
  labs(x = "Date", y = "Temperature [°C]", color = "Site", title = "Airolo",)

air_temp

ggsave("plots/field data/air_temp_HOBO.png", plot = air_temp , width = 10, height = 6, dpi = 150)

oli_temp_long <- data_oli_in %>%
  select(!c(ends_with("_lux"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^oli_") %>% str_remove("_temp$"), starts_with("oli")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

oli_temp <- ggplot(oli_temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(oli_temp_long$datetime), to = max(oli_temp_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title    = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Temperature [°C]", title = "Olivone")

oli_temp

ggsave("plots/field data/oli_temp_HOBO.png", plot = oli_temp , width = 10, height = 6, dpi = 150)

bir_temp_long <- data_bir_in %>%
  select(!c(ends_with("_lux"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^bir_") %>% str_remove("_temp$"), starts_with("bir")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

bir_temp <- ggplot(bir_temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(bir_temp_long$datetime), 
                                to = max(bir_temp_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Temperature [°C]", title = "Birmensdorf")

bir_temp

ggsave("plots/field data/bir_temp_HOBO.png", plot = bir_temp , width = 10, height = 6, dpi = 150)

kol_temp_long <- data_kol_in %>%
  select(!c(ends_with("_lux"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^kol_") %>% str_remove("_temp$"), starts_with("kol")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

kol_temp <- ggplot(kol_temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(kol_temp_long$datetime), 
                                to = max(kol_temp_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Temperature [°C]", title = "Kölliken")

kol_temp

ggsave("plots/field data/kol_temp_HOBO.png", plot = kol_temp , width = 10, height = 6, dpi = 150)

alt_temp_long <- data_alt_in %>%
  select(!c(ends_with("_lux"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^alt_") %>% str_remove("_temp$"), starts_with("alt")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

alt_temp <- ggplot(alt_temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(alt_temp_long$datetime), 
                                to = max(alt_temp_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Temperature [°C]", title = "Altstätten")

alt_temp

ggsave("plots/field data/alt_temp_HOBO.png", plot = alt_temp , width = 10, height = 6, dpi = 150)

ehr_temp_long <- data_ehr_in %>%
  select(!c(ends_with("_lux"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^ehr_") %>% str_remove("_temp$"), starts_with("ehr")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

ehr_temp <- ggplot(ehr_temp_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(ehr_temp_long$datetime), 
                                to = max(ehr_temp_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 25),
                     breaks = seq(0, 25, by = 5)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Temperature [°C]", title = "Ehrendingen")

ehr_temp

ggsave("plots/field data/ehr_temp_HOBO.png", plot = ehr_temp , width = 10, height = 6, dpi = 150)

#############################################################################
#light
air_lux_long <- data_air_in %>%
  select(!c(ends_with("_temp"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^air_") %>% str_remove("_lux$"), starts_with("air")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "illuminance") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

air_lux <- ggplot(air_lux_long, aes(x = datetime, y = illuminance, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(air_lux_long$datetime), to = max(air_lux_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 300000),
                     breaks = seq(0, 300000, by = 50000),
                     labels = label_scientific()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = c(0.9, 0.8),                      
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA)) +
  labs(x = "Date", y = "Illuminance [lx]", color = "Site", title = "Airolo",)

air_lux

ggsave("plots/field data/air_lux_HOBO.png", plot = air_lux , width = 10, height = 6, dpi = 150)

oli_lux_long <- data_oli_in %>%
  select(!c(ends_with("_temp"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^oli_") %>% str_remove("_lux$"), starts_with("oli")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "illuminance") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

oli_lux <- ggplot(oli_lux_long, aes(x = datetime, y = illuminance, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(oli_lux_long$datetime), to = max(oli_lux_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 300000),
                     breaks = seq(0, 300000, by = 50000),
                     labels = label_scientific()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Illuminance [lx]", title = "Olivone")

oli_lux

ggsave("plots/field data/oli_lux_HOBO.png", plot = oli_lux , width = 10, height = 6, dpi = 150)

bir_lux_long <- data_bir_in %>%
  select(!c(ends_with("_temp"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^bir_") %>% str_remove("_lux$"), starts_with("bir")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "temperature") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

bir_lux <- ggplot(bir_lux_long, aes(x = datetime, y = temperature, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(bir_lux_long$datetime), 
                                to = max(bir_lux_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 300000),
                     breaks = seq(0, 300000, by = 50000),
                     labels = label_scientific()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Illuminance [lx]", title = "Birmensdorf")

bir_lux

ggsave("plots/field data/bir_lux_HOBO.png", plot = bir_lux, width = 10, height = 6, dpi = 150)

kol_lux_long <- data_kol_in %>%
  select(!c(ends_with("_temp"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^kol_") %>% str_remove("_lux$"), starts_with("kol")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "illuminance") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

kol_lux <- ggplot(kol_lux_long, aes(x = datetime, y = illuminance, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(kol_lux_long$datetime), 
                                to = max(kol_lux_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 300000),
                     breaks = seq(0, 300000, by = 50000),
                     labels = label_scientific()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Illuminance [lx]", title = "Kölliken")

kol_lux

ggsave("plots/field data/kol_lux_HOBO.png", plot = kol_lux , width = 10, height = 6, dpi = 150)

alt_lux_long <- data_alt_in %>%
  select(!c(ends_with("_temp"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^alt_") %>% str_remove("_lux$"), starts_with("alt")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "illuminance") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

alt_lux <- ggplot(alt_lux_long, aes(x = datetime, y = illuminance, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(alt_lux_long$datetime), 
                                to = max(alt_lux_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 300000),
                     breaks = seq(0, 300000, by = 50000),
                     labels = label_scientific()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Illuminance [lx]", title = "Altstätten")

alt_lux

ggsave("plots/field data/alt_lux_HOBO.png", plot = alt_lux , width = 10, height = 6, dpi = 150)

ehr_lux_long <- data_ehr_in %>%
  select(!c(ends_with("_temp"), date_only, time_only)) %>%
  rename_with(~ str_remove(., "^ehr_") %>% str_remove("_lux$"), starts_with("ehr")) %>%
  pivot_longer(cols = -datetime, names_to = "site", values_to = "illuminance") %>%
  mutate(datetime = as.POSIXct(datetime),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")))

ehr_lux <- ggplot(ehr_lux_long, aes(x = datetime, y = illuminance, color = site)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_plot) +
  scale_x_datetime(breaks = seq(from = min(ehr_lux_long$datetime), 
                                to = max(ehr_lux_long$datetime), length.out = 6),
                   labels = date_format("%d-%b", locale = "en_US"),
                   timezone = "UTC") +
  scale_y_continuous(limits = c(0, 300000),
                     breaks = seq(0, 300000, by = 50000),
                     labels = label_scientific()) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0),
        legend.position = "none") +
  labs(x = "Date", y = "Illuminance [lx]", title = "Ehrendingen")

ehr_lux

ggsave("plots/field data/ehr_lux_HOBO.png", plot = ehr_lux, width = 10, height = 6, dpi = 150)
