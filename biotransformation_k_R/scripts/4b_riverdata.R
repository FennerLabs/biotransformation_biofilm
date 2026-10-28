library(tidyverse)
library(glue)
library(ggpattern)
library(ggpubr)
library(readxl) 
library(janitor) 
library(RColorBrewer)
library(patchwork)

rm(list = ls())

# Set working directory
setwd("/biotransformation_k_R/")

colors_plot <- c("up1" = "chartreuse2",
                 "up2" = "seagreen4",  
                 "dw1" = "navy",
                 "dw2" = "deepskyblue4", 
                 "dw3" = "mediumturquoise")

loc_colors <- c("air" = "blue1", 
                "oli" = "blue4", 
                "bir" = "seagreen3", 
                "kol" = "seagreen4", 
                "alt" = "orange", 
                "ehr" = "orange3")

#################################################################################
#river data
data_in <- read_xlsx("data/riverParameter_summary.xlsx") %>% 
  clean_names() %>% 
  mutate(parameter = paste(parameter, " [", einheit, "]", sep = "")) %>%
  select(!c(angabe, einheit)) %>%
  t() %>% 
  as.data.frame() %>% 
  {colnames(.) <- .[1, ]   # Set the first row as column names
  .[-1, ]     } %>%  # Remove the first row
  rownames_to_column(var = "camp_loc_site") %>%
  filter(camp_loc_site != "loq_value") %>%
  mutate(across(!c(camp_loc_site), ~ round(as.numeric(.), 3))) %>%
  mutate(campaign = word(camp_loc_site, 1, sep = "\\_"),
         location = word(camp_loc_site, 2, sep = "\\_"),
         site = word(camp_loc_site, 3, sep = "\\_")) %>%
  mutate(site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")),
         campaign = factor(campaign, levels = c("fie1", "fie2", "fie3")),
         location = factor(location, levels = c("air", "oli", "bir", "kol", "alt", "ehr")))
  
colnames(data_in) <- c("camp_loc_site", "Conductivity [µS/cm 20°C]", "pH", "Alkalinity [mmol/L]", "Total Hardness [mmol/L]", "Sodium [mg/L]", 
                       "Magnesium [mg/L]", "Calcium [mg/L]", "Potassium [mg/L]", "Fluoride [mg/L]", "Chloride [mg/L]", "Bromide [mg/L]",
                       "Nitrate [mg/L]", "Sulfate [mg/L]", "Ammonium [µg/L]", "Nitrite [µg/L]", "Phosphate [µg/L]", "TP [µg/L]",
                       "Silicic Acid [mg/L]", "DOC [mg/L]", "TOC [mg/L]", "TN [mg/L]", "TIC [mg/L]", "campaign", "location", "site")

data_long <- data_in %>%
  pivot_longer(cols = c("Conductivity [µS/cm 20°C]", "pH", "Alkalinity [mmol/L]", "Total Hardness [mmol/L]", "Sodium [mg/L]", 
             "Magnesium [mg/L]", "Calcium [mg/L]", "Potassium [mg/L]", "Fluoride [mg/L]", "Chloride [mg/L]", "Bromide [mg/L]",
             "Nitrate [mg/L]", "Sulfate [mg/L]", "Ammonium [µg/L]", "Nitrite [µg/L]", "Phosphate [µg/L]", "TP [µg/L]",
             "Silicic Acid [mg/L]", "DOC [mg/L]", "TOC [mg/L]", "TN [mg/L]", "TIC [mg/L]"),
             names_to = "parameter",
             values_to = "value") %>% 
  group_by(parameter) %>%
  mutate(maxvalue = max(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(parameter = as.factor(parameter))

ggplot(data_long, aes(x = site, y = value, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = colors_plot) +
  facet_grid(parameter ~ location, scales = "free_y", switch = "y") + 
  labs(y = "river parameter") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),  
    axis.ticks.y = element_line(), 
    strip.text.y.left = element_text(angle = 0, hjust = 0.5),  
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    strip.placement = "outside",  
    strip.text.x = element_text(size = 10))  

ggsave("plots/riverdata.png", plot = last_plot(), width = 8, height = 12, units = "in", dpi = 150)
