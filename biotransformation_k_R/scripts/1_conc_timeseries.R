library(tidyverse)
library(readxl)
library(janitor)

rm(list = ls())

# Set working directory
setwd("/biotransformation_k_R/")

sub_out <- c("Allopurinol", "Asulam", "Benserazide", "Clotrimazole", "Canagliflozin", "Dasatinib", "Dolutegravir", "Gemfibrozil",
            "Hydroxychloroquine", "Ibuprofen", "Lovastatin", "Mesalazine", "Metformin", "Neohesperidin dihydrochalcone", "Nintedanib",
            "Olanzapine", "Piracetam", "Sertraline", "Simvastatin", "Spironolactone", "Terbinafine", "Regorafenib", "Acemetacin", 
            "Aliskiren", "Apremilast", "Atorvastatin", "Bisacodyl", "Dexlansoprazole", "Diltiazem", "Edoxaban", "Irbesartan", "Indomethacin", 
            "Kresoxim-methyl", "Neotame", "Omeprazole", "Tolperisone", "Valaciclovir", "Amlodipine", "Atovaquone", "Bupropion", 
            "Clomipramine", "Clozapine", "Diosmin", "Duloxetine", "Efavirenz", "Elvitegravir", "Entacapone", "Fluoxetine", "Ticlopidine", 
            "Trimipramine","Apixaban", "Atomoxetine", "Climbazole", "Cyromazine", "Dextromethorphan", "Diphenhydramine", "Enzalutamide", 
            "Fenfluramine", "Fipronil", "Flonicamid", "Flufenamic acid", "Irgarol", "Phenytoin", "Propranolol", "Tebufenozide", "Sucralose")

#########################################################################
# (1) Read .csv file
# Add clusternumber to table
data_in <- read_csv("data/BF_tot_conc.csv") %>% 
  clean_names()  
  
data <- data_in %>%
  filter(type == "bt") %>%          #ct_series only for biotransformation (bt) experiments
  mutate(location = factor(location, levels = c("air", "oli", "bir", "kol", "alt", "ehr"))) %>%
  group_by(substance, loc_batch) %>%
  mutate(conc_max = max(conc_nmol, na.rm = TRUE)) %>%
  mutate(conc_norm = conc_nmol / conc_max) %>%
  ungroup() %>%
  group_by(substance) %>%
  complete(location, time_point) %>%
  fill(x4lettercode, loc_batch, site, replicate, .direction = "downup") %>%
  ungroup() %>%
  select(substance, x4lettercode, location, time_point, replicate, site, loc_batch, conc_nmol, conc_norm)

#pull substances of interest
substance_selection <- data %>%
  filter(!substance %in% sub_out) %>%
  group_by(substance)%>%
  slice(1)%>%
  ungroup()%>%
  pull(substance)

################################################################################
# (2) Concentration-time Series Plot

colors_plot <- c("up1" = "chartreuse2", "up2" = "seagreen4", 
                 "dw1" = "navy", "dw2" = "deepskyblue4", "dw3" = "mediumturquoise")

for (i in 1:length(substance_selection)) {
  cmp <- data %>%
    filter(substance == substance_selection[i])
  
  empty_locations <- cmp %>%
    group_by(location) %>%
    summarise(is_empty = all(is.na(conc_norm))) %>%
    filter(is_empty) %>%
    pull(location)
  
  label_data <- cmp %>%
    filter(location %in% empty_locations) %>%
    distinct(location) %>%
    mutate(time_point = 48, conc_norm  = 0.7, label = "no data available")
  
  title_text <- cmp %>%
    distinct(substance, x4lettercode) %>%
    mutate(text = paste(substance, "-", x4lettercode)) %>%
    pull(text)
    
  cmp_conc_plot <- ggplot(cmp, aes(x = time_point, y = conc_norm)) +
    geom_line(aes(color = site, linetype = replicate)) +
    geom_point(aes(color = site, shape = replicate)) +
    geom_text(data = label_data, aes(label = label), color = "grey30", size = 4.5, fontface = "italic") +
    labs(x = "Time [h]", y = "Concentration normated [nM]", color = "") +
    scale_x_continuous(breaks = c(0, 24, 48, 72, 96)) +
    scale_y_continuous(breaks = seq(0, 1.2, 0.2), limits = c(0, 1.3)) +
    facet_wrap(~location, nrow = 1, ncol = 6, drop = FALSE) +
    theme_minimal() +
    scale_color_manual(values = colors_plot) +
    ggtitle(title_text)
  
  ggsave(
    paste0("plots/ct_series/", substance_selection[i], ".jpg"),
    cmp_conc_plot,
    device = "jpeg",
    dpi = 300,
    units = "px",
    width = 8000,
    height = 1200
  )
}