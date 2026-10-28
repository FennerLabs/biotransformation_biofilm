library(tidyverse)
library(glue)
library(ggpattern)
library(ggpubr)
library(readxl) 
library(janitor) 
library(RColorBrewer)
library(rstatix)

rm(list = ls())

# Set working directory
setwd("/biotransformation_k_R/")

colors_plot <- c("up1_r1" = "chartreuse2", "up1_r2" = "#4CAF3C", 
                 "up2_r1" = "seagreen4", "up2_r2" = "#20543A", 
                 "dw1_r1" = "navy","dw1_r2" = "#2B3E8F", 
                 "dw2_r1" = "#005073", "dw2_r2" = "deepskyblue4", 
                 "dw3_r1" = "#008078", "dw3_r2" = "mediumturquoise")

loc_colors <- c("air" = "blue1", 
                "oli" = "blue4", 
                "bir" = "seagreen3", 
                "kol" = "seagreen4", 
                "alt" = "orange", 
                "ehr" = "orange3")


data_in <- read_xlsx("data/biodata_allcampaigns.xlsx") %>% 
  clean_names() %>%
  mutate(location = word(loc_batch, 1, sep = "\\_"),
         site = word(loc_batch, 2, sep = "\\_"),
         replicate = word(loc_batch, 3, sep = "\\_"),
         batch = str_c(word(loc_batch, 2, sep = "\\_"), word(loc_batch, 3, sep = "\\_"), sep = "_")) %>%
  filter(!(site == "asb")) %>%
  mutate(site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")),
         campaign = factor(campaign, levels = c("fie1", "fie2", "fie3")),
         location = factor(location, levels = c("air", "oli", "bir", "kol", "alt", "ehr")),
         replicate = factor(replicate, levels = c("r1", "r2")), 
         batch = factor(batch, levels = c("up1_r1", "up1_r2", 
                                     "up2_r1", "up2_r2", 
                                     "dw1_r1","dw1_r2", 
                                     "dw2_r1", "dw2_r2", 
                                     "dw3_r1", "dw3_r2")))

#################################################################################
# VSS
data_vss <- data_in %>%
  select(-c(tss_g_l_0, od_1, icc_n_l_1_0, tss_g_l_48, icc_n_l_1_48, tss_g_l_96, icc_n_l_1_96,  chl_a_mg_m_l_1)) %>%
  mutate(across(c(afdm_g_l_0, afdm_g_l_48, afdm_g_l_96), ~ round(as.numeric(.), 3))) %>%
  pivot_longer(cols = starts_with("afdm_g_l_"), names_to = "time_point", values_to = "afdm") %>%
  mutate(time_point = factor(word(time_point, 4, sep = "\\_"), levels = c(0, 48, 96))) %>%
  rename(VSS = afdm)

ggplot(data_vss, aes(x = batch, y = VSS, fill = batch, pattern = time_point)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(width = 0.9), 
                   color = "black",
                   pattern_fill = "black") +
  facet_wrap(~location, scales = "free", ncol = 2) +
  scale_fill_manual(values = colors_plot) +
  scale_pattern_manual(values = c("0" = "none", "48" = "circle", "96" = "crosshatch")) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), limits = c(0, 1.5)) +
  theme_minimal() +
  labs(x = "Reactor",
       y = "VSS [g/L]",
       fill = "Reactor",
       pattern = "Time Point") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25))

ggsave("plots/biodata/VSS.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 150)
#################################################################################
# TSS
data_tss <- data_in %>%
  select(-c(afdm_g_l_0, od_1, icc_n_l_1_0, afdm_g_l_48, icc_n_l_1_48, afdm_g_l_96, icc_n_l_1_96,  chl_a_mg_m_l_1)) %>%
  mutate(across(c(tss_g_l_0, tss_g_l_48, tss_g_l_96), ~ round(as.numeric(.), 3))) %>%
  pivot_longer(cols = starts_with("tss_g_l_"), names_to = "time_point", values_to = "TSS") %>%
  mutate(time_point = factor(word(time_point, 4, sep = "\\_"), levels = c(0, 48, 96)))

ggplot(data_tss, aes(x = batch, y = TSS, fill = batch, pattern = time_point)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(width = 0.9), 
                   color = "black",
                   pattern_fill = "black") +
  facet_wrap(~location, scales = "free", ncol = 2) +
  scale_fill_manual(values = colors_plot) +
  scale_pattern_manual(values = c("0" = "none", "48" = "circle", "96" = "crosshatch")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6), limits = c(0, 6)) +
  theme_minimal() +
  labs(x = "Reactor",
       y = "TSS [g/L]",
       fill = "Reactor",
       pattern = "Time Point") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25))

ggsave("plots/biodata/TSS.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 150)

#################################################################################
#FC
data_fc <- data_in %>%
  select(-c(afdm_g_l_0, od_1, tss_g_l_0, afdm_g_l_48, tss_g_l_48, afdm_g_l_96, tss_g_l_96,  chl_a_mg_m_l_1)) %>%
  mutate(across(c(icc_n_l_1_0, icc_n_l_1_48, icc_n_l_1_96), ~ round(as.numeric(.), 0))) %>%
  pivot_longer(cols = starts_with("icc_n_l_1_"), names_to = "time_point", values_to = "icc_nl") %>%
  mutate(time_point = factor(word(time_point, 5, sep = "\\_"), levels = c(0, 48, 96)))

ggplot(data_fc, aes(x = batch, y = icc_nl, fill = batch, pattern = time_point)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(width = 0.9), 
                   color = "black",
                   pattern_fill = "black") +
  facet_wrap(~location, scales = "free", ncol = 2) +
  scale_fill_manual(values = colors_plot) +
  scale_pattern_manual(values = c("0" = "none", "48" = "circle", "96" = "crosshatch")) +
  scale_y_continuous(breaks = c(0, 100, 150, 200, 250), limits = c(0, 300)) +
  theme_minimal() +
  labs(x = "Batch",
       y = "ICC [cells/nl]",
       fill = "Batch",
       pattern = "Time Point") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25))

ggsave("plots/biodata/FC_absolute.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 150)

data_fc_rel <- data_fc %>%
  select(batch, icc_nl, location, time_point) %>%
  group_by(batch, location) %>%
  mutate(icc_nl_rel = icc_nl / icc_nl[time_point == 0])

ggplot(data_fc_rel, aes(x = batch, y = icc_nl_rel, fill = batch, pattern = time_point)) +
  geom_bar_pattern(stat = "identity", position = position_dodge(width = 0.9), 
                   color = "black",
                   pattern_fill = "black") +
  facet_wrap(~location, scales = "free", ncol = 2) +
  scale_fill_manual(values = colors_plot) +
  scale_pattern_manual(values = c("0" = "none", "48" = "circle", "96" = "crosshatch")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 120)) +
  theme_minimal() +
  labs(x = "Batch",
       y = "ICC growth [cells/nl]",
       fill = "Batch",
       pattern = "Time Point") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25))

ggsave("plots/biodata/FC_relative.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 150)
####################################################################################
#Chl_a
data_chla <- data_in %>%
  select(c(campaign, loc_batch, location, site, replicate, batch, chl_a_mg_m_l_1)) %>%
  mutate(chla = round(1000*chl_a_mg_m_l_1, 1)) 
  
ggplot(data_chla, aes(x = batch, y = chla, fill = batch)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~location, scales = "free", ncol = 2) +
  scale_fill_manual(values = colors_plot) +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50), limits = c(0, 55)) +
  theme_minimal() +
  labs(x = "Batch",
       y = "Chl_a [µg/mL]",
       fill = "Batch") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25))

ggsave("plots/biodata/Chla.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 150)

####################################################################################
#OD
data_od <- data_in %>%
  select(c(campaign, loc_batch, location, site, replicate, batch, od_1)) %>%
  mutate(od_1 = ifelse(location %in% c("air", "oli"), od_1*2, od_1)) %>%
  mutate(od_1 = round(od_1, 2))

ggplot(data_od, aes(x = batch, y = od_1, fill = batch)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~location, scales = "free", ncol = 2) +
  scale_fill_manual(values = colors_plot) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5), limits = c(0, 2.6)) +
  theme_minimal() +
  labs(x = "Batch",
       y = "OD685",
       fill = "Batch") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5),
        panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25))

ggsave("plots/biodata/OD685.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 150)
