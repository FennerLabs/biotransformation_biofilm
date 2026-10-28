library(tidyverse)        
library(glue)
library(dendextend)
library(ggpubr)
library(ComplexHeatmap)
library(readxl)
library(janitor)
library(RColorBrewer)
library(rstatix)
library(circlize)
library(grid)

rm(list = ls())

# Set working directory
setwd("/public/")

################################################################################
# Read and add PassiveSampling table 
df_PS <- read_csv("data/Result_influent_effluent_passiveSampling_complete.csv") %>% 
  clean_names() %>%
  rename(x4lettercode = lettercode4, 
         substance = name) %>%
  filter(type == "PassiveSampling") %>% 
  mutate(location = word(sample, 2, sep = "_"),
         site = word(sample, 3, sep = "_"),
         loc_site = paste(location, site, sep = "_")) %>%
  group_by(substance, sample)%>%
  mutate(meanPS_ng_l = mean(c(ps_extract_ng_l_r1, ps_extract_ng_l_r2), na.rm = TRUE),
         sdPS_ng_l = sd(c(ps_extract_ng_l_r1, ps_extract_ng_l_r2), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(campaign = case_when(location %in% c("air", "oli") ~ "fie1",
                              location %in% c("bir", "kol") ~ "fie2",
                              TRUE ~ "fie3")) %>%
  select(uchem_id, x4lettercode, meanPS_ng_l, loc_site, sdPS_ng_l, campaign, location, site) %>%
  mutate(across(c(meanPS_ng_l, sdPS_ng_l), ~ ifelse(. == 0, NA, .)))

df_PS_blanks <- df_PS %>%
  distinct(uchem_id, x4lettercode) %>%           
  crossing(loc_site = c("air_blank", "oli_blank", "bir_blank", "kol_blank", "ehr_blank", "alt_blank")) %>%
  mutate(sdPS_ng_l   = NA,
         location    = word(loc_site, 1, sep = "_"),
         site = word(loc_site, 2, sep = "_")) %>%
  mutate(campaign = case_when(location %in% c("air", "oli") ~ "fie1",
                              location %in% c("bir", "kol") ~ "fie2",
                              TRUE ~ "fie3")) %>%
  mutate(meanPS_ng_l = case_when(uchem_id == 303 & campaign == "fie1" ~ 54.89,
                                 uchem_id == 303 & campaign == "fie2" ~ 8.14,
                                 uchem_id == 303 & campaign == "fie3" ~ 17.88,
                                 x4lettercode == "DICL" & campaign == "fie1" ~ 0.51,
                                 x4lettercode == "DICL" & campaign == "fie2" ~ 0.47,
                                 x4lettercode == "DICL" & campaign == "fie3" ~ 1.01,
                                 uchem_id == 205 & campaign == "fie1" ~ 2.39,
                                 uchem_id == 205 & campaign == "fie2" ~ 0.44,
                                 uchem_id == 205 & campaign == "fie3" ~ 2.45,
                                 uchem_id == 213 & campaign == "fie1" ~ 26.47,
                                 uchem_id == 213 & campaign == "fie2" ~ 35.45,
                                 uchem_id == 213 & campaign == "fie3" ~ 13.92,
                                 TRUE ~ NA))

df_names <- read.csv("data/namelist_basic.csv", sep = ";") %>%
  clean_names() %>%
  select(substance, x4lettercode, uchem_id)

site<- c("blank", "up1", "up2", "dw1", "dw2", "dw3")
location <- c("air","oli","bir","kol","alt","ehr")
loc_order <- unlist(lapply(location, paste, site, sep = "_"))

data <- df_names %>%
  left_join(bind_rows(df_PS, df_PS_blanks), by = c("uchem_id", "x4lettercode")) %>%
  group_by(substance) %>%
  mutate(PS_norm = meanPS_ng_l/max(meanPS_ng_l, na.rm = TRUE), 
         sdPS_norm = sdPS_ng_l/max(meanPS_ng_l, na.rm = TRUE)) %>%
  ungroup()%>%
  select(!c(sdPS_ng_l, sdPS_norm)) %>%
  mutate(site = factor(site, levels = c("blank", "up1", "up2", "dw1", "dw2", "dw3")), 
         location = factor(location, levels = c("air","oli","bir","kol","alt","ehr"))) %>%
  mutate(code = paste(substance, " - ", x4lettercode)) %>%
  mutate(code = factor(code, levels = sort(unique(code))), 
         loc_site = factor(loc_site, levels = loc_order)) %>%
  arrange(loc_site)

# Make matrix
matrix_PS <- data %>%
  select(code, loc_site, PS_norm) %>%
  pivot_wider(names_from = loc_site, values_from = PS_norm) %>%
  column_to_rownames("code") %>%
  as.matrix()

# Draw heatmap
ht <- Heatmap(matrix_PS,
              name = "Conc. norm.",
              col = colorRamp2(c(0, 1), c("white", "darkred")),
              na_col = "white",
              border_gp = gpar(col = "black", lty = 1, lwd = 2),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = TRUE,   
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              heatmap_legend_param = list(direction = "horizontal", title_position = "leftcenter"))

ht

png("plots/PS_heatmap.png", 
    width = 12, height = 15, units = "in", res = 150)
draw(ht, heatmap_legend_side = "bottom")

decorate_heatmap_body("Conc. norm.", {
  site <- sapply(strsplit(colnames(matrix_PS), "_"), `[`, 2)
  n <- ncol(matrix_PS)
  
  idx_blank <- which(site == "blank")
  idx_up2   <- which(site == "up2")
  idx_dw3   <- which(site == "dw3")
  
  draw_line_after <- function(i, gp) {
    x <- i / n  
    grid.lines(x = unit(c(x, x), "npc"),
               y = unit(c(0, 1), "npc"),
               gp = gp)
  }
  for(i in idx_blank) draw_line_after(i, gpar(col="green4", lty=2, lwd=2))
  for(i in idx_up2)   draw_line_after(i, gpar(col="black", lty=2, lwd=2))
  for(i in idx_dw3)   draw_line_after(i, gpar(col="black", lty=1, lwd=2))
})

dev.off()
