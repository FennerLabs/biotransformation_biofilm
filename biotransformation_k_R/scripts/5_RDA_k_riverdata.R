library(tidyverse)
library(vegan)
library(glue)
library(readxl) 
library(janitor) 
library(stringr)
library(tibble)

rm(list = ls())

set.seed(123)

# Set working directory
setwd("/biotransformation_k_R/")

# rate constant data frame
df_k_in <- read_csv("data/bf_tot_all_kderi.csv") %>% 
  clean_names() %>% 
  filter(!substance %in% "Indomethacin") %>% 
  select(loc_batch, substance, k) %>%
  pivot_wider(id_cols = loc_batch, names_from = substance, values_from = k) %>%
  column_to_rownames(var = "loc_batch") 

names(which(colSums(is.na(df_k_in)) > 0)) #check, which MPs have NAs

df_k <- df_k_in %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>% #replaces NA by the median of each the columns
  log1p() %>% #log(1 + x) avoids log(0)
  decostand(method="standardize") %>% # Standardize each substance (column) to mean=0, sd=1
  as.matrix()

# environmental data frame
df_riverT_in <- read_tsv("data/riverT_summary.tsv") %>%
  rename(loc_site = site, 
         "Temperature" = mean_temp) %>%
  select("loc_site", "Temperature") %>%
  mutate(Temperature = case_when(loc_site == "kol_up2" ~ 16.46,
                                 loc_site == "oli_dw3" ~ 8.99,
                                 TRUE ~ Temperature)) #replace NA values by nearest site
#biomass
df_biodata <- read_xlsx("data/biodata_allcampaigns.xlsx") %>% 
  clean_names() %>%
  mutate(location = word(loc_batch, 1, sep = "\\_"), 
         site = word(loc_batch, 2, sep = "\\_"),
         replicate = word(loc_batch, 3, sep = "\\_")) %>%
  rename(tss_0 = tss_g_l_0, tss_48 = tss_g_l_48, tss_96 = tss_g_l_96,
         afdm_0 =  afdm_g_l_0, afdm_48 = afdm_g_l_48,  afdm_96 = afdm_g_l_96,
         icc_0 = icc_n_l_1_0, icc_48 = icc_n_l_1_48, icc_96 = icc_n_l_1_96) %>%
  select(loc_batch, tss_0, afdm_0, icc_0) %>%
  mutate(tss_0 = case_when(loc_batch == "oli_up1_r1" ~ tss_0[match("oli_up1_r2", loc_batch)], 
                           loc_batch == "oli_dw3_r1" ~ tss_0[match("oli_dw3_r2", loc_batch)], 
                           TRUE ~ tss_0),
         afdm_0 = case_when(loc_batch == "oli_up1_r1" ~ afdm_0[match("oli_up1_r2", loc_batch)], 
                            loc_batch == "oli_dw3_r1" ~ afdm_0[match("oli_dw3_r2", loc_batch)], 
                            TRUE ~ afdm_0)) %>%
  mutate(across(c(tss_0, afdm_0, icc_0), ~ round(as.numeric(.x), 3))) %>%
  rename(TSS = tss_0,
         VSS = afdm_0,
         ICC = icc_0)

#river parameters
df_riverPara_in <- read_xlsx("data/riverParameter_summary.xlsx") %>%
  clean_names() %>%
  mutate(parameter = paste(parameter, "[", einheit, "]")) %>%
  select(-c(angabe, einheit)) %>%
  t() %>%
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  {df <- .
    df[] <- lapply(df, as.numeric)
    df[] <- lapply(seq_along(df), function(j) { # Replace NAs with LOQ values (first row of each column)
      col <- df[[j]]
      loq <- df[1, j]   # first row = LOQ
      ifelse(is.na(col), loq, col)
    })
  as.data.frame(df)} %>%
  rownames_to_column(var = "camp_loc_site") %>%
  mutate(across(-camp_loc_site, ~ round(as.numeric(.), 3))) %>%
  mutate(location = word(camp_loc_site, 2, sep = "_"),
         site = word(camp_loc_site, 3, sep = "_"),
         site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")),
         location = factor(location, levels = c("air", "oli", "bir", "kol", "alt", "ehr"))) %>%
  filter(camp_loc_site != "loq_value") %>%
  mutate(loc_site = paste(location, site, sep = "_")) %>%
  select(-c(camp_loc_site, location, site)) %>%
  set_names(c("Conductivity", "pH", "Alkalinity", "Tot_Hardness",
              "Sodium", "Magnesium", "Calcium", "Potassium", "Fluoride",
              "Chloride", "Bromide", "Nitrate", "Sulfate", "Ammonium",
              "Nitrite", "Phosphate", "TP", "Silicic_Acid", "DOC",
              "TOC", "TN", "TIC", "loc_site")) %>%
  select(-Bromide)

df_env_in <- left_join(df_riverPara_in, df_riverT_in, by = "loc_site")

df_env <- df_env_in %>%
  bind_rows(df_env_in, .id = "replicate") %>%   # duplicate and mark copy
  mutate(loc_batch = paste0(loc_site, if_else(replicate == 1, "_r1", "_r2"))) %>%
  select(-c(replicate, loc_site)) %>%
  left_join(df_biodata, by = "loc_batch") %>%
  mutate(Location = word(loc_batch, 1, sep = "_"),
         Site = ifelse(str_detect(loc_batch, "up"), "up", "dw"),
         Camp = case_when(Location %in% c("air", "oli") ~ "fie1",
                          Location %in% c("bir", "kol") ~ "fie2",
                          TRUE ~ "fie3")) %>%
  column_to_rownames(var = "loc_batch") %>%
  {df <- .
  num_cols <- df %>% select(-Site, -Location, -Camp) %>% as.data.frame()
  cat_cols <- df %>% select(Site, Location, Camp)
  num_cols <- scale(num_cols) %>% as.data.frame()
  cat_cols <- cat_cols %>% mutate(across(everything(), factor))
  bind_cols(num_cols, cat_cols)
  } %>%
  select(!c(TSS, VSS, ICC))

# Run RDA
rda_model <- rda(df_k ~ ., data = df_env)

# Inspect results
summary(rda_model)           # summary of RDA
anova(rda_model)             # test overall significance
anova_axis <- anova(rda_model, by="axis")  # test significance per RDA axis
print(anova_axis)
anova_term <- anova(rda_model, by="term")  # test significance per environmental variable
print(anova_term)

#to find, for which the models alias
X_full <- model.matrix(~ ., data=df_env)  # all predictors
qr(X_full)$rank
ncol(X_full)
alias(lm(df_k[,1] ~ ., data=df_env))      # show which columns are linearly dependent

# Adjusted R²
RsquareAdj(rda_model)

# Visualization of RDA simple
plot(rda_model, scaling=2)  # scaling=2 emphasizes distances between samples

#Visualization custommade
colors_loc = c("air" = "#0000FF",  "oli" = "#00008B",  "bir" = "#43CD80", 
               "kol" = "#2E8B57", "alt" = "#FFA500",  "ehr" = "#CD6600")

site_sym = c("up1" = 3, "up2" = 7, "dw1" = 16, "dw2" = 17, "dw3" = 15)

# Extract scores
site_scores <- scores(rda_model, display = "sites", scaling = 2)
species_scores <- scores(rda_model, display = "species", scaling = 2)
env_scores <- scores(rda_model, display = "bp", scaling = 2)

df_meta <- df_env %>% 
  rownames_to_column("loc_batch") %>%
  mutate(replicate = ifelse(str_detect(loc_batch, "r1"), "r1", "r2"),
         location  = word(loc_batch, 1, sep = "_"),
         site      = word(loc_batch, 2, sep = "_"))

# Extract constrained eigenvalues (RDA axes)
rda_eig <- summary(rda_model)$cont$importance
rda_var <- rda_eig[2, 1:2] * 100  

#saving plot 
png(filename = "plots/RDA_k.png", width = 12, height = 12, units = "in", res = 300) #adjust name of plot!

plot(rda_model, scaling = 2, type = "n", 
     xlab = paste0("RDA1 (", round(rda_var[1],1), "%)"), 
     ylab = paste0("RDA2 (", round(rda_var[2],1), "%)"))

# Plot base site symbols (all samples)
points(site_scores[,1], site_scores[,2],
       col = colors_loc[df_meta$location],
       pch = site_sym[df_meta$site],
       cex = 1.5)

# Add circle around replicate r2
is_r2 <- df_meta$replicate == "r2"
points(site_scores[is_r2,1], site_scores[is_r2,2],
       col = colors_loc[df_meta$location][is_r2], 
       pch = 1,   
       cex = 2.1,     
       lwd = 1.8)

# only choose for significant (p < 0.05) environmental factors
sig_terms <- rownames(anova_term)[anova_term$`Pr(>F)` < 0.05]
env_sel <- env_scores[rownames(env_scores) %in% sig_terms, , drop = FALSE]

# Draw arrows with scaling 5.5
arrows(0, 0, env_sel[,1]*5.5, env_sel[,2]*5.5,         #adjust scaling: for all k. 5.5 for k>-2: 3
       col = "black", length = 0.1, lwd = 1.2)   

# Add arrow labels
text(env_sel[,1] *5.5* 1.1,  #adjust scaling: for all k. 5.5 for k>-2: 3
     env_sel[,2] *5.5* 1.1,
     labels = rownames(env_sel),
     col = "black", cex = 0.8, font = 2)

par(xpd = NA)  # allow drawing outside plot

# Locations
legend(x = 4, y = 2,
       legend = names(colors_loc),
       col = colors_loc, pch = 16, pt.cex = 1.5,
       bty = "n", title = "Location")

# Sites
legend(x = 4, y = 0,
       legend = names(site_sym),
       col = "black", pch = site_sym, pt.cex = 1.5,
       bty = "n", title = "Site")

# Replicate r2
legend(x = 4, y = -2,
       legend = "r2",
       col = "black", pch = 1, pt.cex = 2, pt.lwd = 1.8,
       bty = "n", title = "Replicate")

dev.off()

###################################################################################
#plot R2 values
df_r2 <- anova_term %>%
  as.data.frame() %>%
  rownames_to_column(var = "variable") %>%
  filter(variable != "Residual") %>%
  mutate(r2_percent = (Variance / sum(Variance)) * 100,
         significance = cut(`Pr(>F)`,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                            labels = c("***", "**", "*", ""),
                            right = TRUE)) %>%
  arrange(r2_percent) %>%
  mutate(variable = factor(variable, levels = variable))  

#saving plot 
png(filename = "plots/RDA_r2.png", width = 12, height = 12, units = "in", res = 300) #adjust name of plot!

bar_positions <- barplot(df_r2$r2_percent,
                         names.arg = df_r2$variable,
                         horiz = TRUE,
                         las = 1,
                         col = "skyblue",
                         xlab = "Variance explained [%]",
                         cex.names = 0.8,
                         xlim = c(0, 16))

text(x = df_r2$r2_percent + max(df_r2$r2_percent) * 0.02,
     y = bar_positions,
     labels = df_r2$significance,
     cex = 1.2,
     adj = 0)

dev.off()   