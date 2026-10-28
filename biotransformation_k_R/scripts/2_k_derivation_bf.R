library(tidyverse)
library(ggpubr)
library(readxl) 
library(janitor) 

rm(list = ls())

# Set working directory
setwd("/biotransformation_k_R/") #update the working directory

#########################################################################
# (1) Read .csv file

data_in <- read_csv("data/BF_tot_conc.csv") %>% 
  clean_names() %>%
  select(uchem_id, substance, smiles, cas, x4lettercode, location, site, type, replicate, time_point, conc_nmol) %>%
  mutate(site_replicate = paste(site, replicate, sep="_"),
         batch = paste(site, replicate, type, sep="_"),
         loc_site = paste(location, site, sep = "_"),
         loc_batch = paste(location, site, replicate, sep = "_"),
         site_replicate_type_time_point = paste(site, replicate, type, time_point, sep ="_"),
         ln_conc_nmol = log(conc_nmol))

df_selection_manual <- read.csv("data/substance_exclusion_manual_afterct_BF.csv", sep = ";") %>%
  clean_names() %>%
  select(substance, substance_class, own_istd, log_koc_4, abiotic_transformation, cationic_at_p_h8, 
         bioaccumulation_anna, sorption_visual, iontrapping_visualdigi_bf, possible_normalization_bt_sc_or_watch_out_sorption_flag_difference_bt_sc_large,
         final_exlcusion) %>%
  rename(final_excl = final_exlcusion,
         sorption_flag = possible_normalization_bt_sc_or_watch_out_sorption_flag_difference_bt_sc_large)

data <- left_join(data_in, df_selection_manual, by = "substance") %>% 
  group_by(substance, loc_batch, type, batch) %>%
  arrange(factor(time_point, levels = c(-1, 0, 1, 2, 4, 8, 24, 48, 96))) %>%
  mutate(background = case_when(first(conc_nmol > 0) ~ TRUE, 
                                first(is.na(conc_nmol)) ~ FALSE,
                                TRUE ~ FALSE)) %>%
  ungroup() %>%
  filter(final_excl == "FALSE",
         !(time_point == "-1")) %>%
  select(!c(site_replicate, batch, site_replicate_type_time_point, final_excl))

################################################################################
#sorption labeling
#label compounds showing more than 30% sorption (below 7nM) within 24h (t=0-24) in type "sc"
bf_sorp <- data %>%
  filter(!(time_point %in% c(48, 72, 96)) & type == "sc")%>%
  group_by(loc_batch, substance) %>%
  mutate(min_conc = min(conc_nmol, na.rm = TRUE)) %>%
  mutate(sorp_control = case_when(min_conc < 7 ~ TRUE,
                          TRUE ~ FALSE)) %>%
  slice(1) %>%
  ungroup() %>%
  select(substance, sorp_control, loc_batch)
  
bf_tot <-  left_join(data, bf_sorp, by = c("substance", "loc_batch"), multiple ="all") %>%
  unique() %>%
  filter(!(type == "sc"))

################################################################################
#removes substances per site, which do not have any or just one data point(s)
bf_tot_error <- bf_tot %>%
  group_by(substance, location, loc_batch, replicate) %>%
  mutate(n = sum(!is.na(conc_nmol))) %>%
  ungroup() %>%
  mutate(comp_loc_batch = paste(substance, loc_batch, sep = "_"),
         comp_loc = paste(substance, location, sep = "_"))

###################################################################################
#manual correction of certain values
bf_tot_error_manualcorrection <- bf_tot_error %>%
  mutate(conc_nmol = case_when(loc_batch == "air_up1_r2" & time_point == "48" ~ NA_real_,
                               substance == "Acesulfame" & loc_batch == "oli_dw2_r1" & time_point %in% c("8", "48") ~ NA_real_, 
                               substance == "Acetamiprid" & loc_batch == "bir_up2_r1" & time_point == "4" ~ NA_real_,
                               substance == "Acetamiprid" & loc_batch == "alt_dw1_r1" & time_point == "24" ~ NA_real_,
                               substance == "Alachlor" & loc_batch == "bir_dw3_r2" & time_point == "4" ~ NA_real_,
                               substance == "Albuterol" & loc_batch == "air_up1_r2" & time_point == "24" ~ NA_real_,
                               substance == "Albuterol" & loc_batch == "bir_up2_r1" & time_point == "2" ~ NA_real_,
                               substance == "Capecitabine" & loc_batch %in% c("ehr_dw2_r2", "ehr_dw3_r1") & time_point == "2" ~ NA_real_,
                               substance == "Carbamazepine" & loc_batch == "bir_up2_r1" & time_point == "8" ~ NA_real_,
                               substance == "Carbetamide" & loc_batch == "kol_up2_r2" & time_point == "0" ~ NA_real_,
                               substance == "Chlorotoluron" & loc_batch == "kol_dw1_r2" & time_point == "0" ~ NA_real_,
                               substance == "Clomazone" & loc_batch == "alt_dw1_r2" & time_point == "4" ~ NA_real_,
                               substance == "Crotamiton" & loc_batch == "air_dw2_r2" & time_point %in% c("0", "2") ~ NA_real_,
                               substance == "Darunavir" & loc_batch == "oli_up2_r2" & time_point == "0" ~ NA_real_,
                               substance == "Doxylamine" & loc_batch == "kol_dw2_r1" & time_point == "4" ~ NA_real_,
                               substance == "Enzalutamide" & loc_batch == "alt_dw1_r2" & time_point == "2" ~ NA_real_,
                               substance == "Etodolac" & loc_batch %in% c("alt_dw1_r1", "alt_dw1_r2") & time_point == "24" ~ NA_real_,
                               substance == "Etodolac" & loc_batch == "alt_dw2_r1" & time_point %in% c("2", "4") ~ NA_real_,
                               substance == "Fenhexamid" & loc_batch == "kol_dw2_r1" & time_point == "48" ~ NA_real_,
                               substance == "Flufenacet" & loc_batch == "alt_dw1_r2" & time_point == "2" ~ NA_real_,
                               substance == "Ganciclovir" & loc_batch %in% c("kol_dw2_r1", "kol_dw2_r2") & time_point == "4" ~ NA_real_,
                               substance == "Ganciclovir" & loc_batch == "air_dw1_r2" & time_point == "24" ~ NA_real_,
                               substance == "Ganciclovir" & loc_batch == "bir_up1_r2" & time_point == "1" ~ NA_real_,
                               substance == "Ganciclovir" & loc_batch == "ehr_dw3_r1" & time_point == "8" ~ NA_real_,
                               substance == "Ketamine" & loc_batch == "oli_up1_r1" & time_point == "72" ~ NA_real_,
                               substance == "Ketamine" & loc_batch == "kol_dw2_r1" & time_point == "96" ~ NA_real_,
                               substance == "Lamivudine" & loc_batch %in% c("air_dw1_r1", "air_up2_r2") & time_point == "2" ~ NA_real_,
                               substance == "Lamivudine" & loc_batch == "bir_up1_r2" & time_point == "8" ~ NA_real_,
                               substance == "Lamivudine" & loc_batch == "bir_up1_r1" & time_point == "4" ~ NA_real_,
                               substance == "Lamivudine" & loc_batch == "kol_dw2_r1" & time_point == "96" ~ NA_real_,
                               substance == "Lisinopril" & loc_batch == "air_dw1_r2" & time_point == "24" ~ NA_real_,
                               substance == "Mirabegron" & loc_batch == "bir_up2_r1" & time_point == "1" ~ NA_real_,
                               substance == "Mirtazapine" & loc_batch %in% c("bir_dw1_r1", "bir_dw2_r1") & time_point == "4" ~ NA_real_,
                               substance == "Mirtazapine" & loc_batch == "bir_up2_r1" & time_point == "2" ~ NA_real_,
                               substance == "Oseltamivir" & loc_batch == "alt_dw1_r2" & time_point == "2" ~ NA_real_,
                               substance == "Prometryn" & loc_batch == "alt_dw1_r2" & time_point == "2" ~ NA_real_,
                               substance == "Propachlor" & loc_batch == "oli_dw3_r1" & time_point == "0" ~ NA_real_,
                               substance == "Rosuvastatin" & loc_batch == "air_dw1_r1" & time_point == "2" ~ NA_real_,
                               substance == "Rosuvastatin" & loc_batch == "air_dw1_r2" & time_point == "24" ~ NA_real_,
                               substance == "Rosuvastatin" & loc_batch == "kol_up2_r2" & time_point == "4" ~ NA_real_,
                               substance == "Saccarin" & loc_batch %in% c("air_dw1_r2", "alt_dw1_r1") & time_point == "24" ~ NA_real_,
                               substance == "Tebutam" & loc_batch %in% c("oli_dw1_r1", "oli_dw3_r2") & time_point == "2" ~ NA_real_,
                               substance == "Tramadol" & loc_batch == "bir_dw2_r1" & time_point == "1" ~ NA_real_,
                               TRUE ~ conc_nmol))%>%
  mutate(ln_conc_nmol = log(conc_nmol)) 
  
###################################################################################
bf_tot_exclKderi <- bf_tot_error_manualcorrection %>%
  group_by(substance, loc_batch) %>%
  slice(1)%>%
  filter(n < 2) %>%
  ungroup() %>%
  select(!c(time_point))

###################################################################################
bf_error <- bf_tot_error %>%
  filter(n >= 2)

# Derivation of k for each replicate
# Model
model <- function(df) lm(ln_conc_nmol ~ time_point, data = df, na.action=na.omit)
lm_summary <- function(mod) summary(mod)
lm_slope <- function(mod) coef(mod)["time_point"]
lm_intercept <- function(mod) coef(mod)["(Intercept)"]
lm_R2 <- function(mod) summary(mod)$r.squared
lm_p <- function(mod) summary(mod)$coefficients["time_point",4]

bf_k <- bf_error %>%
  group_by(substance, x4lettercode, loc_batch) %>%
  reframe(slope = lm_slope(model(pick(ln_conc_nmol, time_point))),
          intercept = lm_intercept(model(pick(ln_conc_nmol, time_point))),
          R2 = lm_R2(model(pick(ln_conc_nmol, time_point))),
          p = lm_p(model(pick(ln_conc_nmol, time_point)))) %>%
  ungroup() %>%
  mutate(k = -slope, # 1st order rate constant [1/h]
         t05 = log(2)/k) # half-life [h]

bf_k_full <- left_join(bf_error, bf_k,
                   by = c("substance", "x4lettercode", "loc_batch")) %>%
  group_by(substance, comp_loc, loc_batch) %>%
  slice(1) %>%
  ungroup () %>%
  select(!("time_point"))

missing_cols <- setdiff(names(bf_k_full), names(bf_tot_exclKderi))

for (col in missing_cols) {
  bf_tot_exclKderi[[col]] <- NA
}

bf_tot_exclKderi <- bf_tot_exclKderi[names(bf_k_full)]

bf_combined <- rbind(bf_k_full, bf_tot_exclKderi) %>%
  select(!c("ln_conc_nmol", "loc_site", "type", "conc_nmol", "comp_loc")) %>%
  rename(n_time_points = n)

substanceselection <- unique(bf_combined$substance)

df_missing_loc <- read_csv("data/namelist_location_bf_full_167.csv") %>%
  clean_names() %>%
  filter(substance %in% substanceselection) %>%
  mutate(loc_batch = paste(location, site, replicate, sep = "_")) %>%
  select(!c(x1, location, replicate)) 

sitesname <- df_missing_loc %>%
  group_by(site) %>%
  slice(1) %>%
  pull(site)

bf_combined <- bf_combined %>%
  filter(site %in% sitesname)

bf_combined_fullloc <- left_join(df_missing_loc, bf_combined, by = c("substance", "loc_batch", "site")) %>%
  group_by(substance) %>%
  select(!comp_loc_batch) %>%
  mutate(across(c(cas, smiles, x4lettercode, substance_class, own_istd, log_koc_4, 
                  abiotic_transformation, cationic_at_p_h8, bioaccumulation_anna, 
                  iontrapping_visualdigi_bf, sorption_visual, sorption_flag),
                ~ ifelse(is.na(.), first(.[!is.na(.)]), .))) %>%
  mutate(location = word(loc_batch, 1, sep = "\\_"),
         replicate = word(loc_batch, 3, sep = "\\_"))

#label persistent (kmin) and kmax substances
#label compounds showing less than 10% dissipation (c(t)/c0 = e^-kt = 9/10, t=96h, exclude if k<= -ln(0.9)/96 [h])
#label compounds over kmax -> At least 3 time points with known concentration (more than 10%): k>= -ln(0.1)/2 [h]
bf_label1 <- bf_combined_fullloc %>%
  mutate(kmin = ifelse(k < 0.005, TRUE, FALSE),
         kmax = ifelse(k > 2.303 | n_time_points < 3, TRUE, FALSE),
         nonvalid = ifelse(k > 0.005 & R2 < 0.8 & p > 0.05, TRUE, FALSE))

###################################################################################
#manual correction of k values based on visual ct consultation
bf_label <- bf_label1 %>%
  mutate(across(c(R2, p, k, n_time_points, slope, intercept, t05), ~ case_when(substance == "Carbamazepine" & loc_batch %in% c("bir_dw3_r2", "bir_dw3_r1", "bir_up2_r2") ~ NA_real_,
                                                                               substance == "Cilastatin" & loc_batch %in% c("oli_up1_r2", "bir_dw3_r1", "bir_dw3_r2") ~ NA_real_,
                                                                               substance == "Doxylamine" & loc_batch %in% c("air_dw3_r1", "air_dw3_r2") ~ NA_real_,
                                                                               substance == "Enzalutamide" & loc_batch %in% c("alt_dw2_r1", "alt_dw1_r1") ~ NA_real_,
                                                                               substance == "Ganciclovir" & location == "oli" ~ NA_real_,
                                                                               substance == "Metolachlor" & loc_batch == "bir_dw1_r1" ~ NA_real_,
                                                                               substance == "Mirabegron" & location %in% c("kol", "air") ~ NA_real_,
                                                                               substance == "Mirtazapine" & location %in% c("ehr", "air") ~ NA_real_,
                                                                               substance == "Napropamide" & loc_batch == "kol_dw2_r1" ~ NA_real_,
                                                                               substance == "Naproxen" & loc_batch %in% c("air_dw2_r2", "ehr_dw3_r2") ~ NA_real_,
                                                                               substance == "Naproxen" & location == "kol" ~ NA_real_,
                                                                               substance == "Normorphine" & location == "oli" ~ NA_real_,
                                                                               substance == "Omeprazole" & location == "air" ~ NA_real_,
                                                                               substance == "Paracetamol" & location == "air" ~ NA_real_,
                                                                               substance == "Pargyline" & location == "oli" ~ NA_real_,
                                                                               substance == "Pseudoephedrine" &location == "air" ~ NA_real_,
                                                                               substance == "Quetiapine" & location == "air" ~ NA_real_,
                                                                               substance == "Rivaroxaban" & loc_batch %in% c("bir_dw3_r1", "bir_dw3_r2") ~ NA_real_,
                                                                               substance == "Saccharin" & loc_batch %in% c("kol_up1_r2", "kol_up2_r1") ~ NA_real_,
                                                                               substance == "Tenofovir" & location == "oli" ~ NA_real_,
                                                                               substance == "Terbuthylazine" & loc_batch == "kol_dw2_r1" ~ NA_real_,
                                                                               TRUE ~ .))) %>%
  mutate(across(c(R2, p, slope, intercept, t05), ~ ifelse(kmax == TRUE, NA_real_, .)))

write_csv(bf_label, "data/bf_tot_all_kderi.csv")


