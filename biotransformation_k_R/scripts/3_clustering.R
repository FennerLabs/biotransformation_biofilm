library(tidyverse)
library(glue)
library(dendextend)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(readxl) 
library(janitor) 
library(RColorBrewer)
library(rstatix)
library(grid)
library(stringr)
library(cluster)
library(ggrepel)

set.seed(123)

rm(list = ls())

# Set working directory
setwd("/biotransformation_k_R/") #update the working directory

#functions
to_annot_sub <- function(df, rownames = "x4lettercode") {
  df_ann <- df %>%
    select(-all_of(rownames)) %>%
    as.data.frame()
  rownames(df_ann) <- df %>% pull(rownames)
  return(df_ann)
}

to_annot_locbatch <- function(df, rownames = "loc_batch") {
  df_ann <- df %>%
    select(-all_of(rownames)) %>%
    as.data.frame()
  rownames(df_ann) <- df %>% pull(rownames)
  return(df_ann)
}

to_matrix <- function(df, rownames = "x4lettercode") {
  df_mat <- df %>%
    select(-all_of(rownames)) %>%
    as.matrix()
  rownames(df_mat) <- df %>% pull(rownames)
  return(df_mat)
}

#define kmin
kminnum <- 0.005
#define kmax: more than 90% dissipation over two time points (1h)
kmaxnum <- -log(0.1)/1
#define validation criteria
R2val <- 0.80
pval <- 0.05

#################################################################################
#substance categories based on DOI:10.1016/j.watres.2025.125283
persistent <- c("Ifosfamide", "Crotamiton", "Mesotrione", "Diuron", "Amisulpride", "Metoxuron", "Atrazine", 
                "Bentazon", "Carbamazepine", "Carbendazim", "Cyclophosphamide", "Lamotrigine", "Monuron", 
                "Olmesartan", "Perindopril", "Primidone", "Simeton", "Vildagliptin", "Dimefuron", "Terbuthylazine", 
                "Torsemide", "Tebufenozide", "Isoproturon", "Terbutryn", "Prometryn", "Flonicamid", "Napropamide", 
                "4-Chloro-N,N-dimethylbenzamide", "Enzalutamide")

cometa <- c("Ranitidine", "Pargyline", "Metolachlor", "Trazodone", "Metoprolol", "Flufenacet", "Dimethenamid",
            "Metoclopramide", "Quetiapine", "Mexiletine", "Diphenhydramine", 
            "Vorinostat", "Emtricitabine", "Atenolol", "Alachlor", "Mirtazapine", "Sulfathiazole", "Fenhexamid",
            "Pioglitazone", "Capecitabine", "Pseudoephedrine", "Deprenyl", "Cilastatin", "Propachlor",
            "Flufenamic acid",  "Clomethiazole", "Sulfapyridine",
            "Sulfadiazine", "Fenoxycarb", "Sulfamethoxazole", "Rivastigmine", "Lamivudine",
            "Trinexapac-ethyl", "Acetamiprid", "Sulfamethazine", "Venlafaxine", "Bromoxynil", "Sacubitril",
            "Bezafibrate", "Sotalol", "Codeine",
            "Icaridin", "Rivaroxaban", "Naloxone", "Abacavir", "Naloxegol", "Doxylamine", 
            "Rufinamide", "Albuterol", "Mirabegron", "Apalutamide", "Sitagliptin", "Fludioxonil")

highinf <- c("Caffeine", "Levetiracetam",  "Saccharin", "Valsartan", "Cyclamate", "Paracetamol",  
             "Lisinopril", "Mefenamic acid", "Tenofovir")

cata <- c("Tramadol", "Ketamine", "Normorphine", "Gliclazide", "Morphine", "Naproxen", "Ganciclovir", "Darunavir",
          "Ibuprofen", "Levamisole", "Losartan", "Picoxystrobin", "Chlorotoluron", "Mycophenolic acid", "Carbetamide",
          "Moclobemide", "Etodolac", "Sucralose", "Tebutam", "Iprovalicarb", "Fexofenadine", "Benzophenone",
          "Clopidogrel carboxylic acid", "Clomazone",  "Diazepam", "DEET",  "Furosemide", "Pravastatin", "Ketoprofen",
          "Mecoprop", "Clofibric acid", "Diclofenac", "Pirfenidone", "Acesulfame", "Rosuvastatin") 

#################################################################################
# Read and write tables
data_sub_class <- read_csv2("data/namelist_cluster.csv") %>% 
  clean_names() %>%
  select(substance, substance_class) 

data_in <- read_csv("data/bf_tot_all_kderi.csv") %>% 
  clean_names() %>%
  select(-substance_class) %>%
  left_join(data_sub_class, by = "substance") %>%
  select(!c("slope", "intercept", "t05")) %>%
  mutate(neg_k = ifelse(k < 0, TRUE, FALSE),
         classification = case_when(substance_class == "Pharmaceutical" ~ "PH",
                                    substance_class == "n.a." ~ "HI",
                                    substance_class == "Personal care" ~ "HI",
                                    substance_class == "Fungicide" ~ "PE-F",
                                    substance_class == "Herbicide" ~ "PE-H",
                                    substance_class == "Artificial sweetener" ~ "AS",
                                    substance_class == "Insecticide" ~ "PE-I",
                                    substance_class == "Antibiotic" ~ "AB", 
                                    TRUE ~ NA),
         substance_class = case_when(substance_class == "n.a." ~ "Household and Industrial Chemical",
                                     substance_class == "Personal care" ~ "Household and Industrial Chemical",
                                     TRUE ~ substance_class),)

df_names <- read.csv("data/namelist_basic.csv", sep = ";") %>%
  clean_names()

#substance categories based on control experiments described in DOI:10.1016/j.watres.2025.125283
#substances with measurement issues
mi_sub <- c("Allopurinol", "Asulam", "Benserazide", "Clotrimazole", "Canagliflozin", "Dasatinib", "Dolutegravir", "Gemfibrozil",
            "Hydroxychloroquine", "Ibuprofen", "Lovastatin", "Mesalazine", "Metformin", "Neohesperidin dihydrochalcone", "Nintedanib",
            "Olanzapine", "Piracetam", "Sertraline", "Simvastatin", "Spironolactone", "Terbinafine", "Regorafenib")

#substances with abiotic degradation
ab_sub <- c("Acemetacin", "Aliskiren", "Apremilast", "Atorvastatin", "Bisacodyl", "Dexlansoprazole", "Diltiazem", 
            "Edoxaban", "Irbesartan", "Kresoxim-methyl", "Neotame", "Omeprazole", "Tolperisone", "Valaciclovir")

#sustances with sorption below LOQ (not measurably detectable)
sorp_sub_LOQ <- c("Amlodipine", "Atovaquone", "Bupropion", "Clomipramine", "Clozapine", "Diosmin", "Duloxetine", 
                  "Efavirenz", "Elvitegravir", "Entacapone", "Fluoxetine", "Ticlopidine", 
                  "Trimipramine") 

#substances with high sorption and not differentiable from biotransformation
sorp_sub_noBT <- c("Apixaban", "Atomoxetine", "Climbazole", "Cyromazine", "Dextromethorphan", "Diphenhydramine", "Enzalutamide", 
                   "Fenfluramine", "Fipronil", "Flonicamid", "Flufenamic acid", "Irgarol", "Phenytoin",
                   "Propranolol", "Tebufenozide", "Sucralose")

#bioaccumulative substances in biofilms 
bioacc_sub <- c("Trazodone", "Sulfasalazine")

#substances where sorption takes place, but the biotransformation is distinguishable
sorpflag_sub <- c("Apalutamide", "Benzophenone", "Deprenyl", "Febuxostat", "Fenfluramine", "Fenhexamid", "Fenoxycarb", 
                  "Flufenacet", "Levamisole", "Mefenamic acid", "Metoclopramide", "Metolachlor", "Mexiletine", "Mirtazapine", 
                  "Moclobemide", "Morphine", "Naloxegol", "Naloxone", "Normorphine", "Picoxystrobin", "Prometryn", "Quetiapine", 
                  "Ranitidine", "Rivaroxaban", "Sitagliptin", "Tenofovir", "Terbuthylazine", "Terbutryn", "Trazodone", 
                  "Trinexapac-ethyl", "Vorinostat")

bt_sub <- df_names %>%
  filter(!substance %in% c(mi_sub, ab_sub, sorp_sub_LOQ, sorp_sub_noBT)) %>%
  pull(substance)

#load data for river parameters describing each biofilm sample (site specific)
df_riverPara_in <- read_xlsx("data/riverParameter_summary.xlsx") %>% 
  clean_names() %>% 
  mutate(parameter = paste(parameter, " [", einheit, "]", sep = "")) %>%
  select(!c(angabe, einheit)) %>%
  t() %>% 
  as.data.frame() %>% 
  {colnames(.) <- .[1, ]   
  .[-1, ]     } %>%  
  rownames_to_column(var = "camp_loc_site") %>%
  mutate(across(!c(camp_loc_site), ~ round(as.numeric(.), 3))) %>%
  mutate(location = word(camp_loc_site, 2, sep = "\\_"),
         site = word(camp_loc_site, 3, sep = "\\_")) %>%
  mutate(site = factor(site, levels = c("up1", "up2", "dw1", "dw2", "dw3")),
         location = factor(location, levels = c("air", "oli", "bir", "kol", "alt", "ehr"))) %>%
  filter(!camp_loc_site == "loq_value") %>%
  mutate(loc_site = paste(location, site, sep = "_")) %>%
  select(!c(camp_loc_site, location, site))

#translation for the colnames
colnames(df_riverPara_in) <- c("Conductivity [µS/cm 20°C]", "pH", "Alkalinity [mmol/L]", "Total Hardness [mmol/L]", "Sodium [mg/L]", 
                       "Magnesium [mg/L]", "Calcium [mg/L]", "Potassium [mg/L]", "Fluoride [mg/L]", "Chloride [mg/L]", "Bromide [mg/L]",
                       "Nitrate [mg/L]", "Sulfate [mg/L]", "Ammonium [µg/L]", "Nitrite [µg/L]", "Phosphate [µg/L]", "TP [µg/L]",
                       "Silicic Acid [mg/L]", "DOC [mg/L]", "TOC [mg/L]", "TN [mg/L]", "TIC [mg/L]", "loc_site")
  
df_riverPara <- df_riverPara_in %>%
  select(!c("Fluoride [mg/L]","Chloride [mg/L]", "Bromide [mg/L]", "Phosphate [µg/L]")) %>%
  mutate(across(-c(loc_site, pH), ~ round(. / max(., na.rm = TRUE), 2)), #normate to values between 0 and 1, round to 2 digits
         pH = round(pH/14, 2)) 

df_riverT <- read_tsv("data/riverT_summary.tsv") %>%
  rename(loc_site = site, 
         "Temperature" = mean_temp) %>%
  select("loc_site", "Temperature") %>%
  mutate(Temperature = case_when(loc_site == "kol_up2" ~ 16.46,
                                 loc_site == "oli_dw3" ~ 8.99,
                                 TRUE ~ Temperature)) #replace NA values by nearest site (HOBO loggers got lost during field campaign)

#load data for passive samplers
df_ps <- read_csv("data/Result_influent_effluent_passiveSampling_complete.csv") %>% 
  clean_names() %>%
  rename(x4lettercode = lettercode4, 
         substance = name) %>%
  filter(type == "PassiveSampling") %>% 
  mutate(location = word(sample, 2, sep = "_"),
         site = word(sample, 3, sep = "_"),
         loc_site = paste(location, site, sep = "_")) %>%
  group_by(substance, sample) %>%
  mutate(meanPS_ng_l = mean(c(ps_extract_ng_l_r1, ps_extract_ng_l_r2), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(meanPS_ng_l = ifelse(meanPS_ng_l == 0, NA, meanPS_ng_l)) %>%
  select(uchem_id, x4lettercode, meanPS_ng_l, loc_site, location, site) %>%
  group_by(x4lettercode) %>%
  mutate(PS_norm = meanPS_ng_l/max(meanPS_ng_l, na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(site = ifelse(grepl("^up", site), "up", "dw")) %>%
  group_by(site, x4lettercode) %>%
  mutate(mean_ps_site_norm = mean(PS_norm, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  select(x4lettercode, site, mean_ps_site_norm) %>%
  mutate(mean_ps_site_norm = ifelse(mean_ps_site_norm == "NaN", NA, mean_ps_site_norm))

#load biomass data to describe biomass in each BF sample
df_biodata <- read_xlsx("data/biodata_allcampaigns.xlsx") %>% 
  clean_names() %>%
  mutate(location = word(loc_batch, 1, sep = "\\_"), 
         site = word(loc_batch, 2, sep = "\\_"),
         replicate = word(loc_batch, 3, sep = "\\_")) %>%
  rename(tss_0 = tss_g_l_0, tss_48 = tss_g_l_48, tss_96 = tss_g_l_96,
         afdm_0 =  afdm_g_l_0, afdm_48 = afdm_g_l_48,  afdm_96 = afdm_g_l_96,
         icc_0 = icc_n_l_1_0, icc_48 = icc_n_l_1_48, icc_96 = icc_n_l_1_96) %>%
  select(location, loc_batch, site, replicate, tss_0, afdm_0, icc_0, chl_a_mg_m_l_1, od_1) %>%
  mutate(tss_0 = case_when(loc_batch == "oli_up1_r1" ~ tss_0[match("oli_up1_r2", loc_batch)], 
                           loc_batch == "oli_dw3_r1" ~ tss_0[match("oli_dw3_r2", loc_batch)], 
                           TRUE ~ tss_0),
         afdm_0 = case_when(loc_batch == "oli_up1_r1" ~ afdm_0[match("oli_up1_r2", loc_batch)], 
                            loc_batch == "oli_dw3_r1" ~ afdm_0[match("oli_dw3_r2", loc_batch)], 
                            TRUE ~ afdm_0),
         icc_0 = ifelse(loc_batch %in% c("oli_up1_r1", "alt_dw3_r2"), NA, icc_0)) %>%
  mutate(across(c(tss_0, afdm_0, icc_0), ~ round(as.numeric(.x), 3))) %>%
  rename(chl_a = chl_a_mg_m_l_1,
         od = od_1) %>%
  mutate(od = ifelse(location %in% c("air", "oli"), od*2, od))
  
data <- data_in %>%
  filter(!substance == "Indomethacin") %>% #Indomethacin is a TP of the co-spiked Acemetacin
  select(x4lettercode,substance, loc_batch, site, location, replicate, substance_class, bioaccumulation_anna, 
         sorp_control, r2, p, k, kmin, kmax, classification, sorption_flag) %>%
  filter(substance %in% bt_sub) %>%
  mutate(loc_site = paste(location, site, sep = "_")) %>%
  mutate(nonvalid = ifelse(r2 < 0.8 & p > 0.05 & k > 0.005, TRUE, FALSE)) %>%
  mutate(kori = k,
         k = round(case_when(kmin ~ kminnum,
                       kmax ~ kmaxnum, 
                       nonvalid~ NA, 
                       TRUE ~ kori), 3)) %>%
  mutate(across(c(r2, p), ~ ifelse(kmin | kmax | nonvalid, NA, .))) %>%
  left_join(df_riverT, by = "loc_site") %>%
  left_join(df_biodata %>% select(loc_batch, afdm_0, tss_0, icc_0, chl_a, od), by = "loc_batch") %>%
  mutate(kafdm = k/afdm_0,
         kicc = k/icc_0,
         ktss = k/tss_0,
         kod = k/od) %>%
  mutate(logk = log10(k),
         logkafdm = log10(kafdm),
         logkicc = log10(kicc),
         logktss = log10(ktss),
         logkod = log10(kod)) %>%
  group_by(substance) %>%
  mutate(meank_sub = mean(logk, na.rm = TRUE),
         stdvk_sub = sd(logk, na.rm = TRUE),
         meankafdm_sub = mean(logkafdm, na.rm = TRUE), 
         stdvkafdm_sub = sd(logkafdm, na.rm = TRUE),
         meankicc_sub = mean(logkicc, na.rm = TRUE),
         stdvkicc_sub = sd(logkicc, na.rm = TRUE),
         meanktss_sub = mean(logktss, na.rm = TRUE),
         stdvktss_sub = sd(logktss, na.rm = TRUE)) %>%
  ungroup()%>%
  arrange(meank_sub) %>%  
  mutate(kspeed = case_when(meank_sub <= log10(0.005) ~ "kmin",
                            meank_sub > log10(0.005) & meank_sub < log10(0.02) ~ "slow",
                            meank_sub >= log10(0.02) & meank_sub < log10(0.1) ~ "medium",
                            meank_sub >= log10(0.1) ~ "fast")) %>%
  left_join(df_riverPara, by = "loc_site") %>%
  mutate(site = ifelse(grepl("^up", site), "up", "dw")) %>%
  left_join(df_ps, by = c("site", "x4lettercode")) %>%
  mutate(degr_class = case_when(substance %in% cata ~ "cata",
                                substance %in% cometa ~ "cometa",
                                substance %in% persistent ~ "persistent",
                                substance %in% highinf ~ "highinf",
                                TRUE ~ NA))

###############################################################################
# plotting in heatmap
# Customize the color palette
colors_kmean <- colorRamp2(c(min(data$meank_sub, na.rm=TRUE), max(data$meank_sub, na.rm=TRUE)),c("grey100", "#003000"))
colors_k <- colorRamp2(c(min(data$logk, na.rm=TRUE), max(data$logk, na.rm=TRUE)),c("grey100", "#003000"))
colors_env <- colorRamp2(c(0, 1), c("cornsilk", "dodgerblue3"))
colors_sub_class = c("Antibiotic" = "#9370DB", "Pharmaceutical" = "#4B0082", "Artificial sweetener" = "#87CEEB", 
                    "Household and Industrial Chemical" = "#3333FF", "Fungicide" = "#9ACD32", 
                    "Insecticide" = "#228B22", "Herbicide" = "#669966")
colors_updw = c("up" = "#008000", "dw" = "#3A5FCD")
colors_sd = colorRamp2(c(min(data$stdvk_sub, na.rm=TRUE), max(data$stdvk_sub, na.rm=TRUE)),c("grey100", "darkred"))
colors_ps = colorRamp2(c(min(data$mean_ps_site_norm, na.rm=TRUE), max(data$mean_ps_site_norm, na.rm=TRUE)), c("grey100", "aquamarine4"))
colors_temp = colorRamp2(c(min(data$Temperature, na.rm=TRUE), max(data$Temperature, na.rm=TRUE)), c("grey100", "#00008B"))
colors_loc = c("air" = "#0000FF",  "oli" = "#00008B",  "bir" = "#43CD80", "kol" = "#2E8B57", "alt" = "#FFA500",  "ehr" = "#CD6600")
site_sym = c("up1" = 3, "up2" = 7, "dw1" = 16, "dw2" = 17, "dw3" = 15)
colors_site <- c("up1_r1" = "#7FFF00", "up1_r2" = "#4CAF3C",
                 "up2_r1" = "#2E8B57", "up2_r2" = "#20543A",
                 "dw1_r1" = "#000080","dw1_r2" = "#2B3E8F",
                 "dw2_r1" = "#005073", "dw2_r2" =  "#104E8B",
                 "dw3_r1" = "#008078", "dw3_r2" = "#48D1CC")
colors_singlesite <- c("up1" = "#7FFF00", "up2" = "#2E8B57", "dw1" = "#000080", "dw2" = "#104E8B", "dw3" = "#48D1CC")
colors_sorption = c("TRUE" = "#000080", "FALSE" = "grey100")
colors_degr_class  <- c("highinf" = "#00008B", "cometa" = "green4", "cata" = "#FF1493", "persistent" = "grey")

#annotation for substances
ann_sub_colors <- list(meank_sub = colors_kmean,
                       stdvk_sub = colors_sd,
                       substance_class = colors_sub_class,
                       ps_up = colors_ps,
                       ps_dw = colors_ps,
                       sorption_flag = colors_sorption,
                       degr_class = colors_degr_class) 

#annoations for locbatch
ann_locbatch_colors <- list(location = colors_loc,
                            site = colors_updw,
                            Temperature = colors_temp,
                            pH = colors_env,
                            chl_a = colors_env,
                            afdm_0 = colors_env,
                            icc_0 = colors_env,
                            "Conductivity [µS/cm 20°C]" = colors_env,
                            "Alkalinity [mmol/L]" = colors_env,
                            "Total Hardness [mmol/L]" = colors_env,
                            "Sodium [mg/L]" = colors_env,
                            "Magnesium [mg/L]" = colors_env,
                            "Calcium [mg/L]" = colors_env,
                            "Potassium [mg/L]" = colors_env,
                            "Sulfate [mg/L]" = colors_env,
                            "Nitrate [mg/L]" = colors_env,
                            "Ammonium [µg/L]" = colors_env, 
                            "Nitrite [µg/L]" = colors_env,
                            "TP [µg/L]" = colors_env,
                            "Silicic Acid [mg/L]" = colors_env,
                            "DOC [mg/L]" = colors_env,
                            "TOC [mg/L]" = colors_env,
                            "TN [mg/L]" = colors_env,
                            "TIC [mg/L]" = colors_env)
  
df_annot_sub <- data %>%
  select(x4lettercode, meank_sub, stdvk_sub, substance_class, mean_ps_site_norm, site, sorption_flag) %>% #add: if degr_class should be annotated: , degr_class) %>%
  distinct() %>%
  mutate(site = glue("ps_{site}")) %>%
  pivot_wider(id_cols = c(x4lettercode, meank_sub, stdvk_sub, substance_class, sorption_flag), #add: if degr_class should be annotated: , degr_class),
              names_from = site, values_from = mean_ps_site_norm) %>%
  to_annot_sub()

df_annot_locbatch <- data %>%
  select(c("loc_batch", "location", "site", "Temperature", "pH", "chl_a", "afdm_0", "icc_0", "Conductivity [µS/cm 20°C]", "Alkalinity [mmol/L]", "Total Hardness [mmol/L]", "Sodium [mg/L]",
         "Magnesium [mg/L]","Calcium [mg/L]", "Potassium [mg/L]", "Sulfate [mg/L]", "Nitrate [mg/L]", "Ammonium [µg/L]",  
         "Nitrite [µg/L]", "TP [µg/L]", "Silicic Acid [mg/L]", "DOC [mg/L]", "TOC [mg/L]", "TN [mg/L]", "TIC [mg/L]")) %>%
  mutate(across(c(chl_a, afdm_0, icc_0), ~ round(. / max(., na.rm = TRUE), 2))) %>%
  distinct()  %>%
  to_annot_locbatch() 

df_cluster <- data %>%
  mutate(logk = log10(k)) %>% 
  select(x4lettercode, loc_batch, logk) %>%
  pivot_wider(id_cols = c(x4lettercode), names_from = loc_batch, values_from = logk) %>%
  to_matrix() 

# Cluster rows
row_dist <- dist(df_cluster, method = "euclidean")
row_clust <- hclust(row_dist, method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_clust))
row_dend <- as.dendrogram(row_clust)

# Cluster columns
col_dist <- dist(t(df_cluster), method = "euclidean")
col_clust <- hclust(col_dist, method = "complete")
col_order <- order.dendrogram(as.dendrogram(col_clust))
col_dend <- as.dendrogram(col_clust)

# Reorder matrix
df_cluster_ordered <- df_cluster[row_order, col_order] 

top_ha <- HeatmapAnnotation(df = df_annot_locbatch, 
                            col = ann_locbatch_colors,
                            na_col = "white", 
                            annotation_name_gp = gpar(fontsize = 10))

row_ha <- rowAnnotation(df = as.data.frame(df_annot_sub),
                        col = ann_sub_colors,
                        na_col = "white", 
                        annotation_name_gp = gpar(fontsize = 10))

ht <-Heatmap(df_cluster,
            col = colors_k,
            top_annotation = top_ha,
            left_annotation = row_ha,
            cluster_rows = row_dend,
            cluster_columns = col_dend,
            show_row_names = TRUE,
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 12),
            heatmap_legend_param = list(title = "k"),
            na_col = "grey90")

ht

png("plots/heatmap.png", width = 20, height = 25, units = "in", res = 300)
draw(ht)
dev.off()

######################################################################################
#PCA of logk
PCAdf_logk <- data %>%
  #filter(meank_sub > -2) %>% #enable, if you want to filter out persistent substances: meank_sub -2 (n=42)
  mutate(logk = log10(k)) %>%
  select(x4lettercode, loc_batch, logk) %>%
  pivot_wider(id_cols = c(x4lettercode), names_from = loc_batch, values_from = logk) %>%
  rowwise() %>%
  mutate(across(-x4lettercode, ~ ifelse(is.na(.), median(c_across(-x4lettercode), na.rm = TRUE), .))) %>% #replaces NA by the median of each the row (median of substance)
  ungroup()

PCAmatrix_logk <- PCAdf_logk %>%
  to_matrix() %>%
  t() %>%
  { .[, apply(., 2, sd, na.rm = TRUE) > 0]}

pca_result <- prcomp(PCAmatrix_logk, scale. = TRUE)  

pca_data <- as.data.frame(pca_result$x)
pca_data$Substance <- rownames(pca_data)

pca_variance_ratio <- pca_result$sdev^2 / sum(pca_result$sdev^2)  

pc1_var <- round(pca_variance_ratio[1] * 100, 1)
pc2_var <- round(pca_variance_ratio[2] * 100, 1)

df_meta <- PCAmatrix_logk %>%
  as.data.frame() %>%
  rownames_to_column("loc_batch") %>%
  select(loc_batch) %>%
  mutate(replicate = ifelse(str_detect(loc_batch, "r1"), "r1", "r2"),
         location  = word(loc_batch, 1, sep = "_"),
         site = word(loc_batch, 2, sep = "_"))

pca_data <- as.data.frame(pca_result$x) %>%
  rownames_to_column("loc_batch") %>%
  left_join(df_meta, by = "loc_batch")

ggplot(pca_data, aes(x = PC1, y = PC2, color = location, shape = site)) +
  geom_point(size = 3) +
  geom_point(data = subset(pca_data, replicate == "r2"),
             aes(x = PC1, y = PC2, color = location),
             shape = 1, size = 4.5, stroke = 1) +
  scale_color_manual(values = colors_loc) +
  scale_shape_manual(values = site_sym) +
  theme_minimal() +
  labs(title = "PCA of Rate Constants logk",
       x = paste0("PC1 - ", round(pc1_var, 1), "% variance"),
       y = paste0("PC2 - ", round(pc2_var, 1), "% variance"),
       color = "Location",
       shape = "Site") +
  theme(legend.position = "top")

ggsave("plots/PCA_all_MPs.png", plot = last_plot(), width = 10, height = 10, dpi = 300)

################################################################################
#boxplot
kmean_colors <- scales::seq_gradient_pal("white", "#003000")(seq(0, 1, length.out = 100))
var_logk_colors <- scales::seq_gradient_pal("grey100", "darkred")(seq(0,1, length.out = 100))
degr_class_col  <- c("highinf" = "#00008B", "cometa" = "green4", "cata" = "#FF1493", "persistent" = "grey")

data_box <- data %>%
  mutate(substance_code = paste(substance, x4lettercode, sep = "_")) %>%
  mutate(substance_code = reorder(substance_code, -meank_sub)) %>%
  mutate(degr_class = case_when(substance %in% cata ~ "cata",
                                substance %in% cometa ~ "cometa",
                                substance %in% persistent ~ "persistent",
                                substance %in% highinf ~ "highinf",
                                TRUE ~ NA))

ggplot(data_box, aes(x = substance_code, y = logk)) +
  geom_boxplot(fill = "grey", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = location), width = 0.2, alpha = 1, size = 1.5) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "white") + 
  scale_color_manual(values = colors_loc) +
  
  geom_tile(aes(x = substance_code, y = -2.5, fill = stdvk_sub), 
            inherit.aes = FALSE, width = 0.8, height = 0.1) +
  scale_fill_gradientn(name = "variance_logk", colors = var_logk_colors, na.value = "grey") +
  
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = substance_code, y = -2.7, fill = meank_sub),
            inherit.aes = FALSE, width = 0.8, height = 0.1) +
  scale_fill_gradientn(name = "mean_logk", colors = kmean_colors, na.value = "grey") +
  
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = substance_code, y = -2.9, fill = degr_class),
            inherit.aes = FALSE, width = 0.8, height = 0.1) +
  scale_fill_manual(name = "Degradation class", values = degr_class_col, na.value = "grey80") +
  
  labs(x = "Substance", y = "logk", title = "violinplot of logk by substance (mean_k ordered)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("plots/boxplot.png", plot = last_plot(), width = 13.3, height = 8.0, units = "in", dpi = 150)

#kick out persistent subst
data_box2 <- data %>%
  mutate(substance_code = paste(substance, x4lettercode, sep = "_")) %>%
  mutate(x4lettercode = reorder(x4lettercode, -meank_sub)) %>%
  mutate(degr_class = case_when(substance %in% cata ~ "cata",
                                substance %in% cometa ~ "cometa",
                                substance %in% persistent ~ "persistent",
                                substance %in% highinf ~ "highinf",
                                TRUE ~ NA)) %>%
  filter(!(meank_sub < -2.3))

pos_up <- which(levels(data_box2$x4lettercode) == 
                  unique(data_box2$x4lettercode[data_box2$substance == "Fenhexamid"]))
pos_down <- which(levels(data_box2$x4lettercode) == 
                    unique(data_box2$x4lettercode[data_box2$substance == "Icaridin"]))
pos_mid <- mean(c(pos_up, pos_down), na.rm = TRUE)

#changed for presentation (coloring of degradation and substance_code for x4lettercode)
ggplot(data_box2, aes(x = x4lettercode, y = logk)) +
  geom_boxplot(fill = "grey", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = location), width = 0.2, alpha = 1, size = 1.5) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "white") + 
  scale_color_manual(values = colors_loc) +
  
  geom_vline(xintercept = pos_mid, color = "black", linetype = "dashed", size = 1) +
  
  geom_tile(aes(x = x4lettercode, y = -2.7, fill = stdvk_sub), 
            inherit.aes = FALSE, width = 0.8, height = 0.1) +
  scale_fill_gradientn(name = expression(variance_log[10](k[deg])), colors = var_logk_colors, na.value = "grey") +
  
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = x4lettercode, y = -2.5, fill = meank_sub),
            inherit.aes = FALSE, width = 0.8, height = 0.1) +
  scale_fill_gradientn(name = expression(mean_log[10](k[deg])), colors = kmean_colors, na.value = "grey") +
  
  ggnewscale::new_scale_fill() +
  geom_tile(aes(x = x4lettercode, y = -2.9, fill = degr_class),
            inherit.aes = FALSE, width = 0.8, height = 0.1) +
  scale_fill_manual(name = "Degradation class", values = degr_class_col, na.value = "grey80") +
  
  labs(x = "Substance", y = expression(log[10](k[deg]))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("plots/boxplot_persout.png", plot = last_plot(), width = 14, height = 8.0, units = "in", dpi = 150)

#################################################################################
#violinplot for variance in biomass-corrected k
data_biovar <- data %>%
  select(substance, loc_batch, logk, logkafdm, logkicc, logktss, logkod) %>%
  group_by(substance) %>%
  mutate(var_logk = sd(logk, na.rm = TRUE),
         var_logkVSS = sd(logkafdm, na.rm = TRUE),
         var_logkTSS = sd(logktss, na.rm = TRUE),
         var_logkICC = sd(logkicc, na.rm = TRUE),
         var_logkOD = sd(logkod, na.rm = TRUE)) %>%
  slice(1) %>% 
  ungroup() %>%
  select(substance, var_logk, var_logkTSS, var_logkVSS, var_logkICC, var_logkOD) %>%
  pivot_longer(cols = -substance, names_to = "variance_type", values_to = "variance_value") %>%
  mutate(variance_type = str_remove(variance_type, "^var_"),
         variance_type = factor(variance_type, levels = c("logk", "logkTSS", "logkVSS", "logkICC", "logkOD")))

ggplot(data_biovar, aes(x = variance_type, y = variance_value, fill = variance_type)) +
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE) +  
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA) +
  geom_jitter(color = "black", width = 0.1, size = 1, alpha = 0.7, shape = 1) +  
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "white") +  
  labs(x = "Biomass Normalization", y = "Variance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave("plots/biovar.png", plot = last_plot(), width = 12, height = 8, units = "in", dpi = 150)

##################################################
#variances only for degraded substances
data_biovar_deg <- data %>%
  filter(!(meank_sub < -2)) %>%
  select(substance, loc_batch, logk, logkafdm, logkicc, logktss, logkod) %>%
  group_by(substance) %>%
  mutate(var_logk = sd(logk, na.rm = TRUE),
         var_logkVSS = sd(logkafdm, na.rm = TRUE),
         var_logkTSS = sd(logktss, na.rm = TRUE),
         var_logkICC = sd(logkicc, na.rm = TRUE),
         var_logkOD = sd(logkod, na.rm = TRUE)) %>%
  slice(1) %>% 
  ungroup() %>%
  select(substance, var_logk, var_logkTSS, var_logkVSS, var_logkICC, var_logkOD) %>%
  pivot_longer(cols = -substance, names_to = "variance_type", values_to = "variance_value") %>%
  mutate(variance_type = str_remove(variance_type, "^var_"),
         variance_type = factor(variance_type, levels = c("logk", "logkTSS", "logkVSS", "logkICC", "logkOD")))

ggplot(data_biovar_deg, aes(x = variance_type, y = variance_value, fill = variance_type)) +
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE) +  
  geom_boxplot(color = "black", fill = NA, width = 0.2, outlier.shape = NA) +
  geom_jitter(color = "black", width = 0.1, size = 1, alpha = 0.7, shape = 1) +  
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "white") +  
  labs(x = "Biomass Normalization", y = "Variance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave("plots/biovar_deg.png", plot = last_plot(), width = 12, height = 8, units = "in", dpi = 150)

