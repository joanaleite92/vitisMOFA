library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)
library(ggrepel)
library(data.table)
library(fuzzyjoin)
library(purrr)

# --- Assuming df_merged is loaded here ---
load("df_final5_leaf_metbots_noOut_starchCorrect.RData")

df_duplicatas_completas <- result2_noOut_starch2%>%
  # 1. Agrupar por todas as colunas, exceto 'Laboratory'
  group_by(across(-Laboratory)) %>%
  
  # 2. Contar o número de ocorrências em cada grupo (n())
  # 3. Filtrar apenas os grupos onde a contagem é maior que 1
  filter(n() > 1) %>%
  
  # 4. Desagrupar e organizar para fácil inspeção (opcional)
  ungroup() %>%
  arrange(across(-Laboratory))

df_leaf_metabs2023$FW_sample <- factor(
  df_leaf_metabs2023$FW_sample,
  levels = c("S1", "S2", "S3","S4","S5","S6"),
  ordered = TRUE
)

df_leaf_metabs2023_temp <- df_leaf_metabs2023 %>%
  rename(
    Carotenoids__mg_gFM_new = Carotenoids__mg_gFM_,
    Chlorophyll_a__mg_gFM_new = Chlorophyll_a__mg_gFM_,
    Chlorophyll_b__mg_gFM_new = Chlorophyll_b__mg_gFM_,
    Phenols__g_GAE_mL_new = Phenols__g_GAE_mL_,
    ROS_O2__ABS_gFM_new = ROS_O2__ABS_gFM_,
    Starch__mg_gFM_new = Starch__mg_gFM_,
    Sugar__mg_gFM_new = Sugar__mg_gFM_,
    Anthocyanins__mg_gFM_new = Anthocyanins__mg_gFM_
  )
# --- Passo 2: Juntar os Dataframes ---
# Usamos LEFT JOIN para manter TODAS as linhas do DF_Principal.
df_leaf_joined <- result2_noOut_starch2 %>%
  left_join(
    df_leaf_metabs2023_temp,
    by = c("Code_field", "FW_sample") # Chave composta
  )

# --- Passo 3: Preencher os NAs com COALESCE ---
# coalesce(coluna_original, coluna_nova) substitui o NA na coluna_original
# pelo valor correspondente na coluna_nova.

df_leaf_joined_complete <- df_leaf_joined %>%
  mutate(
    # Para a Coluna_X: use o valor de Coluna_X; se for NA, use Coluna_X_novo
    Carotenoids__mg_gFM_ = coalesce(Carotenoids__mg_gFM_, Carotenoids__mg_gFM_new),
    Chlorophyll_a__mg_gFM_ = coalesce(Chlorophyll_a__mg_gFM_, Chlorophyll_a__mg_gFM_new),
    Chlorophyll_b__mg_gFM_ = coalesce(Chlorophyll_b__mg_gFM_,Chlorophyll_b__mg_gFM_new),
    Phenols__g_GAE_mL_ = coalesce(Phenols__g_GAE_mL_,Phenols__g_GAE_mL_new),
    ROS_O2__ABS_gFM_ = coalesce(ROS_O2__ABS_gFM_,ROS_O2__ABS_gFM_new),
    Starch__mg_gFM_ = coalesce(Starch__mg_gFM_,Starch__mg_gFM_new),
    Sugar__mg_gFM_ = coalesce(Sugar__mg_gFM_,Sugar__mg_gFM_new),
    Anthocyanins__mg_gFM_ = coalesce(Anthocyanins__mg_gFM_,Anthocyanins__mg_gFM_new)
  ) %>%
  select(-Carotenoids__mg_gFM_new, -Chlorophyll_a__mg_gFM_new, -Chlorophyll_b__mg_gFM_new, -Phenols__g_GAE_mL_new,
         -ROS_O2__ABS_gFM_new,-Starch__mg_gFM_new,-Sugar__mg_gFM_new,-Anthocyanins__mg_gFM_new)

save(df_leaf_joined_complete, file = "df_final7_leaf_metbots_noOut_starchCorrect_withlags.RData")

result2_noOut_starch2 <- result2_noOut_starch2 %>%
  filter(Test_site != "Greenhouse Porto")

df_leaf_joined_complete <- df_leaf_joined_complete%>%
  filter(!Test_site %in% c("Greenhouse Porto", "Quinta de Vale de Cavalos"))
# 1. Ensure the Timepoint column is correctly ordered (e.g., using date or yday)
# If Timepoint is a Date object, this sort is essential
result2_noOut_starch2 <- result2_noOut_starch2 %>%
  arrange(Code_field, year, FW_sample)

df_leaf_joined_complete <- df_leaf_joined_complete %>%
  arrange(Code_field, year, FW_sample)

# 2. Identify the columns you want to lag
metabolite_cols <- c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                     "Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                     "Sugar__mg_gFM_","Anthocyanins__mg_gFM_")
gene_cols <- c("RCA (Cq)", "LBCY (Cq)", "CHLG (Cq)")
# Combine them into one list for the lagging function
omics_cols <- c(metabolite_cols, gene_cols)

create_lagged_columns_safe <- function(df, columns_to_lag) {
  
  df_lagged <- df %>%
    group_by(Code_field, year) %>%
    mutate(
      # Use across and all_of to apply lag to all columns in the list
      across(
        .cols = all_of(columns_to_lag),
        .fns = ~ lag(.x, n = 1, default = NA),
        .names = "{.col}_lag7d"
      )
    ) %>%
    ungroup()
  
  return(df_lagged)
}

# Run the function
df_lagged_leaf <- create_lagged_columns_safe(result2_noOut_starch2, omics_cols)
df_leaf_joined_complete_lagged <- create_lagged_columns_safe(df_leaf_joined_complete, omics_cols)

save(df_lagged_leaf, file = "df_final6_leaf_metbots_noOut_starchCorrect_withlags.RData")
save(df_leaf_joined_complete_lagged, file = "df_final7_leaf_metbots_noOut_starchCorrect_withlags.RData")

  # Keep the first row found for each unique combination of columns.
  distinct(
    across(-Laboratory), 
    .keep_all = TRUE
  )
########################################
load("df_final2_fruit_metbots.RData")
load("df_final4_leaf_metbots_noOut_starchCorrect.RData")
df_final <- load("df_final_leaf_metbots.RData")

df_climate <- read_excel("SOUTELO.xlsx", sheet = 2)

str(df_climate)
#str(df_dedup)

df_climate$Data <- ymd(df_climate$Data)
df_climate <- df_climate  %>%
  rename(Date = Data)
str(df_climate)

df_climate_21 <- read_excel("Clima_SouteloDoDouro_21_22.xlsx",sheet = 1)
df_climate_21$Data <- ymd(df_climate_21$Data)
df_climate_21 <- df_climate_21  %>%
  rename(Date = Data)
df_climate_22 <- read_excel("Clima_SouteloDoDouro_21_22.xlsx",sheet = 2)
df_climate_22$Data <- ymd(df_climate_22$Data)
df_climate_22 <- df_climate_22  %>%
  rename(Date = Data)
df_climate_23 <- read_excel("Soutelo2023.xlsx",sheet = 1)
df_climate_23$Date <- ymd(df_climate_23$Date)

df_climate_23 <- df_climate_23 %>%
  # Step 2a: Convert the DateTimePrecipitacao column to a proper datetime object
  # and extract just the date part for grouping
  mutate(
    Date = as.Date(Date)
  ) %>%
  # Step 2b: Group the data by the extracted Date
  group_by(Date) %>%
  # Step 2c: Calculate the mean for the desired columns
  dplyr::summarise(
    # Use 'mean()' to calculate the daily average for each column
    `P (mm)` = sum(`Precipitacao (mm)`, na.rm = TRUE),
    `HRmed (%)` = mean(`Humidade Relativa (% RH)`, na.rm = TRUE),
    `HRmax (%)` = max(`Humidade Relativa (% RH)`, na.rm = TRUE),
    `HRmin (%)` = min(`Humidade Relativa (% RH)`, na.rm = TRUE),
    `Tmed (ºC)` = mean(`Temperatura (°C)`, na.rm = TRUE),
    `Tmax (ºC)` = max(`Temperatura (°C)`, na.rm = TRUE),
    `Tmin (ºC)` = min(`Temperatura (°C)`, na.rm = TRUE),
    `Radmed (W/m2)` = mean(`Radiacao (W/m²)`, na.rm = TRUE),
    `Radmax (W/m2)` = max(`Radiacao (W/m²)`, na.rm = TRUE),
    `Radmin (W/m2)` = min(`Radiacao (W/m²)`, na.rm = TRUE),
    `Windmed (km/h)` = mean(`Velocidade Vento (km/h)`, na.rm = TRUE),
    `Windmax (km/h)` = max(`Velocidade Vento (km/h)`, na.rm = TRUE),
    .groups = 'drop' # Ungroup the data after summarizing
  )

df_climate_all <- bind_rows(df_climate_21,df_climate_22,df_climate_23,df_climate)
str(df_climate_all)

df_dedup$Date <- ymd(df_dedup$Date)
str(df_dedup)

df_merged_fruit_metbots$Date <- ymd(df_merged_fruit_metbots$Date)
str(df_merged_fruit_metbots)

# result <- df_dedup %>%
#   rowwise() %>%
#   mutate(
#     mean_P     = mean(df_climate_all$`P (mm)`   [df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_HRmed = mean(df_climate_all$`HRmed (%)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_Tmed  = mean(df_climate_all$`Tmed (ºC)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_Tmin  = mean(df_climate_all$`Tmin (ºC)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_Tmax  = mean(df_climate_all$`Tmax (ºC)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE)
#   ) %>%
#   ungroup()

########## Compute lags and rolling means for several sizes #####
get_climate_features <- function(climate_df, sample_date, var,
                                 roll_windows = c(3,7,14),
                                 lags = c(1,7),
                                 roll_fun = mean) {
  # Ensure Date column is Date class
  if(!inherits(climate_df$Date, "Date")) {
    stop("climate_df$Data must be of class Date")
  }
  
  results <- list()
  
  # --- Rolling windows ---
  for (win in roll_windows) {
    start_date <- sample_date - win
    end_date   <- sample_date - 1
    vals <- climate_df[[var]][climate_df$Date >= start_date & climate_df$Date <= end_date]
    
    feat_name <- paste0(var, "_roll", win, "d_", fun_name(roll_fun))
    results[[feat_name]] <- if (length(vals) == 0) NA_real_ else roll_fun(vals, na.rm = TRUE)
  }
  
  # --- Pure lags (raw values from exact past days) ---
  for (lag in lags) {
    lag_date <- sample_date - lag
    val <- climate_df[[var]][climate_df$Date == lag_date]
    
    feat_name <- paste0(var, "_lag", lag, "d")
    results[[feat_name]] <- if (length(val) == 0) NA_real_ else val
  }
  
  return(as.data.frame(results))
}

# Helper for naming
fun_name <- function(f) {
  if (identical(f, mean)) return("mean")
  if (identical(f, sum))  return("sum")
  if (identical(f, max))  return("max")
  if (identical(f, min))  return("min")
  return("fun")
}

all_feats <- lapply(df_climate_all$Date, function(d)
  cbind(Date = d,
        get_climate_features(df_climate_all, d, var="P (mm)",
                             roll_windows=c(3,7,14), lags=c(1,3,7), roll_fun=mean))
) %>% bind_rows()

all_feats_Tmax <- lapply(df_climate_all$Date, function(d)
  cbind(Date = d,
        get_climate_features(df_climate_all, d, var="Tmax (ºC)",
                             roll_windows=c(3,7,14), lags=c(1,3,7), roll_fun=mean))
) %>% bind_rows()

all_feats_Tmin <- lapply(df_climate_all$Date, function(d)
  cbind(Date = d,
        get_climate_features(df_climate_all, d, var="Tmin (ºC)",
                             roll_windows=c(3,7,14), lags=c(1,3,7), roll_fun=mean))
) %>% bind_rows()

all_feats_Tmed <- lapply(df_climate_all$Date, function(d)
  cbind(Date = d,
        get_climate_features(df_climate_all, d, var="Tmed (ºC)",
                             roll_windows=c(3,7,14), lags=c(1,3,7), roll_fun=mean))
) %>% bind_rows()

all_feats_HRmin <- lapply(df_climate_all$Date, function(d)
  cbind(Date = d,
        get_climate_features(df_climate_all, d, var="HRmin (%)",
                             roll_windows=c(3,7,14), lags=c(1,3,7), roll_fun=mean))
) %>% bind_rows()

all_feats_HRmed <- lapply(df_climate_all$Date, function(d)
  cbind(Date = d,
        get_climate_features(df_climate_all, d, var="HRmed (%)",
                             roll_windows=c(3,7,14), lags=c(1,3,7), roll_fun=mean))
) %>% bind_rows()


# Merge all feature dataframes by Date
all_feats_combined <- reduce(list(all_feats, all_feats_Tmax, all_feats_Tmin, all_feats_Tmed, all_feats_HRmin,all_feats_HRmed),
                             function(x, y) left_join(x, y, by = "Date"))


result <- df_dedup %>%
  left_join(all_feats_combined, by = c("Date" = "Date"))

result2_noOut_starch <- result2_noOut_starch %>%
  left_join(all_feats_combined, by = c("Date" = "Date"))
#result <- df_dedup %>%
  #left_join(all_feats, by = c("Date" = "Date"))

result2_noOut_starch2 <- result2_noOut_starch2 %>%
  left_join(all_feats_combined, by = c("Date" = "Date"))

result2_noOut_starch2 <- result2_noOut_starch2 %>% rename(Water_P = Water_P.x)
result2_noOut_starch2 <- result2_noOut_starch2 %>%
  # 2. Drop all columns that end with .y
  # The 'matches' function allows you to use a regular expression
  select(-matches("\\.x$")) %>%
  
  # 3. Rename all columns that end with .x by removing the suffix
  # rename_with applies a function to column names that match a criteria
  rename_with(
    .fn = ~ stringr::str_replace(., "\\.y$", ""), # Function: replace '.x' at the end of the string
    .cols = matches("\\.y$")                      # Criteria: only apply to columns ending in '.x'
  )
save(result2_noOut_starch2, file = "df_final5_leaf_metbots_noOut_starchCorrect.RData")

result_fruit <- result_fruit %>% left_join(all_feats_combined,by = c("Date" = "Date"))
save(result_fruit, file = "df_final2_fruit_metbots.RData")

#result2_noOut_starch <- result2_noOut_starch %>% mutate(Code_field = factor(Code_field))
result2_noOut_starch <- result2_noOut_starch %>% mutate(Test_site = factor(Test_site))
result2_noOut_starch <- result2_noOut_starch %>% mutate(Cultivar = factor(Cultivar))
result2_noOut_starch <- result2_noOut_starch %>% mutate(year = factor(year))
result2_noOut_starch <- result2_noOut_starch %>% mutate(Laboratory = factor(Laboratory))
result2_noOut_starch$Irrigation_bin <- ifelse(result2_noOut_starch$Irrigation == "I", 1, 0)
#result2_noOut_starch$irrigation_bin <- scale(result2_noOut_starch$irrigation_bin)   # optional: z-score if mixing with other covariates
save(result2_noOut_starch, file = "df_final3_leaf_metbots_noOut_starchCorrect.RData")

result2_noOut_starch2 <- result2_noOut_starch2 %>% mutate(Test_site = factor(Test_site))
result2_noOut_starch2 <- result2_noOut_starch2 %>% mutate(Cultivar = factor(Cultivar))
result2_noOut_starch2 <- result2_noOut_starch2 %>% mutate(year = factor(year))
result2_noOut_starch2 <- result2_noOut_starch2 %>% mutate(Laboratory = factor(Laboratory))
#########################################################################################
save(result2_noOut_starch2, file = "df_final4_leaf_metbots_noOut_starchCorrect.RData")

# result_fruit <- df_merged_fruit_metbots %>%
#   rowwise() %>%
#   mutate(
#     mean_P     = mean(df_climate_all$`P (mm)`   [df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_HRmed = mean(df_climate_all$`HRmed (%)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_Tmed  = mean(df_climate_all$`Tmed (ºC)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_Tmin  = mean(df_climate_all$`Tmin (ºC)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE),
#     mean_Tmax  = mean(df_climate_all$`Tmax (ºC)`[df_climate_all$Data %in% (Date - 7):(Date - 1)], na.rm = TRUE)
#   ) %>%
#   ungroup()

#change irrigation variable from categorical to binary: NI -> 0; I -> 1
result$Irrigation <- factor(ifelse(result$Irrigation == "NI", 0,
                                   ifelse(result$Irrigation == "I", 1, NA)),
                            levels = c(0,1),
                            labels = c("NI","I"))

result_fruit$Irrigation <- factor(ifelse(result_fruit$Irrigation == "NI", 0,
                                   ifelse(result_fruit$Irrigation == "I", 1, NA)),
                            levels = c(0,1),
                            labels = c("NI","I"))

result <- result %>%
  mutate(
    year = factor(year(Date)),
    month = month(Date),
    day = day(Date),
    weekday = wday(Date, label = TRUE),
    yday = yday(Date),  # day of the year
    month_sin = sin(2*pi*month/12),
    month_cos = cos(2*pi*month/12),
    day_sin = sin(2 * pi * yday / 365),
    day_cos = cos(2 * pi * yday / 365)
  )

result_fruit <- result_fruit %>%
  mutate(
    year = factor(year(Date)),
    month = month(Date),
    day = day(Date),
    weekday = wday(Date, label = TRUE),
    yday = yday(Date),  # day of the year
    month_sin = sin(2*pi*month/12),
    month_cos = cos(2*pi*month/12),
    day_sin = sin(2 * pi * yday / 365),
    day_cos = cos(2 * pi * yday / 365)
  )

# convert FW_sample (timepoint) to ordered factor
result$FW_sample <- factor(result$FW_sample, levels = c("S1","S2","S3","S4","S5","S6"), ordered = TRUE)


# convert FW_sample (timepoint) to ordered factor
result_fruit$FW_sample <- factor(result_fruit$FW_sample, levels = c("S1","S2","S3","S4","S5","S6"), ordered = TRUE)

result <- result %>% mutate(ID = factor(ID))
result <- result %>% mutate(Code_field = factor(Code_field))
result <- result %>% mutate(Test_site = factor(Test_site))
result <- result %>% mutate(Cultivar = factor(Cultivar))

result_fruit <- result_fruit %>% mutate(ID = factor(ID))
result_fruit <- result_fruit %>% mutate(Code_field = factor(Code_field))
result_fruit <- result_fruit %>% mutate(Test_site = factor(Test_site))
result_fruit <- result_fruit %>% mutate(Cultivar = factor(Cultivar))
result_fruit$year <- as.factor(result_fruit$year)

result <- result %>% select(-Project, -Crop)
result_fruit <- result_fruit %>% select(-Project, -Crop)

# separate latitude and longitude into two separate columns
result <- result %>%
  separate(Coord, into = c("Longitude", "Latitude"), sep = ";") %>%
  mutate(across(c(Longitude, Latitude), as.numeric))

result_fruit <- result_fruit %>%
  separate(Coord, into = c("Longitude", "Latitude"), sep = ";") %>%
  mutate(across(c(Longitude, Latitude), as.numeric))

#Days since start of each season/year
result <- result %>%
  group_by(year) %>%
  mutate(days_since_year_start = as.numeric(date - min(date))) %>%
  ungroup()

result_fruit <- result_fruit %>%
  group_by(year) %>%
  mutate(days_since_year_start = as.numeric(Date - min(Date))) %>%
  ungroup()

save(result, file = "df_final2_leaf_metbots.RData")

load("PHB_LEAF_DATA.RData")
PHB_LEAF <- PHB_LEAF %>% mutate(ID = factor(ID))
PHB_LEAF <- PHB_LEAF %>% mutate(Code_field = factor(Code_field))
PHB_LEAF <- PHB_LEAF %>% mutate(Test_site = factor(Test_site))
PHB_LEAF <- PHB_LEAF %>% mutate(Cultivar = factor(Cultivar))

result2 <- result %>%
  left_join(PHB_LEAF %>% dplyr::select(ID, Code_field, Test_site, Cultivar,Water_P),
            by = c("ID",  "Code_field", "Test_site","Cultivar"))

##### Correct all pigment values by a factor of x10^-3
# Specify the columns to divide
cols_to_divide <- c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                    "Anthocyanins__mg_gFM_")

# Divide the selected columns by 1000
result2[cols_to_divide] <- result2[cols_to_divide] / 1000
save(result2, file = "df_final3_leaf_metbots.RData")

### detect and remove outliers
target_cols_out <- c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                              "Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                              "Sugar__mg_gFM_","Anthocyanins__mg_gFM_","MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                     "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
                     "LBCY (Cq)","CHLG (Cq)","GAPDH (Cq)")

# Function to flag outliers using IQR rule
remove_outliers <- function(x) {
  if (all(is.na(x))) return(x)  # handle all-NA groups safely
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower <- Q1 - 1.5 * IQR
  upper <- Q3 + 1.5 * IQR
  x[x < lower | x > upper] <- NA
  return(x)
}
# Apply per Test_site × Year × FW_sample
result2_noOut <- result2 %>%
  group_by(Test_site, year, FW_sample) %>%
  mutate(across(all_of(target_cols), remove_outliers)) %>%
  ungroup()

save(result2_noOut, file = "df_final3_leaf_metbots_noOut.RData")
####### FRUIT - join PHB ##########
#result_fruit <- result_fruit %>% select(-Code_field_temp) 

## Turn aoo values for Starch NA if year == 2023 && Test_site == Quinta dos Aciprestes
result2_noOut_starch <- result2_noOut %>%
  mutate(Starch__mg_gFM_ = ifelse(year == 2023 & Test_site == "Quinta dos Aciprestes",
                         NA, Starch__mg_gFM_))
sum(is.na(result2_noOut$Starch__mg_gFM_))
sum(is.na(result2_noOut_starch$Starch__mg_gFM_))

load("df_final3_leaf_metbots_noOut_starchCorrect.RData")


result2_noOut_starch2 <- fuzzy_left_join(
  result2_noOut_starch %>% mutate(Date_minus7 = Date - 8),
  PHB_LEAF,
  by = c("Code_field" = "Code_field",
         "Date_minus7" = "Date"),
  match_fun = list(`==`, function(x, y) abs(as.numeric(x - y)) <= 5)
) %>% rename(
  ID = ID.x,
  Code_field = Code_field.x,
  Test_site = Test_site.x,
  Cultivar = Cultivar.x,
  Date = Date.x,
  FW_sample = FW_sample.x,
  Irrigation = Irrigation.x,
  Water_P_lag7daytol = Water_P.y
) %>% select(-ID.y, -Code_field.y, -FW_sample.y,-Test_site.y, -Date.y, -Crop, -Cultivar.y, -Irrigation.y,
             -Coord, -Project, -join_key)
save(result2_noOut_starch2, file = "df_final4_leaf_metbots_noOut_starchCorrect.RData")


result2_starch <- result2 %>%
  mutate(Starch__mg_gFM_ = ifelse(year == 2023 & Test_site == "Quinta dos Aciprestes",
                                  NA, Starch__mg_gFM_))
sum(is.na(result2$Starch__mg_gFM_))
sum(is.na(result2_starch$Starch__mg_gFM_))
save(result2_starch, file = "df_final3_leaf_metbots_starchCorrect.RData")

meta_cols <- c("FW_sample","Test_site", "Cultivar","Irrigation", "Longitude", "Latitude", "Laboratory",
               "mean_P","mean_Tmed","mean_Tmax","mean_HRmed","Water_P",
              "month_sin","month_cos","day_sin","day_cos","yday")
              
meta_df_leaf_nout_starch <- result2_noOut_starch

# Create a temporary join key for the factor column (without altering original)

PHB_LEAF <- PHB_LEAF %>% mutate(join_key = as.character(Code_field))
result_fruit <- result_fruit %>% mutate(join_key = gsub("\\.[A-Za-z]$", "", as.character(Code_field)))


PHB_LEAF <- PHB_LEAF %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))


PHB_LEAF <- PHB_LEAF %>%
  mutate(year = year(Date),
         month = month(Date) ,
         day = yday(Date))
##### important #####
PHB_LEAF <- PHB_LEAF %>%
  mutate(Water_P = ifelse(Water_P == 0, NA, Water_P))

# Restrict to July–Sept
PHB_season <- PHB_LEAF %>%
  filter(month %in% c(7, 8, 9)) %>%
  group_by(Code_field, year, month) %>%
  summarise(
    mean_PHB = mean(Water_P, na.rm = TRUE),
    min_PHB  = min(Water_P, na.rm = TRUE),
    max_PHB  = max(Water_P, na.rm = TRUE),
    .groups = "drop"
  )


# Assume you already have phb_season from your summarise()
# Columns: PlantID, Year, Month, mean_PHB, min_PHB, max_PHB

# Convert Month to ordered factor so Jul → Aug → Sep are in order
PHB_season <- PHB_season %>%
  mutate(month = factor(month, levels = c(7, 8, 9)))

# Filter to July–September
phb_season_daily <- PHB_LEAF %>%
  filter(month(Date) %in% 7:9)


# Plot usando DOY no eixo x
ggplot(phb_season_daily, aes(x = day, y = Water_P, color = factor(year))) +
  geom_point(alpha = 0.5) +     # pontos diários
  #geom_line(alpha = 0.7) +      # linha conectando os pontos
  labs(
    title = "Daily dynamics BWP (July-Sep",
    x = "Day of the year",
    y = "BWP (MPa)",
    color = "Year"
  ) +
  theme_minimal(base_size = 14)

# Calcular média (e SD) por Ano × DOY
phb_daily_mean <- phb_season_daily %>%
  group_by(year, day) %>%
  summarise(
    mean_PHB = mean(Water_P, na.rm = TRUE),
    sd_PHB   = sd(Water_P, na.rm = TRUE),
    n        = n(),
    se_PHB   = sd_PHB / sqrt(n),   # erro padrão da média
    .groups = "drop"
  )

# Plot média diária com linhas
ggplot(phb_daily_mean, aes(x = day, y = mean_PHB, color = factor(year))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_PHB - se_PHB, ymax = mean_PHB + se_PHB),
                width = 0.3, alpha = 0.4) +
  labs(
    title = "Daily dynamics BWP (averages per year, July-September)",
    x = "Day of the year",
    y = "Average BWP (MPa)",
    color = "Year"
  ) +
  theme_minimal(base_size = 14)

# Aggregate across plants (per Year × Month), with SD for error bars
phb_summary <- PHB_season %>%
  group_by(year, month) %>%
  summarise(
    mean_PHB = mean(mean_PHB, na.rm = TRUE),
    sd_PHB   = sd(mean_PHB, na.rm = TRUE),
    .groups = "drop"
  )

# Plot with error bars
ggplot(phb_summary, aes(x = month, y = mean_PHB, group = year, color = factor(year))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_PHB - sd_PHB, ymax = mean_PHB + sd_PHB), width = 0.2) +
  labs(
    title = "Seasonal PHB dynamics (July–September)",
    x = "Month",
    y = "Mean PHB (MPa)",
    color = "Year"
  ) +
  theme_minimal(base_size = 14)


# aggregate across whole season (July–Sept)
phb_whole_season <- PHB_season %>%
  group_by(Code_field, year) %>%
  summarise(
    mean_PHB_season = mean(mean_PHB, na.rm = TRUE),
    min_PHB_season  = min(min_PHB, na.rm = TRUE),
    .groups = "drop"
  )


phb_summary2 <- phb_whole_season %>%
  group_by(year) %>%
  summarise(mean_PHB = mean(mean_PHB_season, na.rm = TRUE),
            sd_PHB   = sd(mean_PHB_season, na.rm = TRUE),
            .groups = "drop")

ggplot(phb_summary2, aes(x = year, y = mean_PHB)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_PHB - sd_PHB, ymax = mean_PHB + sd_PHB),
                width = 0.2) +
  labs(title = "Seasonal mean PHB across years (all fields)",
       x = "Year", y = "Mean PHB (MPa)") +
  theme_minimal(base_size = 14)

phb_whole_season_filtered <- phb_whole_season %>%
  filter(year %in% 2021:2024)

library(RColorBrewer)  # for nice color palettes

# Generate a unique color for each field
library(viridis)

# 1. Create large palette
fields <- unique(phb_whole_season_filtered$Code_field)
n_fields <- length(fields)
palette <- viridis(n_fields, option = "turbo")
names(palette) <- fields

# 2. Prepare labels (optional: top 20)
# Generate 334 distinct colors
palette <- viridis(n_fields, option = "turbo")
names(palette) <- fields

# Prepare labels for only top 20 final values
labels_df <- phb_whole_season_filtered %>%
  group_by(Code_field) %>%
  slice_max(year, n = 1) %>%
  slice_max(mean_PHB_season, n = 20)

# Plot
ggplot(phb_whole_season_filtered, aes(x = year, y = mean_PHB_season, color = Code_field, group = Code_field)) +
  geom_line(size = 0.5) +
  geom_point(size = 0.7) +
  geom_text_repel(
    data = labels_df,
    aes(label = Code_field),
    nudge_x = 0.2,
    direction = "y",
    hjust = 0,
    segment.color = NA,
    size = 2.5
  ) +
  scale_color_manual(values = palette) +
  labs(title = "Seasonal mean PHB across years (2021–2024)",
       x = "Year", y = "Mean PHB (MPa)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")



ggplot(phb_whole_season_filtered, aes(x = year, y = mean_PHB_season, color = Code_field, group = Code_field)) +
  geom_line() +
  geom_point() +
  labs(title = "Seasonal mean PHB across years (2021–2024)",
       x = "Year", y = "Mean PHB (MPa)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Identify which columns from PHB_LEAF are new (not overlapping with result_fruit, except join_key)
new_cols <- setdiff(names(PHB_LEAF), c(names(result_fruit), "join_key"))

# Perform the join by date - same day as the leaf
result_fruit2 <- result_fruit %>%
  left_join(select(PHB_LEAF, join_key, Date, all_of(new_cols)), by = c("join_key", "Date")) 
# %>% select(-join_key)  # drop helper column

####################################################
# 1) Aggregate PHB_LEAF by month

PHB_LEAF_monthly <- PHB_LEAF %>%
  mutate(year = as.numeric(as.character(year))) %>%  # convert factor year to numeric
  group_by(join_key, year, month) %>%                # group by join_key, year, month
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    phb_year = year,
    phb_month = month
  )
PHB_LEAF_monthly <- PHB_LEAF_monthly %>%
  select(-year, -month)

# 2) Prepare fruit dataset with previous month/year
result_fruit2 <- result_fruit2 %>%
  mutate(year = as.numeric(as.character(year))) %>%
  mutate(
    month_prev = ifelse(month == 1, 12, month - 1),
    year_prev  = ifelse(month == 1, year - 1, year)
  )

#3) Identify PHB_LEAF columns to join
new_cols <- setdiff(names(PHB_LEAF_monthly), c(names(result_fruit), "join_key", "year", "month"))

# Each month gets previous month average, July fruits get July PHB (proxy)
result_fruit_A <- result_fruit2 %>%
  mutate(
    # Use July PHB itself if month == 7
    month_use = ifelse(month == 7, month, month_prev),
    year_use  = ifelse(month == 7, year, year_prev)
  ) %>%
  left_join(
    select(PHB_LEAF_monthly, join_key, year, month, all_of(new_cols)),
    by = c("join_key", "year_use" = "year", "month_use" = "month")
  )
########################################
#Each month gets previous month average, July fruits get NA
result_fruit_B <- result_fruit2 %>%
  left_join(
    select(PHB_LEAF_monthly, join_key, year, month, all_of(new_cols)),
    by = c("join_key", "year_prev" = "year", "month_prev" = "month")
  ) %>%
  select(-join_key)

########################################
# Cumulative PHB up to harvest month

result_fruit_cumulative <- result_fruit2 %>%
  rowwise() %>% 
  mutate(
    across(all_of(new_cols),
           ~ mean(PHB_LEAF[[cur_column()]][
             PHB_LEAF$join_key == join_key &
               lubridate::year(PHB_LEAF$Date) == year &
               PHB_LEAF$Date < Date
           ], na.rm = TRUE)
    )
  ) %>%
  ungroup()
########################################
##Join tha point value at 7 days before harvesting

#Create shifted date column in fruit data
# Add column for 7 days before harvest
result_fruit2 <- result_fruit2 %>%
  mutate(Date_minus7 = Date - days(7))

#Join PHB data on the shifted date: Date_minus7 = harvest date − 7 days.

# Select only needed PHB columns
new_cols <- setdiff(names(PHB_LEAF), c("join_key", "Date"))

result_fruit2_7daybef <- result_fruit2 %>%
  left_join(
    PHB_LEAF %>% select(join_key, Date, all_of(new_cols)),
    by = c("join_key", "Date_minus7" = "Date")
  )

##########

# Nearest measurement within a tolerance (±X days) - Rolling join 

# Convert to data.table
setDT(result_fruit2)
setDT(PHB_LEAF)

# Make sure Date is Date
result_fruit2[, Date := as.Date(Date)]
PHB_LEAF[, Date := as.Date(Date)]

# Add target date: 7 days before harvest
result_fruit2[, Date_minus7 := Date - 7]

# Rolling join: get nearest PHB date to Date_minus7 (backward or nearest)
#roll = "nearest" finds the closest PHB measurement.
#You can also use roll = -7 to enforce backward-only (e.g., no “future” data leaks).
#If you want a maximum tolerance, e.g. only accept if within ±2 days:
result_fruit2_7daysbefappr <- PHB_LEAF[result_fruit2,
                                  on = .(join_key, Date = Date_minus7),
                                  roll = "nearest"]

#OR use fuzzyjoin
library(fuzzyjoin)

#tried tolerance 2,3 and 4 and no difference
result_fruit2 <- fuzzy_left_join(
  result_fruit2 %>% mutate(Date_minus7 = Date - 8),
  PHB_LEAF,
  by = c("join_key" = "join_key",
         "Date_minus7" = "Date"),
  match_fun = list(`==`, function(x, y) abs(as.numeric(x - y)) <= 5)
)

names(result_fruit2) <- gsub("\\.x$", "", names(result_fruit2))
result_fruit2 <- result_fruit2 %>%
  rename(Water_P_lag7daytol = Water_P.y)
result_fruit2 <- result_fruit2 %>%
  select(-ends_with(".y"))
result_fruit2 <- result_fruit2 %>% select(-Crop) 
result_fruit2 <- result_fruit2 %>% select(-Coord) 
result_fruit2 <- result_fruit2 %>% select(-Project) 
result_fruit2 <- result_fruit2 %>% select(-join_key)

save(result_fruit2, file = "df_final2_fruit_metbots.RData")

load("df_final2_fruit_metbots.RData")
