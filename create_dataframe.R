## make a full  outer  join so we can join all datasets by id, but keeping the common columns unchanged and adding the different ones
library(dplyr)
library(purrr)

library(readxl)

df_ctr <- read_excel("Gene Expression_Omicbots_Cts_2022_2023_2024.xlsx", sheet = "L_GADPH")

# 2. Clean: Keep ONLY the 5 necessary columns and fix the decimal
df_ctr <- df_ctr %>%
  # Select only the columns that exist in your big dataframe
  dplyr::select(Date, Code_field, FW_sample, Irrigation, Laboratory, `GAPDH (Cq)`) %>%
  # Convert commas to dots and make numeric
  dplyr::mutate(
    # 1. Add the "AC_" prefix to Code_field
    Code_field = paste0("AC_", Code_field),
    
    # 2. Convert Cq: Replace comma with dot and make numeric
    `GAPDH (Cq)` = as.numeric(gsub(",", ".", `GAPDH (Cq)`)),
    
    # 3. Convert Date: 
    # (If Excel imported it as a number, use origin = "1899-12-30")
    Date = as.Date(Date),
    
    # 4. Irrigation to Factor
    Irrigation = as.factor(Irrigation),
    
    # 5. FW_sample to Ordered Factor
    FW_sample = factor(FW_sample, 
                       levels = c("S1", "S2", "S3", "S4", "S5", "S6"), 
                       ordered = TRUE)
  )

df_leaf_joined_complete_lagged <- df_leaf_joined_complete_lagged %>%
  rows_update(df_ctr, by = c("Date", "Code_field", "FW_sample", "Irrigation","Laboratory"))

df_leaf_joined_complete_lagged <- df_leaf_joined_complete_lagged %>%
  mutate(
    # We apply the logic across all three target columns
    across(
      c("RCA (Cq)", "LBCY (Cq)", "CHLG (Cq)"), 
      ~ ifelse(year(Date) == 2023 & FW_sample %in% c("S2", "S3"), NA, .)
    )
  )
save(df_leaf_joined_complete_lagged, file = "df_final8_leaf_metbots_noOut_starchCorrect_withlags.RData")
############################################
df_leaf_metabs2023 <- read_excel("Analises_lab_Douro_Aciprestes_2023.xlsx", sheet = "Sheet1",na = c("", "-"))
df_leaf_metabs2023 <- df_leaf_metabs2023 %>%
  mutate(
    Code_field = paste0("AC_", Code_field)
    
    # Alternativa usando stringr:
    # Code_field = str_c("AC_", Code_field) 
  )

df_leaf_spectra2023 <- read.csv(
  "LEAVES_FIELD_MB1_2023.csv", 
  text = lines[-1], 
  header = FALSE,
  #row.names = 1,
  encoding = "UTF-8"
)

new_headers <- as.character(round(as.numeric(df_leaf_spectra2023[1, ]),3))
#clean_names <- make.names(new_headers[-1], unique = TRUE)
colnames(df_leaf_spectra2023) <- new_headers
df_leaf_spectra2023 <- df_leaf_spectra2023[-c(1), -c(1)]
#######################################################

new_col_names <- c("Code_field", "FW_sample", "Date")
N_cols <- ncol(df_leaf_spectra2023)
colnames(df_leaf_spectra2023)[(N_cols - 2):N_cols] <- new_col_names

dim(df_leaf_spectra2023)
dim(df_leaf_metabs2023)

#verirication
id_col_name1 <- "Code_field"
id_col_name2 <- "FW_sample"

chaves_spectra <- df_leaf_spectra2023 %>% 
  select(!!sym(id_col_name1), !!sym(id_col_name2))

chaves_metabs <- df_leaf_metabs2023 %>% 
  select(!!sym(id_col_name1), !!sym(id_col_name2))

amostra_extra_df <- anti_join(chaves_metabs, chaves_spectra, by = c(id_col_name1, id_col_name2))

#Mostrar a linha completa no dataframe de metabólitos
linha_extra_completa <- df_leaf_metabs2023 %>%
  inner_join(amostra_extra_df, by = c(id_col_name1, id_col_name2))

linha_extra_spectra <- df_leaf_spectra2023 %>%
  inner_join(amostra_extra_df, by = c(id_col_name1, id_col_name2))
#######################################################

dataset_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

datasets <- mget(ls()[sapply(ls(), function(x) is.data.frame(get(x)))], envir = .GlobalEnv)
#Check if ID column exist in every dataframe and standardize the name or ignore that dataframe if it doesn't

has_id <- sapply(datasets, function(df) "ID" %in% colnames(df))
# Get the names of datasets that do NOT have "ID"
#names(has_id)[!has_id]

# Put them into a list
datasets_to_join <- dataset_names[has_id]

datasets_to_join <- lapply(datasets_to_join, get)

df_merged <- Reduce(function(x, y) full_join(x, y, by = "ID"), datasets_to_join)

result <- df_merged %>%
  group_by(Code_field) %>%
  summarise(n_unique = n_distinct(Date))

head(result)

# Keep only those that have both "Leaf" and "Metabots" in their names
target_dfs <- dataset_names[grepl("Leaf", dataset_names) & grepl("Metbots", dataset_names)]
# Retrieve them as a list
target_dfs_list <- mget(target_dfs)

feature_idx <- c(2:10, 12:300)

check_features <- function(target_dfs_list, feature_idx) {
  # Take the feature part of the first df
  ref <- target_dfs_list[[1]][, c(1, feature_idx)]
  
  for (i in 2:length(target_dfs_list)) {
    test <- target_dfs_list[[i]][, c(1, feature_idx)]
    
    # Merge by ID to align rows
    merged <- merge(ref, test, by = "ID", suffixes = c(".ref", ".test"))
    
    # Prepare feature column pairs
    feature_pairs <- paste0(names(ref)[feature_idx], c(".ref", ".test"))
    
    # Create logical matrix comparing each feature
    diffs <- sapply(seq_along(feature_idx), function(j) {
      merged[[paste0(names(ref)[feature_idx[j]], ".ref")]] != 
        merged[[paste0(names(ref)[feature_idx[j]], ".test")]]
    })
    
    # Ensure diffs is a matrix
    if (is.null(dim(diffs))) diffs <- matrix(diffs, nrow = nrow(merged))
    
    # Check rows where any feature differs
    rows_with_diff <- apply(diffs, 1, any)
    
    if (any(rows_with_diff)) {
      bad_ids <- merged$ID[rows_with_diff]
      cat("Inconsistencies found in dataframe", i, "for IDs:\n")
      print(bad_ids)
    } else {
      cat("Dataframe", i, ": all features consistent ✅\n")
    }
  }
}
# Merge them all together (by common columns)

df_merged_leaf_metbots <- Reduce(function(x, y) {
 full_join(x, y, by = intersect(names(x), names(y)))
}, target_dfs_list)

df_merged_leaf_metbots <- Reduce(function(x, y) {
  full_join(x, y, by = "ID")
}, target_dfs_list)

#merge only features first
target_cols <- c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_","MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                 "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)","LBCY (Cq)","CHLG (Cq)","GAPDH (Cq)","GA3 (µg/mgFM)","BAP (µg/mgFM)","IAA (µg/mgFM)",
                 "SA (µg/mgFM)","ABA (µg/mgFM)","Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_","Sugar__mg_gFM_","Anthocyanins__mg_gFM_")
# Feature columns = all except ID and target columns
feature_cols_list <- lapply(target_dfs_list, function(df) {
  # Only exclude target columns that exist in this dataframe
  cols_to_exclude <- intersect(names(df), c("ID", target_cols))
  setdiff(names(df), cols_to_exclude)
})

common_features <- Reduce(intersect, feature_cols_list)

feature_dfs <- lapply(target_dfs_list, function(df) {
  df[, c("ID", common_features)]
})

df_features_merged <- reduce(feature_dfs, full_join, by = c("ID", common_features))
save(df_features_merged, file = "df_features_leaf_metbots_merged.RData")

# If you want overlapping columns with same names merged instead of duplicated:
df_merged_leaf_metbots <- df_merged_leaf_metbots %>%
  mutate(across(ends_with(".x"), ~ {
    base <- sub("\\.x$", "", cur_column())
    ycol <- paste0(base, ".y")
    if (ycol %in% names(.)) {
      coalesce(.x, .[[ycol]])
    } else {
      .x
    }
  })) %>%
  select(-ends_with(".y")) %>%
  rename_with(~ sub("\\.x$", "", .x), ends_with(".x"))

save(df_merged_leaf_metbots, file = "df_merged_leaf_metbots.RData")

colnames(df_merged_leaf_metbots)
## Count NAs per column
colSums(is.na(df_merged_leaf_metbots))

# Check if there are any duplicates in the ID column
any(duplicated(df_merged_leaf_metbots$ID))
# TRUE if duplicates exist, FALSE if not

# Get the duplicated rows
df[duplicated(df_merged_leaf_metbots$ID) | duplicated(df_merged_leaf_metbots$ID, fromLast = TRUE), ]

# Get only the duplicated IDs (unique)
unique(df_merged_leaf_metbots$ID[duplicated(df_merged_leaf_metbots$ID)])

dup_ids <- df_merged_leaf_metbots$ID[duplicated(df_merged_leaf_metbots$ID) | duplicated(df_merged_leaf_metbots$ID, fromLast = TRUE)]

df_duplicates <- df_merged_leaf_metbots[df_merged_leaf_metbots$ID %in% dup_ids, ]

id_to_check <- df_duplicates$ID[1]
id_to_check <- "5312"#"5406"

rows <- df_duplicates[df_duplicates$ID == id_to_check, ]

# Identify columns where values differ
cols_with_diff <- sapply(rows, function(col) length(unique(col)) > 1)
# Print columns that differ
names(rows)[cols_with_diff]

#Check which dataframes are contributing the duplicated IDs
#Strategy is adding a source column before merging

# Add a source label to each df in the list
# Add a Source column using the dataframe names in the list
feature_dfs_tagged <- Map(function(df, nm) {
  df %>% mutate(Source = nm)
}, target_dfs_list, names(target_dfs_list))

# Keep only ID + features + Source
feature_dfs_tagged <- lapply(feature_dfs_tagged, function(df) {
  df[, c("ID", common_features, "Source")]
})


# Merge everything
features_merged_tagged <- bind_rows(feature_dfs_tagged)

## check the duplicates
dup_ids <- features_merged_tagged %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  arrange(ID)

dup_ids

#drop duplicates only when both ID and Laboratory are the same (i.e., same sample measured in the same lab).
df_dedup <- df_merged_leaf_metbots[!duplicated(df_merged_leaf_metbots[c("ID", "Laboratory")]), ]

## Keep duplicates across labs (same ID, different Laboratory), but remove only the “true duplicates” (same ID + same lab + same features)… then we’d expand the distinct() to include more columns.
#f the same ID appears in different labs → ✅ keep both.
#If the same ID + same lab but features differ → ✅ keep both (they may represent different measurements).
#If the same ID + same lab + features are identical → ❌ drop duplicates.
feature_cols <- c("Feature1", "Feature2", "Feature3")  # adjust
df_dedup <- df[!duplicated(df[c("ID", "Laboratory", feature_cols)]), ]

# as above, but the first condition changes from keep to drop, the other two conditions remain the same
df_dedup <- df %>%
  # First, drop exact duplicates of ID + Lab + features
  distinct(ID, Laboratory, !!!syms(feature_cols), .keep_all = TRUE) %>%
  # Then, drop all duplicates of ID across labs, keeping the first occurrence
  distinct(ID, .keep_all = TRUE)

save(df_dedup, file = "df_final.RData")


################################################################################################
###FRUIT
###########################################################

# Keep only those that have both "Fruit" and "Metabots" in their names
target_dfs_fruit <- dataset_names[grepl("Fruit", dataset_names) & grepl("Metbots", dataset_names)]
# Retrieve them as a list
target_dfs_list_fruit <- mget(target_dfs_fruit)

feature_idx_fruit <- c(2:10, 12:300)

# Merge them all together (by common columns)

#RUN THIS ONE INSTEAD OF THE ONE BELOW
df_merged_fruit_metbots <- Reduce(function(x, y) {
  full_join(x, y, by = intersect(names(x), names(y)))
}, target_dfs_list_fruit)

df_merged_fruit_metbots <- Reduce(function(x, y) {
  full_join(x, y, by = "ID")
}, target_dfs_list_fruit)

target_cols_fruit <- c("Alpha_Amino_Nitrogen__mg_L","Ammoniacal_Nitrogen__mg_L",
                      "Assimilable_Nitrogen__mg_L","Brix_lab__%_","Brix_tom__%_",
                      "CHS (Cq)", "PAL (Cq)","PAL1 (Cq)","UFGT1 (Cq)","UFGT (Cq)", "GAPDH (Cq)")

feature_cols_list_fruit <- lapply(target_dfs_list_fruit, function(df) {
  # Only exclude target columns that exist in this dataframe
  cols_to_exclude <- intersect(names(df), c("ID", target_cols_fruit))
  setdiff(names(df), cols_to_exclude)
})

common_features_fruit <- Reduce(intersect, feature_cols_list_fruit)

feature_dfs_fruit <- lapply(target_dfs_list_fruit, function(df) {
  df[, c("ID", common_features_fruit)]
})

df_features_merged_fruit <- reduce(feature_dfs_fruit, full_join, by = c("ID", common_features_fruit))
save(df_merged_fruit_metbots, file = "df_fruit_metbots_merged.RData") 

df_dedup_fruit <- df_merged_fruit_metbots[!duplicated(df_merged_fruit_metbots[c("ID", "Laboratory")]), ]
#summary(df_dedup_fruit)
na_counts <- colSums(is.na(df_dedup_fruit))
na_counts

