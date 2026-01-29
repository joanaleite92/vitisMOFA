#install packages if they do not exist yet
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")

#load packages
library(readxl)
library(dplyr)
library(tidyr)

#Set your parameters to read the file
file_path <- "your_data_file.xlsx" # Change to your file name
sheet_name_number <- 1                    
has_header <- TRUE                 

#load to dataframe
df <- read_excel(file_path, sheet = sheet_name, col_names = has_header)

# outlier detection functions
# IQR
is_outlier_iqr <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  return(x < lower | x > upper)
}

# z-score 
is_outlier_zscore <- function(x, threshold = 3) {
  z_scores <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  return(abs(z_scores) > threshold)
}

# find numeric columns and flag outliers
df_cleaned <- df %>%
  mutate(across(where(is.numeric), 
                list(outlier_iqr = ~is_outlier_iqr(.),
                     outlier_z = ~is_outlier_zscore(.)),
                .names = "{.col}_{.fn}"))


# See rows where the first numeric column has an IQR outlier
print(df_cleaned %>% filter(df_cleaned[[ncol(df) + 1]] == TRUE))

#see outliers only
outliers_list <- df_cleaned %>%
  pivot_longer(cols = where(is.numeric) & !contains("_outlier_"), 
               names_to = "Column", 
               values_to = "Value") %>%
  filter(get(paste0(Column, "_outlier_iqr")) == TRUE | 
           get(paste0(Column, "_outlier_z")) == TRUE) %>%
  select(Column, Value) %>%
  distinct()

print(outliers_list)
