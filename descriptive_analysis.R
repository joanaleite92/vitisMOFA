library(ggplot2)
library(dplyr)
library(corrplot)
library(reshape2)
library(tidyr)
library(ggrepel)
library(corpcor)
library(RColorBrewer)

#1. Summarize numeric and categorical variables
#2. Visualize temporal trends (dates, timepoints, years)
#3. Explore correlations between predictors and targets
#4. Check distributions and outliers
#5. Perform multivariate exploration (PCA, clustering)
#6. Inspect relationships between climate averages and metabolites
#7. Check experimental design balance (timepoint × irrigation × year)


#1. Basic summary statistics

##Numeric variables (spectral bands, climate averages, metabolites, gene expression)
summary(result)
sapply(result[, sapply(result, is.numeric)], function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE), min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE)))

##Categorical variables (timepoints, irrigation, year)
table(result$FW_sample)
table(result$Irrigation)
table(result$year)

### Metabolites vs Gene Expression###
## correlation matrices



metabs_cols <- c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                 "Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                 "Sugar__mg_gFM_","Anthocyanins__mg_gFM_")
genes_cols <- c("MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
                "LBCY (Cq)","CHLG (Cq)","GAPDH (Cq)")
metabs <- result2 %>% select(all_of(metabs_cols))
genes <- result2 %>% select(all_of(genes_cols))
genes_expr <- 2^(-genes)

cor_mat <- cor(metabs, genes, method = "spearman", use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", 
         tl.col = "black", tl.cex = 0.7,    # text labels
         number.cex = 0.7, addCoef.col = "black")  # add correlation values

ggsave("cor_matrix_spearman_metabsVsgenes.pdf")


###### FRUIT ##########
ggplot(result_fruit2_7daybef, aes(x = Water_P.y , y = `Brix_tom__%_`)) +
  geom_point(alpha = 0.3) +                          # scatter points
  geom_smooth(method = "gam", formula = y ~ s(x),    # non-linear smoother
              se = TRUE, color = "red") +
  labs(
    x = "PHB (Water_P)",
    y = "Brix or Score",
    title = "Scatter plot with non-linear trend"
  ) +
  theme_minimal()

metabs_cols_fruit <- c("Alpha_Amino_Nitrogen__mg_L_","Ammoniacal_Nitrogen__mg_L_",
"Assimilable_Nitrogen__mg_L_","Brix_lab__%_","Brix_tom__%_")

genes_cols_fruit <- c("UFGT (Cq)","PAL (Cq)","CHS (Cq)","UFGT1 (Cq)","PAL1 (Cq)")

metabs_fruit <- result_fruit %>% select(all_of(metabs_cols_fruit))
genes_fruit <-result_fruit  %>% select(all_of(genes_cols_fruit))
genes_expr_fruit <- 2^(-genes_fruit)
# Keep only rows where both metabolites and genes are non-missing
rows_ok <- complete.cases(result_fruit[metabs_cols_fruit], result_fruit[genes_cols_fruit])
sum(rows_ok)  # how many complete rows?

metabs_fruit_clean <- result_fruit[rows_ok, metabs_cols_fruit]
genes_fruit_clean  <- result_fruit[rows_ok, genes_cols_fruit]
###fruit
cor_mat_metabs <- cor(metabs_fruit, method = "spearman", use = "pairwise.complete.obs")
cor_mat_genes <- cor(genes_fruit, method = "spearman", use = "pairwise.complete.obs")
cor_mat_genes_expr <- cor(genes_expr_fruit, method = "spearman", use = "pairwise.complete.obs")
corrplot(cor_mat_genes_expr,
         method = "color",
         type = "upper",          # upper triangle only
         tl.col = "black",
         tl.cex = 0.8,
         number.cex = 0.7,
         addCoef.col = "black")

library(pheatmap)
plot_cor_matrix <- function(cor_mat) {
  # 1. Replace NA with 0 (safer for viz)
  cor_mat_clean <- cor_mat
  cor_mat_clean[is.na(cor_mat_clean)] <- 0
  
  # 2. Drop empty rows/cols (all 0 after NA replace)
  cor_mat_clean <- cor_mat_clean[rowSums(cor_mat_clean != 0) > 0, , drop = FALSE]
  cor_mat_clean <- cor_mat_clean[, colSums(cor_mat_clean != 0) > 0, drop = FALSE]
  
  # 3. Check dimensions
  if (nrow(cor_mat_clean) < 2 | ncol(cor_mat_clean) < 2) {
    message("Matrix too small for clustering — showing without dendrograms")
    pheatmap(cor_mat_clean,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             display_numbers = TRUE)
  } else {
    pheatmap(cor_mat_clean,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             display_numbers = TRUE)
  }
}

cor_mat <- cor(metabolites, genes,
               method = "spearman",
               use = "pairwise.complete.obs")

plot_cor_matrix(cor_mat)

# Combine into one numeric matrix
dat_all <- cbind(metabs, genes)
dat_all <- dat_all %>% dplyr::mutate_all(as.numeric)  # ensure numeric
#dat_all <- na.omit(dat_all)  # drop rows with missing values

#Compute partial correlations (shrinkage)
# Spearman correlation first
cor_mat <- cor(dat_all, method = "spearman",use = "pairwise.complete.obs")

# Convert to partial correlation matrix (stable)
pcor_mat <- cor2pcor(cor_mat)


#Pathway-Level Correlation Analysis

#Flav_pathway ↔ Anthocyanins to be stronger than any single gene.
#Phenyl_pathway ↔ Phenols should also improve.
#Chlorophyll_pathway ↔ Chlorophyll a/b might come out clearer.
#PROS negatively correlating with flavonoids)

# Flavonoid / Anthocyanin pathway
flav_genes <- c("MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)",
                "CHS1 (Cq)","DFR (Cq)","UFGT1 (Cq)",
                "MATE1 (Cq)","ABCC1 (Cq)")

# Phenylpropanoid core (general phenolics, PAL, C4H)
phenyl_genes <- c("PAL1 (Cq)","C4H1 (Cq)")

# Flavonols
fls_genes <- c("FLS1 (Cq)")

# Chlorophyll biosynthesis
chloro_genes <- c("RCA (Cq)","LBCY (Cq)","CHLG (Cq)")

# Housekeeping (not expected to correlate)
ref_gene <- c("GAPDH (Cq)")

pathway_scores <- result %>%
  mutate(
    Flav_pathway = rowMeans(select(., all_of(flav_genes)), na.rm = TRUE),
    Phenyl_pathway = rowMeans(select(., all_of(phenyl_genes)), na.rm = TRUE),
    Flavonol_pathway = rowMeans(select(., all_of(fls_genes)), na.rm = TRUE),
    Chlorophyll_pathway = rowMeans(select(., all_of(chloro_genes)), na.rm = TRUE)
  )


#pathway-level expression scores: mean (or PCA1)
pathway_scores <- result %>%
  mutate(
    Flav_pathway = rowMeans(select(., all_of(flav_genes)), na.rm = TRUE),
    Phenyl_pathway = rowMeans(select(., all_of(phenyl_genes)), na.rm = TRUE),
    Flavonol_pathway = rowMeans(select(., all_of(fls_genes)), na.rm = TRUE),
    Chlorophyll_pathway = rowMeans(select(., all_of(chloro_genes)), na.rm = TRUE)
  )



#Correlate pathway scores with metabolites
cor_mat_pathways <- cor(
  pathway_scores %>% select(Flav_pathway, Phenyl_pathway, Flavonol_pathway, Chlorophyll_pathway),
  pathway_scores %>% select(all_of(metabs_cols)),
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_long <- melt(cor_mat_pathways)

ggplot(cor_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       limit = c(-1, 1), name = "Spearman r") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Assuming:
# - pathway_scores contains the Cq values and pathway scores
# - chloro_genes = c("RCA (Cq)", "LBCY (Cq)", "CHLG (Cq)")
# - chlorophyll metabolites are named in metab_cols, e.g., "Chlorophyll_a", "Chlorophyll_b"

# Step 1: Make a logical matrix for non-missing gene values
gene_nonNA <- pathway_scores %>%
  select(all_of(chloro_genes)) %>%
  mutate(across(everything(), ~ !is.na(.)))

metab_nonNA <- pathway_scores %>%
  select(Chlorophyll_a__mg_gFM_, Chlorophyll_b__mg_gFM_) %>%
  mutate(across(everything(), ~ !is.na(.)))

# Overlap with Chlorophyll_a
sample_overlap_a <- colSums(gene_nonNA & metab_nonNA$Chlorophyll_a__mg_gFM_)

# Overlap with Chlorophyll_b
sample_overlap_b <- colSums(gene_nonNA & metab_nonNA$Chlorophyll_b__mg_gFM_)

sample_overlap_a
sample_overlap_b





#2. Temporal trends
##Plot metabolites or gene expression across timepoints and years
# 1. Identify target columns (numeric metabolites/genes)
target_cols <- c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_","MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                 "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)","LBCY (Cq)","CHLG (Cq)","GAPDH (Cq)","GA3 (µg/mgFM)","BAP (µg/mgFM)","IAA (µg/mgFM)",
                 "SA (µg/mgFM)","ABA (µg/mgFM)","Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_","Sugar__mg_gFM_","Anthocyanins__mg_gFM_")
# 2. Reshape to long format
result_long <- result2_noOut %>%
  pivot_longer(
    cols = all_of(target_cols),
    names_to = "Target",
    values_to = "Value"
  )

# 3. Plot all targets as a matrix of boxplots
ggplot(result_long %>% dplyr::filter(!is.na(Value)),
       aes(x = FW_sample, y = Value, color = Test_site)) +
  geom_boxplot() +
  theme_bw() +
  labs(color = "Test_site",
       x = "Timepoint",
       y = "Concentration / Expression") +
  theme(legend.position = "right") +
  facet_wrap(~ Target, scales = "free_y", ncol = 5)  # adjust ncol for layout
ggsave("all_targets_boxplots_perLocation_noOut.pdf", width = 16, height = 12)

ggplot(result_long %>% dplyr::filter(!is.na(Value)),
       aes(x = FW_sample, y = Value, color = year)) +
  geom_boxplot() +
  theme_bw() +
  labs(color = "Year",
       x = "Timepoint",
       y = "Concentration / Expression") +
  theme(legend.position = "right") +
  facet_wrap(~ Target, scales = "free_y", ncol = 5)  # adjust ncol for layout

ggsave("all_targets_boxplots_perYear_noOut.pdf", width = 16, height = 12)

ggsave("all_targets_boxplots.png", width = 16, height = 12,dpi=300)

ggplot(result_long %>% dplyr::filter(!is.na(Value)),
       aes(x = FW_sample, y = Value, color = interaction(Test_site, year))) +
  geom_boxplot() +
  theme_bw() +
  labs(color = "Test_site-Year",
       x = "Timepoint",
       y = "Concentration / Expression") +
  theme(legend.position = "right") +
  facet_wrap(~ Target, scales = "free_y", ncol = 5)

ggsave("all_targets_boxplots_interact_year_site.pdf", width = 16, height = 12)

ggplot(result_long %>% dplyr::filter(!is.na(Value)),
       aes(x = FW_sample, y = Value, color = Test_site, fill = year)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  theme_bw() +
  labs(color = "Test_site", fill = "Year",
       x = "Timepoint",
       y = "Concentration / Expression") +
  theme(legend.position = "right") +
  facet_wrap(~ Target, scales = "free_y", ncol = 5)

ggsave("all_targets_boxplots_yearFill_siteShapecolor.pdf", width = 16, height = 12)

ggplot(result_long %>% dplyr::filter(!is.na(Value)),
       aes(x = year, y = Value, color = Test_site)) +
  geom_boxplot() +
  theme_bw() +
  labs(color = "Test_site",
       x = "Year",
       y = "Concentration / Expression") +
  theme(legend.position = "right") +
  facet_wrap(~ Target, scales = "free_y", ncol = 5)

ggsave("all_targets_boxplots_year_perLocation.pdf", width = 16, height = 12)

#3. Correlation analysis
##Between spectral features and targets, or climate variables and metabolites
numeric_vars <- result[, sapply(result, is.numeric)]
# Keep only numeric columns with >1 non-NA value and non-zero SD
numeric_vars_clean <- numeric_vars %>%
  select(where(function(x) {
    vals <- x[!is.na(x)]
    length(vals) > 2 && sd(vals) != 0
  }))

spectral_vars <- result %>% select(starts_with("wl_"))     # spectral columns
aux_vars  <- result %>% select(329:342) # climate averages
target_vars   <- target_cols

#ggpairs(numeric_vars)

#Compute a correlation matrix: helps identify which spectral regions or climate variables are associated with metabolites or gene expression
#cor(numeric_vars, use="complete.obs")
# "complete.obs": only rows with no NA in any column
# "pairwise.complete.obs": use all available pairs


cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")
cor_matrix[lower.tri(cor_matrix)] <- NA
cor_melt <- melt(cor_matrix, na.rm = FALSE)  # keep NAs
# Keep only upper triangle positions
cor_melt_upper <- cor_melt %>% 
  dplyr::filter(as.numeric(Var1) < as.numeric(Var2))

#spearman
cor_matrix_spearman <- cor(numeric_vars, use = "pairwise.complete.obs", method = "spearman")
cor_matrix_spearman[lower.tri(cor_matrix)] <- NA
cor_melt_spearman <- melt(cor_matrix_spearman, na.rm = FALSE)  # keep NAs
# Keep only upper triangle positions
cor_melt_upper_spearman <- cor_melt_spearman %>% 
  dplyr::filter(as.numeric(Var1) < as.numeric(Var2))
ggplot(cor_melt_spearman, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,na.value="grey90") +
  scale_x_discrete(position = "top") +  # horizontal labels at the top
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(fill = "Correlation")

ggsave("cor_matrix_spearman.pdf", width = 60, height = 40,limitsize=FALSE)
#corrplot(cor_matrix, method = "color", type = "upper", 
#        order = "hclust", 
 #       tl.col = "black", tl.srt = 45)

#4. Distribution plots
## Histograms/density plots for metabolites and gene expression to check normality
ggplot(result_long, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_wrap(~ Target, scales = "free") +  # free scales for different ranges
  theme_bw() +
  labs(x = "Value", y = "Count", title = "Distribution of Metabolites and Gene Expression")
ggplot(result, aes(x = metabolite1)) + geom_histogram(bins = 30, fill="skyblue", color="black")

ggsave("all_targets_histograms.pdf", width = 16, height = 12)

#5. Multivariate analysis
#PCA on spectral data to explore patterns
spectra_matrix <- as.matrix(result[, grep("^wl_", names(result))]) # assuming spectral columns S1, S2...
pca <- prcomp(spectra_matrix, scale.=TRUE)
biplot(pca)
#Clustering to see if samples group by timepoint, irrigation, or year

# Function for clean PCA biplot
plot_pca_biplot <- function(df, spectral_prefix = "wl_", top_n = 20, sample_group = NULL, 
                            scale_factor = 5, color_palette = NULL,na_color = "grey90",
                            na_shape = 17, default_shape = 16) {
  #Automatically selects spectral columns starting with wl_
  #Performs PCA on scaled data
  #Shows only the top N contributing wavelengths
  #Colors samples by an optional grouping variable (year, timepoint, treatment…)
  #Uses geom_text_repel to avoid overlapping wavelength labels
  #Axis labels include % variance explained
  
  # 1. Select spectral columns
  spectral_cols <- grep(paste0("^", spectral_prefix), names(df), value = TRUE)
  spectra_matrix <- as.matrix(df[, spectral_cols])
  
  # 2. PCA
  pca <- prcomp(spectra_matrix, scale. = TRUE)
  
  # 3. Scores (samples)
  scores <- as.data.frame(pca$x[,1:2])
  scores$Sample <- rownames(scores)
  
  # 4. Grouping variable
  is_group_discrete <- FALSE
  if (!is.null(sample_group)) {
    scores$Group <- df[[sample_group]]
    
    # Detect if categorical
    if (is.factor(scores$Group) || is.character(scores$Group)) {
      is_group_discrete <- TRUE
      scores$Group <- as.factor(scores$Group)
    } else {
      is_group_discrete <- FALSE
    }
    
    # Shape factor: NA vs not NA
    scores$ShapeGroup <- factor(ifelse(is.na(scores$Group), "NA", "Not NA"),
                                levels = c("Not NA", "NA"))
  } else {
    scores$ShapeGroup <- factor("Not NA")
  }
  
  # 5. Loadings (variables)
  loadings <- as.data.frame(pca$rotation[,1:2])
  loadings$Wavelength <- rownames(loadings)
  loadings$Contribution <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  
  # 6. Top N loadings
  top_loadings <- loadings %>% arrange(desc(Contribution)) %>% head(top_n)
  
  # 7. Base plot → points first
  p <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Group, shape = ShapeGroup), size = 3)
  
  # 8. Setas e labels por cima
  p <- p +
    geom_segment(data = top_loadings,
                 aes(x = 0, y = 0, xend = PC1*scale_factor, yend = PC2*scale_factor),
                 arrow = arrow(length = unit(0.2,"cm")), color = "darkblue", linewidth = 1.5,
                 inherit.aes = FALSE) +
    geom_text_repel(data = top_loadings,
                    aes(x = PC1*scale_factor, y = PC2*scale_factor, label = Wavelength),
                    color = "darkblue", size = 4, max.overlaps = Inf, force = 2,
                    inherit.aes = FALSE)
  
  # 9. Axis labels and theme
  p <- p +
    theme_bw() +
    labs(
      x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)")
    )
  
  # 10. Color and shape scales
  if (!is.null(sample_group)) {
    if (is_group_discrete) {
      n_levels <- length(unique(scores$Group))
      if (is.null(color_palette)) color_palette <- brewer.pal(min(n_levels, 6), "Set1")
      p <- p + scale_color_manual(values = color_palette, na.value = na_color) +
        scale_shape_manual(values = c("Not NA" = default_shape, "NA" = na_shape)) +
        labs(color = sample_group, shape = "Missing")
    } else {
      if (is.null(color_palette)) color_palette <- colorRampPalette(c("skyblue", "#ff0000"))(6)
      p <- p + scale_color_gradientn(colors = color_palette, na.value = na_color) +
        scale_shape_manual(values = c("Not NA" = default_shape, "NA" = na_shape)) +
        labs(color = sample_group, shape = "Missing")
    }
  } else {
    p <- p + scale_shape_manual(values = c("Not NA" = default_shape, "NA" = na_shape),
                              breaks = "NA", labels = "NA", name = "Missing")
  }
  
  return(p)
}
# Example: spectral columns start with "wl_", top 20 wavelengths labeled, color points by timepoint
plot_pca_biplot(df = result, spectral_prefix = "wl_", top_n = 25, sample_group = "year",na_shape = 4, default_shape = 16)


#### -------####
ggplot(result, aes(x = mean_Tmax, y = Carotenoids__mg_gFM_, color = factor(year))) +
  geom_point(aes(shape = factor(year)), size = 3) +
  geom_line(aes(group = year), linewidth = 1) +   # draws a line connecting points per year
  theme_bw() +
  labs(x = "Max Temperature", y = "Carotenoids (mg/g FM)", color = "Year", shape = "Year")

#6. Time-series / rolling averages
##Plot 7-day climate averages vs metabolite concentrations
ggplot(df, aes(x = mean_Tmed, y = metabolite1, color = timepoint_ord)) + geom_point() + geom_smooth(method="lm")

##Check if recent climate conditions correlate with metabolic profiles

# 7. Cross-tabulations and contingency analysis
table(df$timepoint_ord, df$irrigation)
#uneven sampling or interaction patterns