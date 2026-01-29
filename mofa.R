# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
# Install DESeq2
BiocManager::install("DESeq2")
install.packages("KEGGREST", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("future", dependencies = TRUE)
install.packages("hardhat", dependencies = TRUE)
install.packages("MatrixModels", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("future.apply", dependencies = TRUE)
install.packages("recipes", dependencies = TRUE)
install.packages("lava", dependencies = TRUE)
BiocManager::install(c("biomaRt", "KEGGREST"))

if (!requireNamespace("org.Vv.eg.db", quietly = TRUE))
  BiocManager::install("org.Vv.eg.db")

library(org.Vv.eg.db)

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

library(clusterProfiler)

library(MOFA2)
library(DESeq2)
library(e1071)
library(pheatmap)
library(dplyr)
library(caret)
library(ggplot2)
library(reshape2)
library(tidyr)
library(limma)
library(markdown)
library(ggtext)
library(xfun)


# ---- 1. Gene expression preprocessing ----
preprocess_genes <- function(genes_mat, type = c("counts","cq")) {
  type <- match.arg(type)
  
  if (type == "counts") {
    # RNA-seq raw counts → VST
    dds <- DESeqDataSetFromMatrix(countData = t(genes_mat), 
                                  colData = data.frame(row.names=colnames(genes_mat)), 
                                  design = ~1)
    vsd <- vst(dds)
    mat <- t(assay(vsd))
    
  } else if (type == "cq") {
    # qPCR Cq values
    # Step 1: global mean ΔCq normalization
    sample_means <- rowMeans(genes_mat, na.rm = TRUE)
    delta_cq <- sweep(genes_mat, 1, sample_means, "-")
    
    # Step 2: flip sign so higher = more expression
    expr <- -delta_cq
    
    # Step 3: center & scale - careful with this step, probably not needed???!!!!
    #mat <- scale(expr, center = TRUE, scale = TRUE)
  }
  
  return(as.matrix(expr))
}

# ---- 2. Metabolites preprocessing ----
preprocess_metabolites <- function(metab_mat, skew_threshold = 1) {
  M <- metab_mat
  skews <- apply(M, 2, function(x) e1071::skewness(x, na.rm=TRUE))
  do_log <- skews > skew_threshold & apply(M, 2, function(x) all(x >= 0, na.rm=TRUE))
  
  for (j in which(do_log)) {
    x <- M[, j]
    offset <- min(x[x > 0], na.rm=TRUE)/2
    if (offset <= 0) offset <- 1
    M[, j] <- log(x + offset)
  }
  
  M_z <- scale(M, center = TRUE, scale = TRUE)
  return(as.matrix(M_z))
}

# ---- 3. Spectra preprocessing ----
preprocess_spectra <- function(spec_mat) {
  spec_z <- scale(spec_mat, center = TRUE, scale = TRUE)
  return(as.matrix(spec_z))
}
load("df_final5_leaf_metbots_noOut_starchCorrect.RData")
load("df_final7_leaf_metbots_noOut_starchCorrect_withlags.RData")
load("df_final2_fruit_metbots.RData")
load("df_final6_leaf_metbots_noOut_starchCorrect_withlags.RData")
result2_noOut_starch2 <- result2_noOut_starch2 %>%
  filter(Test_site != "Greenhouse Porto")
result2_noOut_starch2 <- result2_noOut_starch2 %>%
  filter(!Test_site %in% c("Greenhouse Porto", "Quinta de Vale de Cavalos"))

df_leaf_joined_complete <- df_leaf_joined_complete%>%
  filter(!Test_site %in% c("Greenhouse Porto", "Quinta de Vale de Cavalos"))


df <- result_fruit2
resultado <- df %>%
  dplyr::filter(
    Test_site == "Quinta de Vale de Cavalos",
    year(Date) == 2024
  ) %>%
  # Seleciona a coluna específica E as que têm o prefixo
  dplyr::select(Code_field, starts_with("wl_"))

library(stringr)

# 1. Filtrar apenas as amostras que terminam em ".a" para o local pretendido
amostras_a <- result_fruit2 %>%
  filter(str_detect(Test_site, "Quinta dos Aciprestes")) %>%
  filter(str_ends(Code_field, "\\.a"))

# 2. Contar quantos Code_field distintos temos
n_distintos <- n_distinct(amostras_a$Code_field)

print(paste("Encontrados", n_distintos, "Code_fields terminados em 'a'"))

# 3. Ver quais são e o estado dos dados neles
resumo_amostras_a <- amostras_a %>%
  dplyr::select(Code_field,Test_site,FW_sample, Date,
         `GAPDH (Cq)`, `UFGT (Cq)`, `PAL (Cq)`, `CHS (Cq)`,
         `Brix_tom__%_`, `Brix_lab__%_`, `Alpha_Amino_Nitrogen__mg_L_`, `Ammoniacal_Nitrogen__mg_L_`, 
         `Assimilable_Nitrogen__mg_L_`)

print(resumo_amostras_a)

resumo_amostras_a <- amostras_a %>%
  dplyr::select(Code_field,Test_site,FW_sample, Date,
                `GAPDH (Cq)`, `UFGT (Cq)`, `PAL (Cq)`, `CHS (Cq)`,
                `Brix_tom__%_`, `Brix_lab__%_`, `Alpha_Amino_Nitrogen__mg_L_`, `Ammoniacal_Nitrogen__mg_L_`, 
                `Assimilable_Nitrogen__mg_L_`)

print(resumo_amostras_a)


tabela_final_final <- result_fruit2 %>%
  # 1. Filtrar pelo local
  filter(str_detect(Test_site, "Aciprestes")) %>%
  
  # 2. Excluir amostras "d"
  filter(!str_ends(Code_field, "\\.d")) %>%
  
  # 3. Excluir o ano de 2021
  # (Garante que o R percebe que a coluna Date é uma data)
  mutate(Date = as.Date(Date)) %>% 
  filter(year(Date) != 2021) %>%
  
  # 4. Criar o código limpo (remover a, b, c)
  mutate(Code_field_limpo = str_remove(Code_field, "\\.[a-z]$")) %>%
  
  # 5. Agrupar e Consolidar
  group_by(Code_field_limpo, Test_site, FW_sample, Date) %>%
  summarise(
    across(
      c(`GAPDH (Cq)`, `UFGT (Cq)`, `PAL (Cq)`, `CHS (Cq)`,
        `Brix_tom__%_`, `Brix_lab__%_`, `Alpha_Amino_Nitrogen__mg_L_`, 
        `Ammoniacal_Nitrogen__mg_L_`, `Assimilable_Nitrogen__mg_L_`),
      ~ mean(., na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  
  # Limpeza final de NaNs
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

# Ver resultado
print(tabela_final_final, n = Inf)

# Ver quantas sobraram agora
cat("\nTotal de campos após filtrar 2021 e amostras 'd':", nrow(tabela_final_final))


# 1. Definir os grupos de colunas
cols_metab <- c("Alpha_Amino_Nitrogen__mg_L_", "Ammoniacal_Nitrogen__mg_L_", 
                "Assimilable_Nitrogen__mg_L_", "Brix_lab__%_", "Brix_tom__%_")

cols_genes <- c("GAPDH (Cq)", "UFGT (Cq)", "PAL (Cq)", "CHS (Cq)")

# 2. Processar e gerar o diagnóstico de NAs
diagnostico_campos <- result_fruit2 %>%
  # Filtros anteriores (Aciprestes, sem 'd', sem 2021)
  filter(str_detect(Test_site, "Aciprestes")) %>%
  filter(!str_ends(Code_field, "\\.d")) %>%
  mutate(Date = as.Date(Date)) %>% 
  filter(year(Date) != 2021) %>%
  mutate(Code_field_limpo = str_remove(Code_field, "\\.[a-z]$")) %>%
  
  # Agrupar e consolidar médias (a, b, c)
  group_by(Code_field_limpo, Test_site, FW_sample, Date) %>%
  summarise(
    across(all_of(c(cols_metab, cols_genes)), ~ mean(., na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  
  # 3. Calcular as contagens de preenchimento
  mutate(
    # Quantos dos 5 metabolitos NÃO são NA?
    metabolitos_preenchidos = rowSums(!is.na(dplyr::select(., all_of(cols_metab)))),
    
    # Quantos dos 4 genes SÃO NA?
    genes_em_falta = rowSums(is.na(dplyr::select(., all_of(cols_genes))))
  )

# 4. Ver o resumo ordenado pelos que têm mais metabolitos e menos genes
resultado_final <- diagnostico_campos %>%
  dplyr::select(Code_field_limpo, FW_sample, Date, metabolitos_preenchidos, genes_em_falta) %>%
  arrange(desc(metabolitos_preenchidos), genes_em_falta)

print(resultado_final, n = Inf)

# 5. Resumo estatístico para a tua decisão
cat("\n--- RESUMO PARA DECISÃO ---")
cat("\nTotal de campos analisados:", nrow(resultado_final))
cat("\nCampos com TODOS os metabolitos (5/5):", sum(resultado_final$metabolitos_preenchidos == 5))
cat("\nCampos com TODOS os genes em falta (4/4):", sum(resultado_final$genes_em_falta == 4))
cat("\nCampos 'Perfeitos' (5 metabs e 0 genes em falta):", 
    sum(resultado_final$metabolitos_preenchidos == 5 & resultado_final$genes_em_falta == 0))


################################################
df_lagged_leaf <- df_lagged_leaf  %>%
  filter(Test_site != "Quinta de Vale de Cavalos",
  Cultivar != "Moscatel Galego Branco")
#####################################
##### Prepare blocks ######
result = df_leaf_joined_complete_lagged

result = result2_noOut_starch2
result = df_lagged_leaf
result = df_leaf_joined_complete_lagged
result_fruit = result_fruit2

colunas_a_contar <- c("MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                      "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
                      "LBCY (Cq)","CHLG (Cq)","Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                      "Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                      "Sugar__mg_gFM_","Anthocyanins__mg_gFM_") 

# 2. Executar o agrupamento e a contagem
contagem_observacoes <- result2_noOut_starch2 %>%
  # Agrupar pelo local de teste
  group_by(Test_site) %>%
  
  # Criar novas colunas (uma para cada coluna especificada em 'colunas_a_contar')
  summarise(
    # Usar across() para aplicar a mesma função a várias colunas
    across(
      .cols = all_of(colunas_a_contar), # Aplica-se apenas às colunas definidas
      .fns = ~ sum(!is.na(.)),          # A função: Conta o número de valores que NÃO são NA
      .names = "Contagem_{.col}"        # Nomeia as novas colunas como 'Contagem_Coluna_A', etc.
    ),
    .groups = "drop" # Remove o agrupamento
  )

tabela_longa <- contagem_observacoes %>%
  pivot_longer(
    cols = starts_with("Contagem_"),
    names_to = "Variavel",
    values_to = "Contagem"
  ) %>%
  # Limpar os nomes das variáveis (remover o prefixo "Contagem_")
  mutate(
    Variavel = gsub("Contagem_", "", Variavel)
  )

# 2. Pivotar para o formato wide (usa Test_site como cabeçalho de coluna)
tabela_transposta <- tabela_longa %>%
  pivot_wider(
    names_from = Test_site,
    values_from = Contagem
  )

# Visualizar a tabela transposta
print(tabela_transposta)

load("Joana_preto_branco.RData")
spectra_white_black <- as.data.frame(Joana_preto_branco)
# Spectra block
spectra <- as.matrix(result[, grep("^wl_", colnames(result))])
# 1. Sample Data (Replace this with your actual data)
wavelengths <- seq(340, 850, length.out = 288) # X-axis values
# Matrix where each COLUMN is a spectrum (10 spectra, 100 bands each)
spectral_data <- t(spectra)
spectral_data_correct <- spectral_data/spectra_white_black$Branco
spectra_reflectance <- 1 / spectral_data_correct
pseudo_absorbance_spectra <- log10(1 / (spectral_data_correct + 1e-6))
# 2. Plotting the Spectra
matplot(
  x = wavelengths, 
  y = spectral_data_correct, 
  type = "l", 
  lty = 1, 
  col = 1:ncol(spectral_data), 
  xlab = "Wavelength (nm)", 
  ylab = "Intensity/Reflectance/Absorbance", 
  main = "All Spectra Overview"
)

# 3. Add Legend (Optional but recommended)
legend("topright", 
       legend = paste("Spectrum", 1:ncol(spectral_data)), 
       col = 1:ncol(spectral_data), 
       lty = 1, 
       cex = 0.8
)

spectra_fruit <- as.matrix(result_fruit[, grep("^wl_", colnames(result_fruit))])
# Genes block
genes <- as.matrix(result[,c("MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                            "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
                            "LBCY (Cq)","CHLG (Cq)")]) #"GAPDH (Cq)"

genes_vine <- as.matrix(result[,c("RCA (Cq)",
                                  "LBCY (Cq)","CHLG (Cq)")])
genes_vine_lagged <- as.matrix(result[,c("RCA (Cq)_lag7d",
                                  "LBCY (Cq)_lag7d","CHLG (Cq)_lag7d")])
genes <- genes[, !colnames(genes) %in% "GAPDH (Cq)"] #retirar housekeeping

genes_fruit <- as.matrix(result_fruit[,c("UFGT (Cq)", "PAL (Cq)","CHS (Cq)" ,"UFGT1 (Cq)","PAL1 (Cq)")])
genes_fruit <- genes_fruit[, !colnames(genes_fruit) %in% "PAL1 (Cq)"]
dim(genes_fruit)
# Convert Cq -> expression
genes_expr <- 2^(-genes_vine)
genes_expr_lag <- 2^(-genes_vine_lagged)
#genes_expr <- 2^(-genes_vine)
genes_expr_fruit <- 2^(-genes_fruit)

genes_expr <- genes %>%
  mutate(across(-Date, ~ 2^(-.)))

# (optional) keep column names the same, but remove the " (Cq)" suffix
colnames(genes_expr_fruit) <- gsub(" \\(Cq\\)", "", colnames(genes_expr_fruit))
colnames(genes_expr) <- gsub(" \\(Cq\\)", "", colnames(genes_expr))
colnames(genes_expr_lag) <- gsub(" \\(Cq\\)", "", colnames(genes_expr_lag))


# Metabolites block
metabolites <- as.matrix(result[,c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                              "Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                              "Sugar__mg_gFM_","Anthocyanins__mg_gFM_")])

metabolites_lagged <- as.matrix(result[,c("Carotenoids__mg_gFM__lag7d","Chlorophyll_a__mg_gFM__lag7d","Chlorophyll_b__mg_gFM__lag7d",
                                   "Phenols__g_GAE_mL__lag7d","ROS__mmol_gFM__lag7d","ROS_O2__ABS_gFM__lag7d","Starch__mg_gFM__lag7d",
                                   "Sugar__mg_gFM__lag7d","Anthocyanins__mg_gFM__lag7d")])

metabolites_fruit <- as.matrix(result_fruit[,c("Alpha_Amino_Nitrogen__mg_L_","Ammoniacal_Nitrogen__mg_L_",
                                               "Assimilable_Nitrogen__mg_L_","Brix_lab__%_","Brix_tom__%_")])

# Find sample(s) with unusually high values
which(metabolites_fruit[, "Assimilable_Nitrogen__mg_L_"] > 1000)


# Remove outliers - set the value to NA
metabolites_fruit_noout <- metabolites_fruit
metabolites_fruit_noout[metabolites_fruit_noout[, "Assimilable_Nitrogen__mg_L_"] > 1000, "Assimilable_Nitrogen__mg_L_"] <- NA


rownames(genes_vine) <- rownames(result)
rownames(genes_expr) <- rownames(result)
rownames(genes_expr_lag) <- rownames(result)
rownames(metabolites) <- rownames(result)
rownames(metabolites_lagged) <- rownames(result)
rownames(spectra) <- rownames(result)

rownames(genes_fruit) <- rownames(result_fruit)
rownames(genes_expr_fruit) <- rownames(result_fruit)
rownames(metabolites_fruit) <- rownames(result_fruit)
rownames(spectra_fruit) <- rownames(result_fruit)

##### CHECK VARIANCE OF EACH BLOCK #########

## LEAF
# Variance per feature
var_gene_expr <- apply(genes_expr, 2, var, na.rm = TRUE)
var_metabolites <- apply(metabolites, 2, var, na.rm = TRUE)
var_spectra_pcs <- apply(spectra_pcs, 2, var, na.rm = TRUE)
var_covs <- apply(cov_block_select, 2, var, na.rm = TRUE)
# Average variance per block
mean(var_gene_expr)
mean(var_metabolites)
mean(var_spectra_pcs)
mean(var_covs)
## Total variance per block
# Total variance per block
tot_var_genes_expr <- sum(var_gene_expr)
tot_var_metabolites <- sum(var_metabolites)
tot_var_spectra_pcs <- sum(var_spectra_pcs)
tot_var_covs <- sum(var_covs)

# Compare
tot_var_genes_expr
tot_var_metabolites
tot_var_spectra_pcs
tot_var_covs

boxplot(list(Genes = var_gene_expr, Metabolites = var_metabolites, 
             Spectra_pcs = var_spectra_pcs, Covariates = var_covs ),
        ylab = "Variance per feature",
        main = "Variance distribution across blocks")


boxplot(list(Genes = tot_var_genes_expr, Metabolites = tot_var_metabolites, 
             Spectra_pcs = tot_var_spectra_pcs, Covariates = tot_var_covs ),
        ylab = "Total block variance",
        main = "Variance distribution across blocks")

# Variance per feature AFTER Z-score
var_gene_expr_scaled <- apply(genes_scaled, 2, var, na.rm = TRUE)
var_metabolites_scaled <- apply(metabs_scaled, 2, var, na.rm = TRUE)
var_spectra_pcs_scaled <- apply(spectra_pcs_scaled, 2, var, na.rm = TRUE)
var_cov_scaled <- apply(cov_block_scaled, 2, var, na.rm = TRUE)

boxplot(list(Genes = var_gene_expr_scaled, Metabolites = var_metabolites_scaled,
             Spectra_pcs = spectra_pcs_scaled, Covariates = var_cov_scaled ),
        ylab = "Variance per feature",
        main = "Variance distribution across blocks (after z-score)")


## FRUIT
# Variance per feature
var_gene_expr_fruit <- apply(genes_expr_fruit, 2, var, na.rm = TRUE)
var_metabolites_fruit <- apply(metabolites_fruit, 2, var, na.rm = TRUE)
var_spectra_pcs_fruit <- apply(spectra_fruit_pcs, 2, var, na.rm = TRUE)
# Average variance per block
mean(var_gene_expr_fruit)
mean(var_metabolites_fruit)
mean(var_spectra_pcs_fruit)

## Total variance per block
# Total variance per block
tot_var_genes_expr_fruit <- sum(var_gene_expr_fruit)
tot_var_metabolites_fruit <- sum(var_metabolites_fruit)
tot_var_spectra_pcs_fruit <- sum(var_spectra_pcs_fruit)
                                 
# Compare
tot_var_genes_expr_fruit
tot_var_metabolites_fruit
tot_var_spectra_pcs_fruit

boxplot(list(Genes = var_gene_expr_fruit, Metabolites = var_metabolites_fruit, Spectra_pcs = var_spectra_pcs_fruit ),
        ylab = "Variance per feature",
        main = "Variance distribution across blocks")

# Variance per feature AFTER Z-score
var_gene_expr_fruit_scaled <- apply(genes_fruit_scaled, 2, var, na.rm = TRUE)
var_metabolites_fruit_scaled <- apply(metabs_fruit_scaled, 2, var, na.rm = TRUE)
var_spectra_pcs_fruit_scaled <- apply(spectra_fruit_scaled, 2, var, na.rm = TRUE)

boxplot(list(Genes = var_gene_expr_fruit_scaled, Metabolites = var_metabolites_fruit_scaled, Spectra_pcs = spectra_fruit_scaled ),
        ylab = "Variance per feature",
        main = "Variance distribution across blocks (after z-score)")

# Get names of columns that are all NA
na_cols <- names(genes)[sapply(genes, function(x) all(is.na(x)))]

na_cols

# Check if there are constant genes
# Compute standard deviation per gene
gene_sd <- apply(genes, 2, sd, na.rm = TRUE)

# See which genes are exactly constant
constant_genes <- names(gene_sd[gene_sd == 0])

# Summary
length(constant_genes)       # number of constant genes
constant_genes               # their names
summary(gene_sd)             # distribution of variability across genes to check if they have ow variance

# --- PCA on spectra to check redundancy --- #
# If first few PCs explain most variance, strong multicollinearity is present
pca_spec <- prcomp(spectra, center=TRUE, scale.=TRUE)
cumvar <- cumsum(summary(pca_spec)$importance[2,])
plot(cumvar, type="l", main="Cumulative variance explained (Spectra PCA)",
     xlab="PCs", ylab="Cumulative variance")

# block scaling function
block_total_variance_scale <- function(mat){
  #mat <- as.matrix(mat)
  mat <- scale(mat, center=TRUE, scale=TRUE)    # scale features
  total_var <- sum(apply(mat, 2, var, na.rm=TRUE))
  if (total_var > 0) mat <- mat / sqrt(total_var)   # normalize total block variance
  return(mat)
}

block_scale <- function(mat){
  mat <- scale(mat, center=TRUE, scale=TRUE)
  total_var <- sum(apply(mat, 2, var, na.rm=TRUE))
  mat / sqrt(total_var)   # normalize total block variance
}

block_total_variance_scale <- function(mat){
  mat <- scale(mat, center=TRUE, scale=TRUE)  # z-score each feature
  total_var <- sum(apply(mat, 2, var, na.rm=TRUE))
  if (total_var > 0) mat <- mat / sqrt(total_var)  # scale total block variance to 1
  return(mat)
}

eff <- 2         # PCR efficiency (use assay-specific if known)
eps <- 1e-6      # pseudocount
# Convert to relative expression
genes_linear <- eff^(-genes)

# Log-transform for variance stabilization
genes_log <- log2(genes_linear + eps)

### Prepare two alternative datasets ###

# 1) Scaled version (all blocks balanced)
genes_scaled   <- block_total_variance_scale(genes)
metabs_scaled  <- block_total_variance_scale(metabolites)
spectra_scaled <- block_total_variance_scale(spectra)

spectra_scaled_pcs <- block_total_variance_scale(spectra_pcs)
total_var <- c(
  genes     = sum(apply(genes_scaled,   2, var, na.rm=TRUE)),
  metabs    = sum(apply(metabs_scaled,  2, var, na.rm=TRUE)),
  spectra   = sum(apply(spectra_scaled_pcs, 2, var, na.rm=TRUE))
)
print(total_var)

genes_clean <- genes_scaled[complete.cases(genes_scaled), ]
pca_genes <- prcomp(genes_clean, center = FALSE, scale. = FALSE)
genes_pcs <- pca_genes$x[, 1:3]  # take first 3 PCs
##Check if any genes became effectively constant after the transformation:
gene_sd_post <- apply(genes_scaled, 2, sd, na.rm = TRUE)
summary(gene_sd_post)

#genes_scaled   <- scale(genes, center=TRUE, scale=TRUE)
#metabs_scaled  <- scale(metabolites, center=TRUE, scale=TRUE)
#spectra_scaled <- scale(spectra, center=TRUE, scale=TRUE)

coverage_genes <- 1 - colMeans(is.na(genes))
scale_factor <- mean(coverage_genes)  # downweight block if coverage is low
genes_scaled2 <- genes_scaled * scale_factor

data_scaled <- list(
  gene_expression = t(genes_scaled),
  metabolites     = t(metabs_scaled),
  spectra         = t(spectra_scaled)
)

# 2) Spectra reduced to PCs (remove redundancy
# Keep PCs explaining ~90% variance
pca_spec <- prcomp(spectra, center=TRUE, scale.=TRUE)
cumvar <- cumsum(summary(pca_spec)$importance[2,])
k <- which(cumvar >= 0.99)[1]
spectra_pcs <- pca_spec$x[, 1:k, drop=FALSE]

pc_scores_df <- as.data.frame(spectra_pcs)
correlation_matrix <- cor(pc_scores_df)

print("Correlation Matrix of Principal Component Scores:")
print(correlation_matrix)


# 1. Extract the standard deviations (sdev) and square them to get the variance of each PC
pc_variance <- pca_spec$sdev^2

# 2. Calculate the proportion of variance explained by each PC
variance_explained_ratio <- pc_variance / sum(pc_variance)

# 3. Calculate the Cumulative Variance Explained (CVE)
cumulative_variance_explained <- cumsum(variance_explained_ratio)
# 4. Create a DataFrame for easy viewing
variance_summary_df <- data.frame(
  PC = 1:length(pc_variance),
  Variance_Ratio = variance_explained_ratio,
  Cumulative_Variance = cumulative_variance_explained
)

# 5. Print the summary for the first 10 PCs (or more if needed)
print(head(variance_summary_df, 10))

pca_spec_fruit <- prcomp(spectra_fruit, center=TRUE, scale.=TRUE)
cumvar_fruit <- cumsum(summary(pca_spec_fruit)$importance[2,])
g <- which(cumvar_fruit >= 0.99)[1]
spectra_fruit_pcs <- pca_spec_fruit$x[, 1:g, drop=FALSE]

# Reconstruct from first g PCs
#spectra_fruit_recon <- scale(
 # spectra_fruit_pcs %*% t(pca_spec_fruit$rotation[, 1:g]),
  #center = -pca_spec_fruit$center,
  #scale  = 1 / pca_spec_fruit$scale
#)

#all contributions for the first g PCs
loadings <- pca_spec_fruit$rotation[, 1:g]

# Top contributing features for PC1
pc1_loadings <- pca_spec_fruit$rotation[, 1]
pc2_loadings <- pca_spec_fruit$rotation[, 2]
pc3_loadings <- pca_spec_fruit$rotation[, 3]
pc4_loadings <- pca_spec_fruit$rotation[, 4]
pc5_loadings <- pca_spec_fruit$rotation[, 5]
pc6_loadings <- pca_spec_fruit$rotation[, 6]
pc7_loadings <- pca_spec_fruit$rotation[, 7]
pc8_loadings <- pca_spec_fruit$rotation[, 8]
pc9_loadings <- pca_spec_fruit$rotation[, 9]

sort(abs(pc1_loadings), decreasing = TRUE)[1:100]  # top 10
sort(abs(pc1_loadings), decreasing = TRUE)
sort(abs(pc2_loadings), decreasing = TRUE)[1:50]
sort(abs(pc3_loadings), decreasing = TRUE)
sort(abs(pc4_loadings), decreasing = TRUE)[1:100]
sort(abs(pc5_loadings), decreasing = TRUE)
sort(abs(pc6_loadings), decreasing = TRUE)[1:50]
sort(abs(pc7_loadings), decreasing = TRUE)
sort(abs(pc8_loadings), decreasing = TRUE)
sort(abs(pc9_loadings), decreasing = TRUE)

barplot(pca_spec_fruit$rotation[,9],
        las = 2, main = "Feature contributions to PC9")



#### Leaf
loadings <- pca_spec$rotation[, 1:k]

# Top contributing features for PC1
pc1_loadings <- pca_spec$rotation[, 1]
pc2_loadings <- pca_spec$rotation[, 2]
pc3_loadings <- pca_spec$rotation[, 3]
pc4_loadings <- pca_spec$rotation[, 4]
pc5_loadings <- pca_spec$rotation[, 5]
pc6_loadings <- pca_spec$rotation[, 6]
pc7_loadings <- pca_spec$rotation[, 7]
pc8_loadings <- pca_spec$rotation[, 8]

indices_ordenados <- pc1_loadings[order(abs(pc1_loadings), decreasing = TRUE)]
indices_ordenados[1:100]
indices_ordenados <- order(abs(pc8_loadings), decreasing = TRUE)

# 2. Aplicar os índices ordenados ao vetor original (mantendo o sinal)
pc4_loadings_ordenados <- pc8_loadings[indices_ordenados]

# 3. Mostrar os 100 maiores loadings (em valor absoluto) com o sinal original
pc4_loadings_ordenados[1:100]

sort(abs(pc1_loadings), decreasing = TRUE)[1:100]  # top 10
sort(abs(pc1_loadings), decreasing = TRUE)
sort(abs(pc2_loadings), decreasing = TRUE)[1:100]
sort(abs(pc3_loadings), decreasing = TRUE)[1:100]
sort(abs(pc4_loadings), decreasing = TRUE)[1:100]
sort(abs(pc5_loadings), decreasing = TRUE)[1:100]
sort(abs(pc6_loadings), decreasing = TRUE)
sort(abs(pc7_loadings), decreasing = TRUE)
sort(abs(pc8_loadings), decreasing = TRUE)


barplot(pca_spec$rotation[,8],
        las = 2, main = "Feature contributions to PC8")



######## Do it for all PCs facet wrap ###
# Get loadings for the first g PCs
# Get loadings for the first g PCs
loadings <- pca_spec_fruit$rotation[, 1:g]
loadings <- pca_spec$rotation[, 1:k]

# Extract numeric wavelengths from row names (assuming they are like "wl_770.0348")
wavelengths <- as.numeric(sub("wl_", "", rownames(loadings)))

# Convert to data frame
loadings_df <- as.data.frame(loadings)
loadings_df$wavelength <- wavelengths

# Reshape into long format
loadings_long <- melt(loadings_df, id.vars = "wavelength",
                      variable.name = "PC",
                      value.name = "loading")

# Plot as line plot across wavelength
ggplot(loadings_long, aes(x = wavelength, y = loading)) +
  geom_line() +
  facet_wrap(~ PC, scales = "free_y") +
  labs(title = "Spectral loadings for PCs",
       x = "Wavelength (nm)",
       y = "Loading") +
  theme_minimal()
##### do PCA on genes as well)
i <- which(cumvar >= 0.99)[1]
pca_genes <- prcomp(spectra, center=TRUE, scale.=TRUE)
genes_pcs <- pca_genes$x[, 1:k, drop=FALSE]

anyNA(pca_genes$x) 

# Variance of each PC
pc_var <- pca_genes$sdev^2

# Proportion of variance explained
pc_var_explained <- pc_var / sum(pc_var)

pc_var_explained

#k <- min(10, ncol(spectra_scaled))
#spectra_pcs <- pca_spec$x[, 1:k, drop=FALSE]
#pca_spec <- prcomp(spectra_scaled, center=TRUE, scale.=TRUE)


# Scaled model
# Step 1. Wrap your data list into a MOFA object
MOFA_scaled <- create_mofa(data_scaled)

# Step 2. Set options
data_opts  <- get_default_data_options(MOFA_scaled)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(MOFA_scaled)
model_opts$num_factors <- 15

train_opts <- get_default_training_options(MOFA_scaled)
train_opts$maxiter <- 1000
train_opts$convergence_mode <- "fast"

# (optional) Explicit likelihoods
#likelihoods <- c(
 # "gene_expression" = "gaussian",
  #"metabolites"     = "gaussian",
  #"spectra"         = "gaussian"
#)

# Step 3. Prepare the model
MOFA_scaled <- prepare_mofa(
  object            = MOFA_scaled,
  data_options      = data_opts,
  model_options     = model_opts,
  training_options  = train_opts,
  #likelihoods       = likelihoods
)


# Step 4. Train the model
MOFA_scaled <- run_mofa(MOFA_scaled, use_basilisk = TRUE)



#preproc_genes <- caret::preProcess(genes, method = "range")   # range = min–max scaling
#genes_mm <- predict(preproc_genes, genes_scaled)

#preproc_metabs <- caret::preProcess(metabolites, method = "range")   # range = min–max scaling
#metabs_mm <- predict(preproc_metabs, metabs_scaled)

#preproc_spectra_pcs <- caret::preProcess(spectra_pcs, method = "range")   # range = min–max scaling
#spectra_pcs__mm <- predict(preproc_spectra_pcs, spectra_pcs)

minmax_scale <- function(mat){
  apply(mat, 2, function(x){
    rng <- range(x, na.rm = TRUE)
    if(diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  })
}

genes_scaled_mm <- minmax_scale(genes)
metabs_scaled_mm <- minmax_scale(metabolites)
spectra_pcs_mm <-  minmax_scale(spectra_pcs)


# only z-score scaling
genes_scaled <- scale(genes_expr)
genes_lag_scaled <- scale(genes_expr_lag)#/sqrt(ncol(genes_expr))      # center + unit variance
metabs_scaled <- scale(metabolites)#/sqrt(ncol(metabolites))
metabs_lag_scaled <- scale(metabolites_lagged)
spectra_pcs_scaled <- scale(spectra_pcs)#/sqrt(ncol(spectra_pcs)) # same for spectra
#cov_block_scaled <- cov_block_scaled/sqrt(ncol(cov_block_select_cat_num))

genes_scaled <- genes_scaled * sqrt(1 / sum(apply(genes_scaled, 2, var)))
metabs_scaled  <- metabs_scaled * sqrt(1 / sum(apply(metabs_scaled , 2, var)))
spectra_pcs2  <- spectra_pcs * sqrt(1 / sum(apply(spectra_pcs , 2, var)))


# weigh each block (or feature) by the proportion of non-missing values.
# Compute coverage per gene (proportion of non-missing values)
coverage_genes <- 1 - colMeans(is.na(genes))

# Scale each gene individually by its coverage
genes_scaled2 <- sweep(genes_scaled, 2, coverage_genes, `*`)

coverage_metabs <- 1- colMeans(is.na(metabolites))
metabs_scaled2 <-  sweep(metabs_scaled, 2, coverage_metabs, `*`)

coverage_spectra_pcs <- 1- colMeans(is.na(spectra_pcs))
spectra_scaled2 <- sweep(spectra_pcs, 2, coverage_spectra_pcs, `*`)

# Optionally rescale to match variance of other blocks
var_genes <- mean(apply(genes, 2, var, na.rm = TRUE))
var_spectra <- mean(apply(spectra_pcs, 2, var))
var_metabs <- sum(apply(metabs_scaled, 2, var,na.rm = TRUE)) 


target_var <- mean(c(var_genes, var_spectra, var_metabs))  # average total variance


# boost genes by 10-20% so they dominate early factors????
scale_genes <- sqrt((target_var) / var_genes)
scale_spectra <- sqrt(target_var / var_spectra)
scale_metabs <- sqrt(target_var / var_metabs)

#genes_scaled2 <- genes_scaled * sqrt(var_spectra / var_genes)
genes_scaled2 <- genes_scaled * scale_genes
spectra_scaled2 <- spectra_pcs * scale_spectra
metabs_scaled2 <- metabs_scaled * scale_metabs

################
### FRUIT ######
genes_fruit_scaled <- scale(genes_expr_fruit) 
spectra_fruit_scaled <- scale(spectra_fruit_pcs)
metabs_fruit_scaled <- scale(metabolites_fruit)
################

# PCA-reduced model

apply(genes, 2, var,na.rm=TRUE)
apply(genes_scaled, 2, var,na.rm=TRUE)
apply(spectra_pcs_scaled, 2, var)

sum(apply(genes_scaled, 2, var, na.rm=TRUE)) 
sum(apply(metabs_scaled, 2, var, na.rm=TRUE)) 
sum(apply(spectra_pcs_scaled, 2, var, na.rm=TRUE)) 


# Block-wise scaling: divide each block by sqrt(total variance of the block)
genes_block_scaled <- genes_scaled / sqrt(sum(apply(genes_scaled, 2, var, na.rm=TRUE)))
metabs_block_scaled <- metabs_scaled / sqrt(sum(apply(metabs_scaled, 2, var, na.rm=TRUE)))
spectra_block_scaled_pcs <- spectra_pcs_scaled / sqrt(sum(apply(spectra_pcs_scaled, 2, var, na.rm=TRUE)))


sum(apply(genes_block_scaled , 2, var, na.rm=TRUE)) 
sum(apply(metabs_block_scaled, 2, var, na.rm=TRUE)) 
sum(apply(spectra_block_scaled_pcs, 2, var, na.rm=TRUE)) 


data_pcs <- list(
  gene_expression = t(genes_scaled),
  metabolites     = t(metabs_scaled),
  spectra         = t(spectra_pcs)
)
# Step 1. Wrap your data list into a MOFA object
MOFA_pcs <- create_mofa(data_pcs)

# Step 2. Set options
data_opts  <- get_default_data_options(MOFA_pcs)
data_opts$scale_views <- FALSE # disable automatic scaling
#data_opts$block_weights <- c(genes = 1, metabolites = 2, spectra = 2)  # upweight smaller blocks

model_opts <- get_default_model_options(MOFA_pcs)
model_opts$num_factors <- 12
model_opts$spikeslab_weights$use <- FALSE
# Enable sparsity on weights
#model_opts$spikeslab_weights$use <- TRUE
#model_opts$spikeslab_weights$alpha0 <- list(genes = 1e-16, spectra = 1e-2)
#model_opts$spikeslab_weights$beta0 <- list(genes = 10, spectra = 1)

#model_opts$sparsity <- list(factors = 0.5, views = NULL)


train_opts <- get_default_training_options(MOFA_pcs)
train_opts$maxiter <- 2000
train_opts$convergence_mode <- "medium"

# (optional) Explicit likelihoods
#likelihoods <- c(
# "gene_expression" = "gaussian",
#"metabolites"     = "gaussian",
#"spectra"         = "gaussian"
#)

# Step 3. Prepare the model
MOFA_pcs <- prepare_mofa(
  object            = MOFA_pcs,
  data_options      = data_opts,
  model_options     = model_opts,
  training_options  = train_opts
  #likelihoods       = likelihoods
)

# Step 4. Train the model
MOFA_pcs <- run_mofa(MOFA_pcs, use_basilisk = TRUE)

# Extract factor values for all samples
# Get factors (list of length = 1, because factors are shared across views)
factors_list <- get_factors(MOFA_pcs, factors = "all")
# Check names
names(factors_list)
Z <- factors_list$group1  
dim(Z)
class(Z)

#### Diagnostics for FRUIT #####
#Compute total expression per sample: assumes matrix is already in expression scale (2^-ΔCt or 2^-Cq
# Sum across all genes per sample
total_expr <- rowSums(genes_expr, na.rm = TRUE)
total_expr <- rowSums(genes_expr_fruit, na.rm = TRUE)

total_expr_metabs <- rowSums(metabolites, na.rm = TRUE)
total_expr_metabs <- rowSums(metabolites_fruit, na.rm = TRUE)

# Or mean expression per sample
mean_expr <- rowMeans(genes_expr, na.rm = TRUE)
mean_expr <- rowMeans(genes_expr_fruit, na.rm = TRUE)

mean_expr_metab <- rowMeans(metabolites,na.rm = TRUE)
mean_expr_metab <- rowMeans(metabolites_fruit,na.rm = TRUE)

# Basic barplot
barplot(total_expr, main = "Total expression per sample", xlab = "Sample", ylab = "Sum of expression")
barplot(total_expr_metabs, main = "Concentration per sample", xlab = "Sample", ylab = "Sum of expression")

# Or boxplot
boxplot(genes_expr_fruit, main = "Expression distribution per gene", las=2,cex.axis=0.8)

png("metabolite_boxplot_noout.png", width=1600, height=1000, res=150)

# expand bottom margin so labels fit
par(mar=c(12, 4, 4, 2))

boxplot(metabolites_fruit_noout, las=2, cex.axis=0.8,
        main="Metabolite distributions - no outlier")

dev.off()

df_long_metabolites <- melt(as.data.frame(metabolites))

df_long_metabolites_fruit <- melt(as.data.frame(metabolites_fruit))


p <- ggplot(df_long_metabolites, aes(x=variable, y=value)) +
  geom_boxplot() + ggtitle("Concentration distribution per metabolite") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("Boxplot_conc_metabs_fruit.pdf", p, width=20, height=12, dpi=300)


# --- 2. Prepare covariates -------------------------------------------------
# scale numeric covariates (z-score)
#num_vars <- c("mean_P","mean_Tmed","mean_Tmax","mean_HRmed","Water_P",
              #"month_sin","month_cos","day_sin","day_cos","yday")
#num_vars <- result2_noOut_starch %>% select(c(333), 337:379)
#num_idx <- c(333,337:379)
num_idx <- c(337:381)
num_idx_fruit <- c(319:363) 
#result2_noOut_starch[num_vars] <- scale(result2_noOut_starch[num_vars])

# 2C: one-hot encode categorical variables (for regression)
# We'll make a design matrix for regressCovariates (no intercept)
cat_vars <- c("Test_site","Cultivar","year","Irrigation")
# design mat: samples x (sum levels-#levels) binary columns
cov_cat_dummy <- model.matrix(~ 0 + Test_site + Cultivar + year + Irrigation, 
                           data = result2_noOut_starch2,
                           contrasts.arg = lapply(result[c("Test_site","Cultivar","year","Irrigation")],
                                                  function(x) contrasts(x, contrasts = FALSE)))

cov_cat_dummy <- model.matrix(~ 0 + year + Irrigation, 
                              data = result,
                              contrasts.arg = lapply(result[c("year","Irrigation")],
                                                     function(x) contrasts(x, contrasts = FALSE)))

cov_cat_dummy_fruit <-  model.matrix(~ 0 + Test_site + Cultivar + year + Irrigation, 
                               data = result_fruit2,
                               contrasts.arg = lapply(result_fruit2[c("Test_site","Cultivar","year","Irrigation")],
                                                      function(x) contrasts(x, contrasts = FALSE)))
# Make sure they are factors or numeric
result[cat_vars] <- lapply(result[cat_vars], as.factor)
result[num_idx] <- lapply(result[num_idx], as.numeric)
result_fruit2[num_idx_fruit] <- lapply(result_fruit2[num_idx_fruit], as.numeric)
### MOFA with covariates ##########
#data_cov <- list(
#  gene_expression = t(genes_scaled),
#  metabolites     = t(metabs_scaled),
#  spectra         = t(spectra_pcs)
#)


####Include numeric covariates
 ##1 -  include as a separate view (covariate block)

cov_block <- as.matrix(result[, num_idx])# samples x numeric covariates
#cov_block <- cbind(covariates_dummy,cov_block)
rownames(cov_block) <- rownames(result)
colnames(cov_block) <- colnames(cov_block) %>%
  gsub("º", "deg", .) %>%       # replace degree symbol
  gsub("\\.", "_", .) %>%       # replace dots with underscore
  gsub("-", "_", .) %>%         # replace dashes with underscore
  gsub(" ", "_", .)              # replace spaces with underscore

cov_block_fruit <- as.matrix(result_fruit2[, num_idx_fruit])# samples x numeric covariates
colnames(cov_block_fruit) <- colnames(cov_block_fruit) %>%
  gsub("º", "deg", .) %>%       # replace degree symbol
  gsub("\\.", "_", .) %>%       # replace dots with underscore
  gsub("-", "_", .) %>%         # replace dashes with underscore
  gsub(" ", "_", .)              # replace spaces with underscore


cov_block_select <- subset(cov_block, 
                           select = c(yday, month_sin, Water_P,P__mm__roll7d_mean,
                                      P__mm__lag1d,Tmax__degC__roll7d_mean,Tmax__degC__lag1d,
                                      HRmed_____roll7d_mean, Water_P_lag7daytol))

cov_block_select_fruit <- subset(cov_block_fruit, 
                           select = c(yday, month_sin, Water_P,P__mm__roll7d_mean,
                                      P__mm__lag1d,Tmax__degC__roll7d_mean,Tmax__degC__lag1d,
                                      HRmed_____roll7d_mean, Water_P_lag7daytol))

cov_block_select_cat_num <- cbind(cov_block_select,cov_cat_dummy)
cov_block_select_cat_num_fruit <- cbind(cov_block_select_fruit,cov_cat_dummy_fruit)
#colnames(cov_block) <- colnames(result)[num_idx] 
#cov_block <- apply(cov_block, 2, function(x) as.numeric(as.character(x)))
cov_block_scaled <- scale(cov_block_select_cat_num)
cov_block_fruit_scaled <- scale(cov_block_select_cat_num_fruit)

sapply(cov_block_fruit, class)
sapply(cov_cat_dummy, class)
#cov_block_scaled <- cbind(design_mat, cov_block_scaled)

residualize_matrix <- function(mat, design_df) {
  # mat: samples x features
  # design_df: data.frame with categorical covariates
  mat_resid <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  rownames(mat_resid) <- rownames(mat)
  colnames(mat_resid) <- colnames(mat)
  
  for (i in seq_len(ncol(mat))) {
    fit <- lm( mat[, i] ~ . , data = design_df,qr=True)
    mat_resid[, i] <- residuals(fit)
  }
  
  return(mat_resid)
}

colnames(cov_block_scaled) <- colnames(cov_block_scaled) %>%
  gsub("º", "deg", .) %>%       # replace degree symbol
  gsub("\\.", "_", .) %>%       # replace dots with underscore
  gsub("-", "_", .) %>%         # replace dashes with underscore
  gsub(" ", "_", .)              # replace spaces with underscore

#cov_block_scaled <- subset(cov_block_scaled, select = -c(P__mm__lag7d,Irrigation_bin))
cov_block_scaled <- subset(cov_block_scaled, 
                           select = c(yday, month_sin, Water_P,P__mm__roll7d_mean,
                                      P__mm__lag1d,Tmax__degC__roll7d_mean,Tmax__degC__lag1d,
                                      HRmed_____roll7d_mean))

#cov_block_scaled <- cov_block_scaled/sqrt(ncol(cov_block_scaled))
cov_block_scaled2  <- cov_block_scaled * sqrt(1 / sum(apply(cov_block_scaled , 2, var)))

# Compute correlation matrix
cov_corr <- cor(cov_block_select_cat_num, use = "pairwise.complete.obs")  # handles NAs
library(corrplot)
# Optional: visualize as heatmap
corrplot(cov_corr, method = "color", type = "upper", na.label = " ",
         addCoef.col = "black",number.cex = 0.2 )



# Define the columns that need to be removed
cols_to_remove <- c("ROS__mmol_gFM_","ROS__mmol_gFM__lag7d")
columns_to_remove <- c("ROS__mmol_gFM_","Phenols__g_GAE_mL_")
columns_to_remove2 <- c("Phenols__g_GAE_mL__lag7d","ROS__mmol_gFM__lag7d")
cols_to_remove <- c("Test_siteGreenhouse Porto")
# Create a clean data frame by removing these columns
metabs_scaled <- metabs_scaled[, !colnames(metabs_scaled) %in% columns_to_remove]
metabs_scaled <- metabs_scaled[, !colnames(metabs_scaled) %in% cols_to_remove]
metabs_lag_scaled <- metabs_lag_scaled[, !colnames(metabs_lag_scaled) %in% cols_to_remove]
cov_block_scaled <- cov_block_scaled[, !colnames(cov_block_scaled) %in% cols_to_remove]


views <- list(
  genes_expression = t(genes_scaled),
  metabolites = t(metabs_scaled),
  spectra = t(spectra_pcs),
  others = t(cov_block_scaled)
)

views <- list(
  metabolites = t(metabs_scaled),
  metabolites_lagged = t(metabs_lag_scaled),
  spectra = t(spectra_pcs),
  others = t(cov_block_scaled)
)

views <- list(
  genes_expression = t(genes_fruit_scaled),
  metabolites = t(metabs_fruit_scaled),
  spectra = t(spectra_fruit_pcs),
  others = t(cov_block_fruit_scaled)
)

sapply(views, function(v) mean(is.na(v)))
sapply(views, function(v) mean(!is.na(v)))

#vars <- sapply(views, function(x) sum(apply(x, 1, var, na.rm = TRUE)))
#weights <- 1 / sqrt(vars)  # inverso da raiz da variância total
#views_weighted <- mapply(function(x, w) x * w, views, weights, SIMPLIFY = FALSE)

# Scale each view by the fraction of observed values
views_scaled <- lapply(views, function(v) {
  frac_obs <- mean(!is.na(v))
  v * frac_obs  # or v * sqrt(frac_obs) if you want a softer adjustment
})
#views$genes_expression <- views$genes_expression * 0.009

sapply(views_weighted, function(v) mean(is.na(v)))
sapply(views_weighted, function(v) mean(!is.na(v)))

scale_views_by_features_missingness <- function(views,
                                                alpha = 1,       # exponent for completeness (frac_obs^alpha)
                                                beta  = 0.5,     # exponent for n_features (n_feat^beta)
                                                per_feature = FALSE, # if TRUE also reweight each feature by its own frac_obs_row^gamma
                                                gamma = 1,       # exponent for per-feature completeness
                                                keep_na = TRUE,
                                                eps = 1e-8) {
  # views: named list of matrices/data.frames (features x samples)
  scaled <- lapply(views, function(v) {
    mat <- as.matrix(v)
    # effective number of features: features that have at least one observed value
    frac_obs_view <- mean(!is.na(mat))
    n_feat_eff <- sum(rowMeans(!is.na(mat)) > 0)  # count features with any non-NA
    if (n_feat_eff == 0) {
      warning("A view has zero effective features (all NA). Returning original.")
      return(mat)
    }
    
    # compute view-level weight
    # weight increases with completeness, decreases with number of features
    view_weight <- (frac_obs_view + eps)^alpha / ( (n_feat_eff + eps)^beta )
    
    # apply view-level weight
    mat_w <- mat * view_weight
    
    if (per_feature) {
      # per-feature completeness
      frac_obs_feat <- rowMeans(!is.na(mat))
      # avoid zero weights
      feat_weight <- (frac_obs_feat + eps)^gamma
      # multiply each row by its feature weight
      # keep NA positions as NA (unless keep_na==FALSE)
      mat_w <- mat_w * feat_weight
    }
    
    if (!keep_na) {
      # optionally replace NAs by 0 to stabilise algorithms that don't like NA
      mat_w[is.na(mat_w)] <- 0
    }
    mat_w
  })
  
  names(scaled) <- names(views)
  scaled
}

views_scaled_default <- scale_views_by_features_missingness(views, alpha = 1, beta = 0.5)
views_scaled_mild <- scale_views_by_features_missingness(views, alpha = 0.5, beta = 0.5,
                                                         per_feature = TRUE, gamma = 0.5)


# aggressive downweighting of large views
views_scaled_aggressive <- scale_views_by_features_missingness(views, alpha = 1, beta = 1)

# Simple linear weighting by number of features and fraction observed
scale_views_linear <- function(views) {
  lapply(views, function(v) {
    mat <- as.matrix(v)
    n_feat <- nrow(mat)
    frac_obs <- mean(!is.na(mat))
    #print(frac_obs)
    weight <- frac_obs / n_feat
    print(weight)
    mat <- mat * weight
    return mat
  })
}
views_scaled <- scale_views_linear(views)

views_scaled <- lapply(views, function(v) {
   v * (mean(!is.na(v)) / sqrt(nrow(v)))
}) #as.matrix(v) *
#views_scaled <- lapply(views, function(x) x / sqrt(sum(apply(x, 1, var), na.rm = TRUE)))

# Function to scale each view by sqrt(fraction observed) / sqrt(number of features)
scale_views_soft <- function(views) {
  lapply(views, function(v) {
    mat <- as.matrix(v)
    n_feat <- nrow(mat)
    frac_obs <- mean(!is.na(mat))
    weight <- sqrt(frac_obs) / sqrt(n_feat)
    mat * weight  # scaled matrix, NAs remain as-is
  })
}

views_scaled <- scale_views_soft(views)


scale_dominant_block <- function(views, block_name) {
  mat <- as.matrix(views[[block_name]])
  n_feat <- nrow(mat)
  frac_obs <- mean(!is.na(mat))
  weight <- frac_obs / sqrt(n_feat)
  views[[block_name]] <- mat * weight
  views
}

views_scaled <- scale_dominant_block(views, "genes_expression")


# Extract the genes matrix
mat <- t(views$genes_expression)

# Convert to long format
mat_long <- melt(mat, na.rm = TRUE)
colnames(mat_long) <- c("Sample", "Gene", "Expression")

# Basic boxplot per gene
ggplot(mat_long, aes(x = factor(Gene), y = Expression)) +
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Gene", y = "Expression", title = "Distribution of Expression per Gene")

remove_outliers_per_feature <- function(mat) {
  for (feature in 1:nrow(mat)) {
    vals <- mat[feature, ]
    # Skip if all NA
    if (all(is.na(vals))) next
    
    Q1 <- quantile(vals, 0.25, na.rm = TRUE)
    Q3 <- quantile(vals, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    
    lower <- Q1 - 1.5 * IQR_val
    upper <- Q3 + 1.5 * IQR_val
    
    # Replace outliers with NA
    mat[feature, vals < lower | vals > upper] <- NA
  }
  return(mat)
}

# Apply to your genes_expression view
views$genes_expression <- remove_outliers_per_feature(views$genes_expression)

# Compute per-block variance (mean variance across features)
block_vars <- sapply(views, function(mat) {
  # Variance per feature (row), ignoring NAs
  feature_vars <- apply(mat, 1, var, na.rm = TRUE)
  # Mean variance of all features in the block
  mean(feature_vars, na.rm = TRUE)
})

# Median variance across all blocks
median_var <- median(block_vars, na.rm = TRUE)

# Compute block weights
block_weights <- sqrt(median_var / block_vars)

# Scale each block
views_scaled <- mapply(function(mat, w) {
  mat * w  # scale the entire block by its weight
}, views, block_weights, SIMPLIFY = FALSE)

views_scaled$genes_expression <- t(scale(t(views$genes_expression), center = TRUE, scale = TRUE))
views_scaled$metabolites <- t(scale(t(views_scaled$metabolites), center = TRUE, scale = TRUE))
views_scaled$spectra <- t(scale(t(views_scaled$spectra), center = TRUE, scale = TRUE))
views_scaled$others <- t(scale(t(views_scaled$others), center = TRUE, scale = TRUE))

# 4. Create MOFA object
mofa_cov <- create_mofa(views)

# Step 2. Set options
data_opts  <- get_default_data_options(mofa_cov)
data_opts$scale_views <- FALSE # disable automatic scaling
#data_opts$block_weights <- c(genes = 1, metabolites = 2, spectra = 2)  # upweight smaller blocks

model_opts <- get_default_model_options(mofa_cov)
model_opts$num_factors <- 3
model_opts$spikeslab_weights$use <- FALSE
# Enable sparsity on weights
#model_opts$spikeslab_weights$use <- TRUE
#model_opts$spikeslab_weights$alpha0 <- list(genes = 1e-16, spectra = 1e-2)
#model_opts$spikeslab_weights$beta0 <- list(genes = 10, spectra = 1)

#model_opts$sparsity <- list(factors = 0.5, views = NULL)


train_opts <- get_default_training_options(mofa_cov)
#train_opts$maxiter <- 3000
#train_opts$
train_opts$convergence_mode <- "medium"

# (optional) Explicit likelihoods
#likelihoods <- c(
# "gene_expression" = "gaussian",
#"metabolites"     = "gaussian",
#"spectra"         = "gaussian"
#)

# Step 3. Prepare the model
mofa_cov<- prepare_mofa(
  object            = mofa_cov,
  data_options      = data_opts,
  model_options     = model_opts,
  training_options  = train_opts
  #likelihoods       = likelihoods
)


# Step 4. Train the model
mofa_cov <- run_mofa(mofa_cov,use_basilisk = TRUE)

saveRDS(mofa_cov, file = "mofa_leaf_results_complete201125.rds")

mofa_cov <- readRDS("mofa_leaf_results_complete201125.rds")
input_data <- metabs_scaled
# 1. Identify which PC columns are entirely NA
na_columns <- which(apply(is.na(input_data), 2, all))

# 2. Check which PC columns have very high NA proportions (e.g., > 90%)
sparse_columns <- which(apply(is.na(input_data), 2, sum) / nrow(input_data) > 0.9)

# 3. Create a list of columns to remove
columns_to_remove <- unique(c(na_columns, sparse_columns))

# 4. Filter the data matrix (and then re-create your MOFA object)
clean_data <- input_data[, -columns_to_remove]

# 2. Extract and Clean the ELBO vector
elbo_vector <- mofa_cov@training_stats$elbo

# CRITICAL STEP: Identify and remove NaN/Inf values
# We keep only values that are not NaN and are finite.
valid_indices <- !is.nan(elbo_vector) & is.finite(elbo_vector)
clean_elbo_vector <- elbo_vector[valid_indices]

# 3. Create a clean data frame for plotting
convergence_df <- data.frame(
  iteration = (1:length(elbo_vector))[valid_indices], # Use the correct iteration index
  ELBO = clean_elbo_vector
)

# 4. Generate the Plot (Single Run)
p <- ggplot(convergence_df, 
            aes(x = iteration, y = ELBO)) +
  
  # Use geom_point to show every recorded value clearly
  geom_point(color = "darkred", size = 1.5) +
  # Use geom_line to connect the points
  geom_line(color = "darkblue", linewidth = 0.8) +
  
  labs(title = paste("MOFA Convergence Check (Valid Iterations:", nrow(convergence_df), ")"), 
       x = "Iteration Number", 
       y = "ELBO (Evidence Lower Bound)") +
  
  scale_y_continuous(labels = scales::scientific) +
  theme_bw()

# 5. Display the Plot
print(p)

# 1. Access the ELBO vector inside the @training_stats slot
elbo_vector <- mofa_cov@training_stats$elbo

# 2. Extract the last element of the vector (the final converged value)
final_elbo_value <- tail(elbo_vector, 1)

# 3. Print the value
print(paste("The Final Converged ELBO Value is:", final_elbo_value))

# 1. Define parameters for the loop
N_RUNS <- 10
FACTORS_TO_TEST <- 3
MAX_ITER <- 5000
TOLERANCE_VALUE <- 0.01

# 2. Initialize storage lists

model_results <- list()
elbo_summary <- data.frame(Run_ID = character(), Final_ELBO = numeric())

# 3. Loop through the runs
for (i in 1:N_RUNS) {
  
  current_seed <- i * 100 
  current_run_id <- paste0("Run_", i)
  cat(paste("\n--- Starting", current_run_id, " (Factors:", FACTORS_TO_TEST, ", Seed:", current_seed, ") ---\n"))
  
  # --- DYNAMIC PREPARATION (Inside Loop - Minimalist) ---
  
  # 4. Get default options
  data_opts <- get_default_data_options(mofa_cov)
  data_opts$scale_views <- FALSE # Keep explicit scaling disabled
  
  model_opts <- get_default_model_options(mofa_cov)
  model_opts$num_factors <- FACTORS_TO_TEST # Set the stable factor count
  model_opts$spikeslab_weights$use <- FALSE
  
  train_opts <- get_default_training_options(mofa_cov)
  train_opts$maxiter <- MAX_ITER
  #train_opts$tolerance <- TOLERANCE_VALUE
  train_opts$convergence_mode <- "slow"
  
  # CRITICAL: We DO NOT define or modify train_opts at all.
  # We will let prepare_mofa use its internal defaults for training.
  
  # CRITICAL: Set the unique seed globally for the next run
  set.seed(current_seed) 
  
  # 5. Prepare the model 
  current_mofa_object <- prepare_mofa(
    object = mofa_cov,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
    # CRITICAL: training_options ARGUMENT IS REMOVED ENTIRELY
  )
  
  # --- TRAINING EXECUTION ---
  
  tryCatch({
    # 6. Train the model 
    current_model <- run_mofa(
      object = current_mofa_object, 
      use_basilisk = TRUE
    )
    
    # Robustly extract the final ELBO, ignoring NaN/Inf values
    elbo_vector <- current_model@training_stats$elbo
    valid_elbo_values <- elbo_vector[!is.nan(elbo_vector) & is.finite(elbo_vector)]
    
    if (length(valid_elbo_values) > 0) {
      final_elbo <- tail(valid_elbo_values, 1)
      
      # Store the successful result and the ELBO value
      model_results[[current_run_id]] <- current_model
      elbo_summary <- rbind(elbo_summary, data.frame(Run_ID = current_run_id, Final_ELBO = final_elbo))
      cat(paste("✅ Success! Final Valid ELBO:", final_elbo, "\n"))
    } else {
      cat("❌ Run failed: No valid ELBO values recorded.\n")
    }
    
  }, error = function(e) {
    cat(paste("❌ Run failed with fatal R error:", e$message, "\n"))
  })
}
### check whether there are all NAs for eevery sample of each feature
apply(t(metabs_scaled), MARGIN = 1, FUN = function(x) {
  +     all(is.na(x))
  + })

# 4. Define parameters for the loop
N_RUNS <- 10
MAX_ITER <- 5000
TOLERANCE_VALUE <- 0.01
FACTORS_TO_TEST <- 5
# 5. Initialize storage lists
# CRITICAL STABILIZATION PARAMETER
# This limits the size of the parameter updates (like a global learning rate).
# A smaller value (e.g., 0.1) can stop large, unstable jumps that cause NaNs.
DELTA_MIN_VALUE <- 0.1
model_results <- list()
elbo_summary <- data.frame(Run_ID = character(), Final_ELBO = numeric())

# 6. Loop through the runs, changing the seed each time
for (i in 1:N_RUNS) {
  
  current_seed <- i * 100 
  cat(paste("--- Starting Run", i, " (Factors: 5, Seed:", current_seed, ") ---\n"))
  

  #--- DYNAMIC PREPARATION (Inside Loop) ---
    
    # 5. Get default options (MUST be inside the loop to reset the object)
    data_opts <- get_default_data_options(mofa_cov)
  data_opts$scale_views <- FALSE # Keep explicit scaling disabled
  
  model_opts <- get_default_model_options(mofa_cov)
  model_opts$num_factors <- FACTORS_TO_TEST # Set the stable factor count
  model_opts$spikeslab_weights$use <- FALSE
  
  train_opts <- get_default_training_options(mofa_cov)
  train_opts$maxiter <- MAX_ITER
  #train_opts$tolerance <- TOLERANCE_VALUE
  train_opts$convergence_mode <- "slow"
  #train_opts$delta_min <- DELTA_MIN_VALUE
  
  # CRITICAL FIX: Set the unique seed dynamically in the training options
  #train_opts$seed <- current_seed 
  
  # CRITICAL FIX: Set the seed globally before the run_mofa call
  set.seed(current_seed)
  # 6. Prepare the model (MUST be inside the loop to apply the new seed)
  current_mofa_object <- prepare_mofa(
    object = mofa_cov,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  # --- TRAINING EXECUTION ---
  
  tryCatch({
    # 7. Train the model (No extra arguments needed here, they are in the object)
    current_model <- run_mofa(
      object = current_mofa_object, 
      use_basilisk = TRUE
    )
    
    # Robustly extract the final ELBO, ignoring NaN/Inf values
    elbo_vector <- current_model@training_stats$elbo
    valid_elbo_values <- elbo_vector[!is.nan(elbo_vector) & is.finite(elbo_vector)]
    
    if (length(valid_elbo_values) > 0) {
      final_elbo <- tail(valid_elbo_values, 1)
      
      # Store the successful result and the ELBO value
      model_results[[current_run_id]] <- current_model
      elbo_summary <- rbind(elbo_summary, data.frame(Run_ID = current_run_id, Final_ELBO = final_elbo))
      cat(paste("✅ Success! Final Valid ELBO:", final_elbo, "\n"))
    } else {
      cat("❌ Run failed: No valid ELBO values recorded.\n")
    }
    
  }, error = function(e) {
    cat(paste("❌ Run failed with fatal R error:", e$message, "\n"))
  })
}

# --- 8. Identify the Best Model ---
cat("\n------------------------------------------------------\n")
if (nrow(elbo_summary) > 0) {
  # Find the row corresponding to the highest ELBO (least negative number)
  best_run_data <- elbo_summary[which.max(elbo_summary$Final_ELBO), ]
  
  # This is your final, most robust model!
  best_mofa_model <- model_results[[best_run_data$Run_ID]] 
  
  cat(paste("🏆 Best Model ID (Highest ELBO):", best_run_data$Run_ID, "\n"))
  cat(paste("Highest Valid ELBO:", best_run_data$Final_ELBO, "\n"))
  cat("\nUse the factor weights from the 'best_mofa_model' object for final interpretation.\n")
} else {
  cat("⚠️ ERROR: All runs failed. Please check your data cleaning steps.\n")
}
cat("------------------------------------------------------\n")

# 1. Retrieve the raw ELBO vector from the best model
elbo_vector <- best_mofa_model@training_stats$elbo

# 2. Check for NaN values
num_nans <- sum(is.nan(elbo_vector))

cat("\n--- Stability Analysis of Best Model (Run 1) ---\n")
cat(paste("Total Iterations Attempted (Length of Vector):", length(elbo_vector), "\n"))
cat(paste("Number of NaN (Crashed) Iterations:", num_nans, "\n"))
cat(paste("Final Converged ELBO Value:", tail(elbo_vector[!is.nan(elbo_vector)], 1), "\n"))
cat("--------------------------------------------------\n")

# 3. Optional: Plot the first 100 values to check the start of convergence
if (length(elbo_vector) > 100) {
  # Print the first 100 ELBO values (or the full vector if small)
  cat("\nFirst 100 ELBO values (for visual inspection of the trace):\n")
  print(elbo_vector)
}
var_explained_metrics_k5 <- calculate_variance_explained(best_mofa_model)
r2_values_per_view_k5 <- var_explained_metrics_k5$r2_total$group1
total_var_explained_k5 <- mean(r2_values_per_view_k5) * 100
cat(paste("Total Variance Explained for K=5:", round(total_var_explained_k5, 2), "%\n"))

r2_values_per_view_k5 <- var_explained_metrics_k5$r2_total[[1]] 

# 3. Calculate the mean R2 across all views and convert to percentage
total_var_explained_k5 <- mean(r2_values_per_view_k5)

cat(paste("Total Variance Explained for K=5:", round(total_var_explained_k5, 2), "%\n"))

plot_variance_explained(best_mofa_model)
#########################################
# ----------------------------------------------------------------------
## 1. Define K-Optimization Parameters
# ----------------------------------------------------------------------

# Define the range of factors to test
K_TESTS <- c(3, 5, 7, 10) 

# Inner Loop Parameters
N_STABILITY_RUNS <- 10  # Number of stability checks per K
MAX_ITER <- 5000
TOLERANCE_VALUE <- 0.01
DELTA_MIN_VALUE <- 0.1 # CRITICAL STABILIZATION PARAMETER

# Storage for results across all K values
k_optimization_results <- data.frame(
  K_Value = integer(), 
  Best_Final_ELBO = numeric(), 
  Total_Variance_Explained = numeric()
)

# ----------------------------------------------------------------------
## 2. OUTER LOOP: Iterate through different K values
# ----------------------------------------------------------------------
for (K_VALUE in K_TESTS) {
  
  cat(paste("\n======================================================\n"))
  cat(paste("STARTING OPTIMIZATION FOR K =", K_VALUE, "\n"))
  cat(paste("======================================================\n"))
  
  # Initialize storage for the stability runs within this K
  inner_model_results <- list()
  inner_elbo_summary <- data.frame(Run_ID = character(), Final_ELBO = numeric())
  
  
  # --- INNER LOOP: Stability Test for current K ---
  for (i in 1:N_STABILITY_RUNS) {
    
    current_seed <- i * 100
    current_run_id <- paste0("K", K_VALUE, "_Run_", i)
    cat(paste("\n--- Starting", current_run_id, " (Seed:", current_seed, ") ---\n"))
    
    # --- Prepare MOFA Object ---
    data_opts <- get_default_data_options(mofa_cov)
    data_opts$scale_views <- FALSE
    
    model_opts <- get_default_model_options(mofa_cov)
    model_opts$num_factors <- K_VALUE # <--- DYNAMIC K VALUE
    model_opts$spikeslab_weights <- list(use = FALSE)
    #model_opts$spikeslab_weights$use <- FALSE
    #tolerance = TOLERANCE_VALUE 
    #delta_min = DELTA_MIN_VALUE
    #train_opts <- list(
     # maxiter = MAX_ITER,
     # convergence_mode = "slow"
    #)

    #train_opts <- get_default_training_options(mofa_cov)
    #train_opts$maxiter <- MAX_ITER
    #train_opts$tolerance <- TOLERANCE_VALUE
    #train_opts$convergence_mode <- "slow"
    #train_opts$delta_min <- DELTA_MIN_VALUE
    
    set.seed(current_seed)
    
    current_mofa_object <- prepare_mofa(
      object = mofa_cov, data_options = data_opts,
      model_options = model_opts#,training_options = train_opts
    )
    
    # --- Run MOFA and Store Best Result ---
    tryCatch({
      current_model <- run_mofa(object = current_mofa_object, use_basilisk = TRUE)
      
      elbo_vector <- current_model@training_stats$elbo
      valid_elbo_values <- elbo_vector[!is.nan(elbo_vector) & is.finite(elbo_vector)]
      
      if (length(valid_elbo_values) > 0) {
        final_elbo <- tail(valid_elbo_values, 1)
        inner_model_results[[current_run_id]] <- current_model
        inner_elbo_summary <- rbind(inner_elbo_summary, data.frame(Run_ID = current_run_id, Final_ELBO = final_elbo))
        cat(paste("✅ Success! Final Valid ELBO:", final_elbo, "\n"))
      } else {
        cat("❌ Run failed: No valid ELBO values recorded.\n")
      }
      
    }, error = function(e) {
      cat(paste("❌ Run failed with fatal R error:", e$message, "\n"))
    })
  }
  
  # --- Summarize Results for the current K ---
  
  # --- Summarize Results for the current K ---
  if (nrow(inner_elbo_summary) > 0) {
    # Find the best stable model for this K
    best_run_data <- inner_elbo_summary[which.max(inner_elbo_summary$Final_ELBO), ]
    best_mofa_model_for_K <- inner_model_results[[best_run_data$Run_ID]]
    
    # Calculate Total Variance Explained (Using the working logic)
    var_explained_metrics <- calculate_variance_explained(best_mofa_model_for_K)
    
    # Extract the numeric R2 vector using numeric indexing
    r2_values_per_view <- var_explained_metrics$r2_total[[1]] 
    
    # Calculate the mean R2 (Total Variance Explained percentage)
    total_variance_explained <- mean(r2_values_per_view) 
    
    # FIX 1: Print the result immediately for confirmation
    cat(paste("Total R2 for K=", K_VALUE, " across all views: ", 
              round(total_variance_explained, 2), "%\n", sep=""))
    
    # Store the results for the outer analysis
    k_optimization_results <- rbind(k_optimization_results, data.frame(
      K_Value = K_VALUE,
      Best_Final_ELBO = best_run_data$Final_ELBO,
      Total_Variance_Explained = total_variance_explained
    ))
    
    cat(paste("\nSummary for K =", K_VALUE, ": Total Variance Explained (R2):", round(total_variance_explained, 2), "%\n"))
    
  } else {
    cat(paste("\nWARNING: All runs for K =", K_VALUE, "failed. Skipping to next K.\n"))
    # If all runs fail, we need to ensure this K value is NOT added to the plot data
  }
}
  
  # --- END OF SUMMARIZE RESULTS BLOCK ---
# ----------------------------------------------------------------------
## 3. Final K-Selection Analysis
# ----------------------------------------------------------------------

cat("\n\n======================================================\n")
cat("K-OPTIMIZATION RESULTS (Variance Explained)\n")
cat("======================================================\n")
print(k_optimization_results)

# Plot the results to find the elbow point
if (nrow(k_optimization_results) > 1) {
  plot(
    x = k_optimization_results$K_Value, 
    y = k_optimization_results$Total_Variance_Explained, 
    type = "b", 
    main = "Variance Explained vs. Number of Factors (K)",
    xlab = "Number of Factors (K)",
    ylab = "Total Variance Explained (%)",
    lwd = 2, col = "blue"
  )
  text(k_optimization_results$K_Value, k_optimization_results$Total_Variance_Explained, 
       labels = paste0(round(k_optimization_results$Total_Variance_Explained, 1), "%"), 
       pos = 4, offset = 0.5)
} else {
  cat("\nCannot plot results: Only one K-value was successfully tested or data is corrupted.\n")
}


#1. Filter the results table for K=7
k7_summary <- k_optimization_results[k_optimization_results$K_Value == 7, ]

# 2. Identify the Run_ID with the highest ELBO for K=7
best_run_id_k7 <- k7_summary$Best_Final_ELBO[which.max(k7_summary$Best_Final_ELBO)]

# 3. CRITICAL: Find the full run ID string (e.g., "K7_Run_3")
# You'll need to look back at the ELBO summary you printed during the run 
# to match the ELBO value to the Run_ID.
# The exact code to extract the ID is complex due to the nested structure, 
# so manual matching is often easiest if you have the console output.
# Example: Find the Run_ID matching best_run_id_k7 in your console output.

# This only works if you ran the K=7 block last and 'inner_model_results' 
# still contains the K=7 models.
best_mofa_trained_k7 <- inner_model_results[["K7_Run_3"]] 

# Now plot the correlation
plot_factor_cor(best_mofa_trained_k7)


K_VALUE <- 3
current_seed <- 300 # Use a consistent seed for reproducibility

data_opts <- get_default_data_options(mofa_cov)
data_opts$scale_views <- FALSE # Use the successful scaling fix
data_opts$use_float32 <- FALSE

model_opts <- get_default_model_options(mofa_cov)
model_opts$num_factors <- K_VALUE
model_opts$spikeslab_weights <- list(use = FALSE)

train_opts <- get_default_training_options(mofa_cov)
train_opts$maxiter <- 5000
#train_opts$tolerance <- 0.0000001
train_opts$convergence_mode <- "slow"  # Forces more iterations
train_opts$seed <- 300                # Change this number to escape the current minimum

#set.seed(current_seed)

# --- 2. Prepare and Train the model outside the loop ---
mofa_prepared_k3 <- prepare_mofa(
  object = mofa_cov, 
  data_options = data_opts, 
  model_options = model_opts,
  training_options = train_opts
  # Use default training options (no train_opts argument)
)

mofa_trained_k3 <- run_mofa(object = mofa_prepared_k3, use_basilisk = TRUE)
plot_variance_explained(mofa_trained_k3) + 
  theme_minimal(base_size = 16)

plot_factor_cor(mofa_trained_k3)
##############################
K_VALUE <- 5
current_seed <- 300 # Use a consistent seed for reproducibility

data_opts <- get_default_data_options(mofa_cov)
data_opts$scale_views <- FALSE # Use the successful scaling fix

model_opts <- get_default_model_options(mofa_cov)
model_opts$num_factors <- K_VALUE
model_opts$spikeslab_weights <- list(use = FALSE)

train_opts <- get_default_training_options(mofa_cov)
#train_opts$maxiter <- 5000
#train_opts$tolerance <- TOLERANCE_VALUE
train_opts$convergence_mode <- "slow"

set.seed(current_seed)

# --- 2. Prepare and Train the model outside the loop ---
mofa_prepared_k5 <- prepare_mofa(
  object = mofa_cov, 
  data_options = data_opts, 
  model_options = model_opts
  #training_options = train_opts
  # Use default training options (no train_opts argument)
)

mofa_trained_k5 <- run_mofa(object = mofa_prepared_k5, use_basilisk = TRUE)

# --- 3. Plot the Factor Correlation (This will now work) ---
cat("\n--- Checking Factor Correlation for K=7 ---\n")
plot_factor_cor(mofa_trained_k5)

plot_variance_explained(mofa_trained_k5)
#############

# --- 1. Define Options for K=7 ---
K_VALUE <- 7
current_seed <- 300 # Use a consistent seed for reproducibility

data_opts <- get_default_data_options(mofa_cov)
data_opts$scale_views <- FALSE # Use the successful scaling fix

model_opts <- get_default_model_options(mofa_cov)
model_opts$num_factors <- K_VALUE
model_opts$spikeslab_weights <- list(use = FALSE)

train_opts <- get_default_training_options(mofa_cov)
#train_opts$maxiter <- 5000
#train_opts$tolerance <- TOLERANCE_VALUE
train_opts$convergence_mode <- "slow"

set.seed(current_seed)

# --- 2. Prepare and Train the model outside the loop ---
mofa_prepared_k7 <- prepare_mofa(
  object = mofa_cov, 
  data_options = data_opts, 
  model_options = model_opts
  #training_options = train_opts
  # Use default training options (no train_opts argument)
)

mofa_trained_k7 <- run_mofa(object = mofa_prepared_k7, use_basilisk = TRUE)

# --- 3. Plot the Factor Correlation (This will now work) ---
cat("\n--- Checking Factor Correlation for K=7 ---\n")
plot_factor_cor(mofa_trained_k7)

plot_variance_explained(mofa_trained_k7)
###########################################
factor_scores_object <- get_factors(mofa_trained_k3)

# 2. Convert to a standard matrix or data frame if it's a list (common MOFA behavior)
if (is.list(factor_scores_object)) {
  # If the model has groups, get_factors() returns a list.
  # We assume you have a single group, so we take the first element (the matrix).
  # If you have multiple groups, use the first group's scores.
  factor_scores_matrix <- factor_scores_object[[1]] 
} else {
  # If it's already a single matrix-like object, convert it explicitly
  factor_scores_matrix <- as.matrix(factor_scores_object)
}

# 3. Calculate the correlation matrix (R)
correlation_matrix <- cor(factor_scores_matrix)
print(correlation_matrix)


# 1. Retrieve the Factor Scores matrix (Z) again, ensuring it's a matrix
factor_scores_object <- get_factors(best_mofa_model)
if (is.list(factor_scores_object)) {
  factor_scores_matrix <- factor_scores_object[[1]] 
} else {
  factor_scores_matrix <- as.matrix(factor_scores_object)
}

cat("Factor Scores matrix ready for clustering.\n")
# 2. Calculate the distance matrix between samples
# We use Euclidean distance across the 3 Factor dimensions.
distance_matrix <- dist(factor_scores_matrix, method = "euclidean")

# 3. Perform Hierarchical Clustering
# 'ward.D2' method minimizes the total within-cluster variance.
hclust_results <- hclust(distance_matrix, method = "ward.D2")

cat("\nHierarchical Clustering performed.\n")
# 4. Cut the dendrogram to assign samples to K=3 clusters
K_CLUSTERS <- 3
sample_clusters <- cutree(hclust_results, k = K_CLUSTERS)

cat(paste("\nSamples assigned to K =", K_CLUSTERS, "clusters.\n"))
print(table(sample_clusters))

# 5. Combine scores with cluster assignment
cluster_df <- data.frame(
  Factor1 = factor_scores_matrix[, "Factor1"],
  Factor2 = factor_scores_matrix[, "Factor2"],
  Factor3 = factor_scores_matrix[, "Factor3"],
  Cluster = as.factor(sample_clusters)
)

# 6. Calculate the Mean Factor Score for each Cluster (Cluster Centroids)
# This defines the "personality" of each cluster
cluster_centroids <- aggregate(. ~ Cluster, data = cluster_df, FUN = mean)

cat("\n## 🎯 Cluster Characterization (Mean Factor Scores)\n")
print(cluster_centroids)


# --- Define the Plotting Function ---
# This function creates a standardized plot for any two factors
create_factor_plot <- function(df, x_factor, y_factor, x_label, y_label, colors, centroids) {
  
  # Selects only the necessary factor columns and Cluster column from the centroids
  current_centroids <- centroids[, c(x_factor, y_factor, "Cluster")]
  
  p <- ggplot(df, aes_string(x = x_factor, y = y_factor, color = "Cluster")) +
    
    # Add points
    geom_point(size = 3.5, alpha = 0.8) +
    
    # Add cluster centroids
    geom_point(data = current_centroids, aes_string(x = x_factor, y = y_factor), 
               shape = 18, size = 6, color = "black") +
    
    # Add 95% confidence ellipses
    stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.1, aes(fill = Cluster)) +
    
    # Add zero lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    
    # Apply custom colors
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    
    # Labels and Title
    labs(
      title = paste("Factor", sub("Factor", "", x_factor), "vs. Factor", sub("Factor", "", y_factor), "Scores"),
      x = x_label,
      y = y_label,
      caption = "Ellipses represent the 95% confidence region for each cluster."
    ) +
    
    # Publication-ready theme
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")
  
  return(p)
}

# ----------------------------------------------------------------------
# 1. Data Setup (Self-Contained)
# ----------------------------------------------------------------------

# 1a. Retrieve the Factor Scores matrix (Z)
factor_scores_object <- get_factors(best_mofa_model)
if (is.list(factor_scores_object)) {
  factor_scores_matrix <- factor_scores_object[[1]] 
} else {
  factor_scores_matrix <- as.matrix(factor_scores_object)
}

# 1b. Create the full plot data frame including all three factors
plot_data_df <- data.frame(
  Sample_ID = rownames(factor_scores_matrix),
  Factor1 = factor_scores_matrix[, "Factor1"], 
  Factor2 = factor_scores_matrix[, "Factor2"], # Now included
  Factor3 = factor_scores_matrix[, "Factor3"], 
  Cluster = as.factor(sample_clusters)        
)

# 1c. Define custom colors
custom_colors <- c(
  "1" = "#1b9e77", # Cluster 1: Mild/Stable (Green)
  "2" = "#d95f02", # Cluster 2: High Growth/High Stress (Orange)
  "3" = "#7570b3"  # Cluster 3: Arrested Growth/Defended (Purple)
)

# 1d. Calculate cluster centroids on ALL three factor scores
cluster_centroids <- aggregate(. ~ Cluster, 
                               data = plot_data_df[, c("Factor1", "Factor2", "Factor3", "Cluster")], 
                               FUN = mean)

# Define labels based on biological interpretation
F1_LABEL <- "Factor 1 (Biomass and Photosynthesis)"
F2_LABEL <- "Factor 2 (Protective Pigments/Defense)"
F3_LABEL <- "Factor 3 (Acute Stress and Senescence)"

# ----------------------------------------------------------------------
# 2. Generate and Print All Plots
# ----------------------------------------------------------------------

# --- Plot 1: Factor 1 vs Factor 3 (Original Plot) ---
p_f1_f3 <- create_factor_plot(
  df = plot_data_df, 
  x_factor = "Factor1", 
  y_factor = "Factor3",
  x_label = F1_LABEL, 
  y_label = F3_LABEL,
  colors = custom_colors,
  centroids = cluster_centroids
)
cat("\n--- Plot 1: Factor 1 vs Factor 3 (Growth vs Acute Stress) ---\n")
print(p_f1_f3)


# --- Plot 2: Factor 1 (Growth) vs Factor 2 (Defense) ---
p_f1_f2 <- create_factor_plot(
  df = plot_data_df, 
  x_factor = "Factor1", 
  y_factor = "Factor2",
  x_label = F1_LABEL, 
  y_label = F2_LABEL,
  colors = custom_colors,
  centroids = cluster_centroids
)

cat("\n--- Plot 2: Factor 1 vs Factor 2 (Growth vs Defense) ---\n")
print(p_f1_f2)

# --- Plot 3: Factor 2 (Defense) vs Factor 3 (Acute Stress) ---
p_f2_f3 <- create_factor_plot(
  df = plot_data_df, 
  x_factor = "Factor2", 
  y_factor = "Factor3",
  x_label = F2_LABEL, 
  y_label = F3_LABEL,
  colors = custom_colors,
  centroids = cluster_centroids
)

cat("\n--- Plot 3: Factor 2 vs Factor 3 (Defense vs Acute Stress) ---\n")
print(p_f2_f3)

###### rerun with k=4 ###########
#1. Rerun Clustering with K=4
# ----------------------------------------------------------------------

# NOTE: Assuming hclust_results object is available from the previous step
# hclust_results <- hclust(distance_matrix, method = "ward.D2")

K_CLUSTERS_NEW <- 4
sample_clusters <- cutree(hclust_results, k = K_CLUSTERS_NEW)

cat(paste("\nNew Sample assignment to K =", K_CLUSTERS_NEW, "clusters.\n"))
print(table(sample_clusters))

# ----------------------------------------------------------------------
# 2. Data Preparation for K=4 Plots
# ----------------------------------------------------------------------

# 2a. Retrieve the Factor Scores matrix (Z)
factor_scores_object <- get_factors(best_mofa_model)
if (is.list(factor_scores_object)) {
  factor_scores_matrix <- factor_scores_object[[1]] 
} else {
  factor_scores_matrix <- as.matrix(factor_scores_object)
}

# 2b. Create the full plot data frame including all three factors
plot_data_df <- data.frame(
  Sample_ID = rownames(factor_scores_matrix),
  Factor1 = factor_scores_matrix[, "Factor1"], 
  Factor2 = factor_scores_matrix[, "Factor2"],
  Factor3 = factor_scores_matrix[, "Factor3"], 
  Cluster = as.factor(sample_clusters)        
)

# 2c. Define custom colors for K=4 (using a palette for better separation)
# We will use four colors from a colorblind-friendly palette
custom_colors <- c(
  "1" = "#1B9E77", # New Cluster 1 (e.g., Very Low Stress)
  "2" = "#D95F02", # New Cluster 2 (e.g., Moderate Stress)
  "3" = "#7570B3", # New Cluster 3 (e.g., High Stress)
  "4" = "#E7298A"  # New Cluster 4 (e.g., Arrested Defense - Purple/Pink)
)

# 2d. Calculate cluster centroids on ALL three factor scores
cluster_centroids <- aggregate(. ~ Cluster, 
                               data = plot_data_df[, c("Factor1", "Factor2", "Factor3", "Cluster")], 
                               FUN = mean)

cat("\n## 🎯 New K=4 Cluster Characterization (Mean Factor Scores)\n")
print(cluster_centroids)


# ----------------------------------------------------------------------
# 3. Plotting Function (Reusable from previous script)
# ----------------------------------------------------------------------

create_factor_plot <- function(df, x_factor, y_factor, x_label, y_label, colors, centroids) {
  
  current_centroids <- centroids[, c(x_factor, y_factor, "Cluster")]
  
  p <- ggplot(df, aes_string(x = x_factor, y = y_factor, color = "Cluster")) +
    geom_point(size = 3.5, alpha = 0.8) +
    geom_point(data = current_centroids, aes_string(x = x_factor, y = y_factor), 
               shape = 18, size = 6, color = "black") +
    stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.1, aes(fill = Cluster)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste("Factor", sub("Factor", "", x_factor), "vs. Factor", sub("Factor", "", y_factor), "Scores (K=4)"),
      x = x_label,
      y = y_label,
      caption = "Ellipses represent the 95% confidence region for each cluster."
    ) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")
  
  return(p)
}


# Define labels based on biological interpretation
F1_LABEL <- "Factor 1 (Biomass and Photosynthesis)"
F2_LABEL <- "Factor 2 (Protective Pigments/Defense)"
F3_LABEL <- "Factor 3 (Acute Stress and Senescence)"

# ----------------------------------------------------------------------
# 4. Generate and Print All K=4 Plots
# ----------------------------------------------------------------------

# --- Plot 1: Factor 1 vs Factor 3 (Growth vs Acute Stress) ---
p_f1_f3_k4 <- create_factor_plot(
  df = plot_data_df, x_factor = "Factor1", y_factor = "Factor3",
  x_label = F1_LABEL, y_label = F3_LABEL, colors = custom_colors, centroids = cluster_centroids
)
cat("\n--- Plot 1: Factor 1 vs Factor 3 (K=4) ---\n")
print(p_f1_f3_k4)


# --- Plot 2: Factor 1 (Growth) vs Factor 2 (Defense) ---
p_f1_f2_k4 <- create_factor_plot(
  df = plot_data_df, x_factor = "Factor1", y_factor = "Factor2",
  x_label = F1_LABEL, y_label = F2_LABEL, colors = custom_colors, centroids = cluster_centroids
)

cat("\n--- Plot 2: Factor 1 vs Factor 2 (K=4) ---\n")
print(p_f1_f2_k4)

# --- Plot 3: Factor 2 (Defense) vs Factor 3 (Acute Stress) ---
p_f2_f3_k4 <- create_factor_plot(
  df = plot_data_df, x_factor = "Factor2", y_factor = "Factor3",
  x_label = F2_LABEL, y_label = F3_LABEL, colors = custom_colors, centroids = cluster_centroids
)

cat("\n--- Plot 3: Factor 2 vs Factor 3 (K=4) ---\n")
print(p_f2_f3_k4)



################################

# Example: factor 1 correlations with total expression
cor(total_expr, Z[,1], use = "complete.obs")
cor(total_expr_metabs, Z[,1], use = "complete.obs")

# Correlation with all factors
apply(Z, 2, function(f) cor(total_expr, f, use = "complete.obs"))
apply(Z, 2, function(f) cor(total_expr_metabs, f, use = "complete.obs"))

###### CORRELATION WITH covarIATES ########
# Combine them for plotting
library(tidyr)

plot_data <- as_tibble(cbind(Z, covariate = result_fruit_7daybef$Water_P.y))

# Convert to long format
plot_long <- plot_data %>%
  pivot_longer(
    cols = starts_with("Factor"),  # all factor columns
    names_to = "Factor",
    values_to = "Score"
  )
#Compute correlations per factor
cor_df <- plot_long %>%
  group_by(Factor) %>%
  summarise(
    r = cor(Score, covariate, use = "complete.obs",method="spearman"),
    ymax = max(Score, na.rm = TRUE),  # y-position for label
    xmin = min(covariate, na.rm = TRUE), # x-position for label
    .groups = "drop"
  ) %>%
  mutate(label = paste0("spearman's r = ", round(r, 2)))

#plot with annotaion of correlations
ggplot(plot_long, aes(x = covariate, y = Score)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_text(
    data = cor_df,
    aes(x = xmin, y = ymax, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 3
  ) +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(x = "Mean BWP cumulative ", y = "Factor Score", title = "MOFA Factors vs Base Water Potential \n (using gene expression levels and excluding PAL1)\n measured in leafs 7 days before harvest") +
  theme_minimal()
### not linear curvas fitting
ggplot(plot_long, aes(x = covariate, y = Score)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "red") +
  geom_smooth(method = "loess", span = 0.75, se = TRUE, color = "blue") +
  geom_text(
    data = cor_df,
    aes(x = xmin, y = ymax, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = 3
  ) +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(
    x = "Mean BWP cumulative",
    y = "Factor Score",
    title = "MOFA Factors vs Base Water Potential \n (using gene expression levels and excluding PAL1)\n measured in leafs 7 days before harvest)") +
  theme_minimal()
############################################

#####FOR CATEGORICAL######
# Convert matrix to data frame
# Convert factor_scores matrix to data frame
factor_scores_df <- as.data.frame(Z)

plot_data <- cbind(factor_scores_df, covariate = result_fruit$Test_site)
#plot_data$covariate_num <- as.numeric(plot_data$covariate)

# --- ANOVA + correlation results ---
anova_results <- data.frame(Factor = character(),
                            ANOVA_p = numeric(),
                            #Pearson_r = numeric(),
                            #Pearson_p = numeric(),
                            stringsAsFactors = FALSE)

for (factor in colnames(factor_scores_df)) {
  # ANOVA
  anova_model <- aov(plot_data[[factor]] ~ covariate, data = plot_data)
  anova_p <- summary(anova_model)[[1]]["Pr(>F)"][1]
  
  # Pearson correlation
 # cor_test <- cor.test(plot_data[[factor]], 
   #                    plot_data$covariate_num, 
      #                 method = "pearson")
  
  # Store results
  anova_results <- rbind(anova_results,
                         data.frame(Factor = factor,
                                    ANOVA_p = anova_p
                                    #Pearson_r = cor_test$estimate,
                                    #Pearson_p = cor_test$p.value)
                                    )
  )
}
print(anova_results)


# --- Pearson correlation results ---
cor_results <- data.frame(Factor = character(),
                          Pearson_r = numeric(),
                          Pearson_p = numeric(),
                          stringsAsFactors = FALSE)

for (factor in colnames(factor_scores_df)) {
  cor_test <- cor.test(plot_data[[factor]], 
                       plot_data$covariate_num, 
                       method = "pearson")
  
  cor_results <- rbind(cor_results,
                       data.frame(Factor = factor,
                                  Pearson_r = cor_test$estimate,
                                  Pearson_p = cor_test$p.value))
}

print(cor_results)
# Convert factor_scores to long format for facet plotting
plot_long <- plot_data %>%
  mutate(sample = rownames(plot_data)) %>%  # optional ID
  pivot_longer(cols = colnames(factor_scores_df), 
               names_to = "Factor", values_to = "Score")

plot_long <- plot_data %>%
  pivot_longer(cols = colnames(factor_scores_df), 
               names_to = "Factor", values_to = "Score")

# We'll put the text at the top of each facet
annotation_df <- plot_long %>%
  group_by(Factor) %>%
  summarize(x = max(covariate_num), 
            y = max(Score), .groups = "drop") %>%
  left_join(cor_results, by = "Factor") %>%
  mutate(label = paste0("r = ", round(Pearson_r, 2), 
                        ", p = ", signif(Pearson_p, 2)))

# --- Plot with regression line and Pearson annotation ---
ggplot(plot_long, aes(x = covariate_num, y = Score)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ Factor, scales = "free_y") +
  geom_text(data = annotation_df, 
            aes(x = x, y = y, label = label), 
            hjust = 1, vjust = 1, size = 3.5, color = "red") +
  labs(x = "Sampling timepoint", y = "Factor Score",
       title = "Factor Scores vs Sampling Timepoint  \n (using gene expression levels and excluding PAL1)") +
  theme_minimal()

ggplot(plot_long, aes(x = covariate, y = Score)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(x = "Test site", y = "Factor Score", title = "Factor Scores vs test site") +
  theme_minimal()

ggplot(plot_long, aes(x = covariate, y = Score)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~ Factor, scales = "free_y") +
  labs(x = "Maximum Temperature (ºC) (averaged across 7 days prior) ", y = "Factor Score") +
  theme_minimal()

ggplot(plot_data, aes(x = mean_Tmax, y = Factor1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "mean_Tmax (ºC)", y = "Factor1 Score (AU)", title = "MOFA Factor1 vs mean Temperature (average of 7 days prior)\n (using gene expression levels and excluding PAL1) ") +
  theme_minimal()

cor.test(plot_data$Factor1, plot_data$mean_Tmax)


#### Diagnostics ###

ve <- get_variance_explained(MOFA_pcs)   # MOFA2 object name you used
# This prints a table you can inspect: which factors have most variance in which view
print(ve$r2_per_factor)                  # or whatever structure get_variance_explained returns

# 2) correlation of Factor1 with gene-level total/mean
factors_mat <- get_factors(mofa_cov, factors = "all")[[1]]

# now grab Factor 1 as numeric
f1 <- factors_mat[,1]

# mean gene signal
mg <- rowMeans(genes_scaled, na.rm=TRUE)
num_detect <- rowSums(!is.na(genes_scaled))

# correlation
cat("cor with mean gene:", cor(f1, mg, use="pairwise.complete.obs"), "\n")
cat("cor with #detected:",  cor(f1, num_detect, use="pairwise.complete.obs"), "\n")

# 3) how many weights per view are active on Factor1
W <- get_weights(MOFA_pcs)   # list with views
count_active <- sapply(W, function(mat){
  ws <- as.matrix(mat)
  sum(abs(ws[,1]) > 1e-6, na.rm=TRUE)   # threshold small
})
print(count_active)

# 4) average absolute weight magnitude per view per factor
avg_abs <- sapply(W, function(mat) colMeans(abs(as.matrix(mat)), na.rm=TRUE))
print(avg_abs)  # each column = factor; rows = views

# 5) PCA comparison: does PC1 of genes align to Factor1?
# Remove rows (genes) with any missing/infinite values
genes_clean <- genes_scaled[apply(is.finite(genes_scaled), 1, all), ]
# Remove constant columns
genes_nonconstant <- genes_clean[, apply(genes_clean, 2, function(x) sd(x, na.rm = TRUE) != 0)]
pca_g <- prcomp(genes_nonconstant, center=TRUE, scale.=TRUE)
cat("cor(F1, genes PC1):", cor(f1, pca_g$x[,1], use="pairwise.complete.obs"), "\n")

#spectra_norm <- spectra / rowSums(spectra, na.rm = TRUE)
#spectra_norm[is.na(spectra_norm)] <- 0  # handle samples with total=0

#Compare variance explained
# Original (already computed)
plot_variance_explained(MOFAobject)

# Scaled
plot_variance_explained(MOFA_scaled)

# PCA-reduced
plot_variance_explained(MOFA_pcs)

#with covariates
plot_variance_explained(mofa_cov)


# Plot top contributing features
plot_weights(
  mofa_cov,
  view = "gene_expression",   # or "metabolites", "spectra"
  factor = 1, 
  nfeatures = 20, 
  scale = TRUE                # rescales weights for comparability
)

# Preprocess depending on input type
# If RNA-seq counts
#genes_proc <- preprocess_genes(genes_mat, type = "counts")

# If qPCR Cq values
#genes_proc <- preprocess_genes(genes, type = "cq")

# Preprocess metabolites & spectra
#metab_proc <- preprocess_metabolites(metabolites)
#spec_proc  <- preprocess_spectra(spectra_norm)

#### Build the list of views/blocks #####

mofa_data_list <- list(
  "gene_expression" = t(genes_proc),
  "metabolites" = t(metab_proc),
  "spectra" = t(spec_proc)
)

MOFAobject <- create_mofa(mofa_data_list)

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE   # <--- herev
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 25  # adjust as needed

likelihoods <- c(
  "gene_expression" = "gaussian",
  "metabolites" = "gaussian",
  "spectra" = "gaussian"
)

train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter <- 1000
train_opts$convergence_mode <- "fast"

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  model_options = model_opts,
  training_options = train_opts
)


MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

#### check which block is causing the correlation that is signal and not biological
#Compute total signal per block
# Row sums per block (samples × total intensity)
total_genes <- rowSums(genes, na.rm = TRUE)
total_metabolites <- rowSums(metabolites, na.rm = TRUE)
total_spectra <- rowSums(spectra, na.rm = TRUE)


#Correlate each factor with totals
factors <- get_factors(MOFA_pcs, factors = "all", as.data.frame = FALSE)
factors <- factors[[1]]   # now samples × factors matrix
# Example: check Factor 1
#The block with the highest correlation is the one driving the warning.
#If correlation > 0.8–0.9, Factor 1 mostly reflects that block’s total magnitude.
cor(factors[,1], total_genes)        # Factor 1 vs gene totals
cor(factors[,1], total_metabolites)  # Factor 1 vs metabolite totals
cor(factors[,1], total_spectra)      # Factor 1 vs spectra totals

factors_mat <- get_factors(MOFAobject, as.data.frame = FALSE)[[1]]

# name rows/columns if not already
rownames(factors_mat) <- rownames(metadata)
colnames(factors_mat) <- paste0("Factor", 1:ncol(factors_mat))

#  Heatmap of factors
pheatmap(factors_mat,
         scale="row",      # center each sample
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         main="MOFA factor values per sample")

orig_total <- rowSums(spectra, na.rm = TRUE)

apply(factors_mat, 2, function(f) cor(f, orig_total, use="complete.obs"))

#  Correlation of factors with covariates
covariates <- metadata[, c("precipitation","humidity","temperature","timepoint","year")]  # adapt names
cor_table <- sapply(1:ncol(factors_mat), function(f) {
  sapply(covariates, function(cov) cor(factors_mat[,f], cov, use="complete.obs"))
})
colnames(cor_table) <- paste0("Factor", 1:ncol(factors_mat))
cor_table <- as.data.frame(t(cor_table))  # factors x covariates
print(cor_table)

views_names(MOFAobject)
#  Top contributing features per factor per block
blocks <- c("gene_expression","metabolites","spectra")
top_n <- 10

for(block in blocks){
  cat("\n===== Block:", block, "=====\n")
  for(f in 1:ncol(factors_mat)){
    loadings <- get_weights(MOFAobject, views=block, factors=f)[[1]]
    loadings_sorted <- sort(abs(loadings[,1]), decreasing=TRUE)  # by absolute contribution
    top_features <- names(loadings_sorted)[1:top_n]
    cat(paste0("Factor ", f, " top features: "), paste(top_features, collapse=", "), "\n")
  }
}

#Heatmap of Factor 1 correlations
# Factor 1, top 10 per block
# MOFA weights
weights_genes <- get_weights(MOFA_pcs, "gene_expression", 7)[[1]]
weights_metabs <- get_weights(MOFA_pcs, "metabolites", 7)[[1]]
weights_spectra <- get_weights(MOFA_pcs, "spectra", 7)[[1]]

# Pick top N indices
top_idx_genes <- order(abs(weights_genes[,1]), decreasing=TRUE)[1:top_n]
top_idx_metabs <- order(abs(weights_metabs[,1]), decreasing=TRUE)[1:top_n]
top_idx_spectra <- order(abs(weights_spectra[,1]), decreasing=TRUE)[1:top_n]

# Subset original matrices by index
genes_data <- genes[, top_idx_genes, drop=FALSE]
metabs_data <- metabolites[, top_idx_metabs, drop=FALSE]
spectra_data <- spectra[, top_idx_spectra, drop=FALSE]

# Combine
combined_data <- cbind(genes_data, metabs_data, spectra_data)
# Remove features with zero variance
combined_data <- combined_data[, apply(combined_data, 2, function(x) sd(x, na.rm=TRUE) > 0)]

# Impute missing values with column mean
impute_na <- function(x){
  if(all(is.na(x))) return(rep(0, length(x)))  # all NA → zeros
  x[is.na(x)] <- mean(x, na.rm=TRUE)
  return(x)
}

combined_data_imputed <- apply(combined_data, 2, impute_na)

# Remove constant columns (sd=0)
combined_data_imputed <- combined_data_imputed[, apply(combined_data_imputed, 2, sd) > 0]


# Compute correlation
cor_matrix <- cor(combined_data_imputed, use="pairwise.complete.obs")

# Pre-cluster
hc <- hclust(dist(cor_matrix))
ordered <- cor_matrix[hc$order, hc$order]

# Mask lower triangle after clustering
ordered[lower.tri(ordered)] <- NA

# Heatmap
pheatmap(ordered,
         cluster_rows=FALSE, cluster_cols=FALSE,  # already clustered
         main="Correlations among top Factor 7 features (clustered)",
         na_col="white")

#check variance explained per block per factor:
plot_variance_explained(MOFAobject)

#### 1. Factor Weights / Loadings
# Extract weights (per view)
weights <- get_weights(mofa_cov, views = "spectra")

# Turn into tidy format
weights_df <- as.data.frame(weights)
weights_df$feature <- rownames(weights_df)
weights_long <- pivot_longer(weights_df, -feature, names_to = "Factor", values_to = "Weight")

weights_long <- weights_long %>%
  mutate(Factor = factor(Factor,
                         levels = unique(Factor[order(as.numeric(gsub("\\D", "", Factor)))])))
# Plot top features per factor
top_weights <- weights_long %>% 
  group_by(Factor) %>%
  slice_max(order_by = abs(Weight), n = 20)

ggplot(top_weights, aes(x = reorder(feature, Weight), y = Weight, fill = Weight > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Factor, scales = "free_y") +
  labs(x = "Feature", y = "Weight", title = "Top MOFA Factor Weights scaling by sqrt(ncols) genes \n(using gene expression, all genes, no outliers)") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue"), guide = FALSE)

ggsave("MOFA_weights_spectra_scaledsqrtncolsgenes.pdf", width=20, height=12, dpi=300)



## plot factors in rows and views in columns
# Example: already have weights_long with feature, Factor, Weight, View
# Make sure Factor and View are factors in the correct order
#### 1 Extract weights for all views
weights <- get_weights(best_mofa_model)  # NULL gets all views
weights <- get_weights(mofa_trained_k3) 

#### 2 Convert weights to long tidy format with View column
# 1 lapply returns a list of data frames
weights_list <- lapply(names(weights), function(view_name) {
  df <- as.data.frame(weights[[view_name]])
  df$feature <- rownames(df)
  df_long <- pivot_longer(df, -feature, names_to = "Factor", values_to = "Weight")
  df_long$View <- view_name
  return(df_long)
})

weights_long <- bind_rows(weights_list)

weights_long <- weights_long %>%
  mutate(
    Factor = factor(as.character(Factor),
                    levels = unique(Factor[order(as.numeric(gsub("\\D", "", as.character(Factor))))])),
    View = factor(View, levels = c("genes_expression","metabolites","spectra","others"))
  )

#  Select top features per Factor × View
top_weights <- weights_long %>%
  group_by(Factor, View) %>%
  slice_max(order_by = abs(Weight), n = 20)

#  Plot: rows = Factor, columns = View
ggplot(top_weights, aes(x = reorder(feature, Weight), y = Weight, fill = Weight > 0)) +
  geom_col() +
  coord_flip() +
  facet_grid(rows = vars(Factor), cols = vars(View), scales = "free_y") +
  labs(x = "Feature", y = "Weight",
       title = "Top MOFA Factor Weights by Factor (rows) and View (columns)") +
  theme_minimal() +
  scale_fill_manual(values = c("lightsalmon", "cornflowerblue"), guide = FALSE) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),     # title
    axis.title = element_text(size = 14, face = "bold"),     # axis titles
    axis.text = element_text(size = 12),                     # axis text
    strip.text = element_text(size = 14, face = "bold")      # facet labels
  )

ggsave("MOFA_weights_factors_blocks_fruit_difxscale_genesvinelag_metabs_modelA_test.pdf", final_plot, width=30, height=30, dpi=300, limitsize=FALSE)

top_weights <- top_weights %>%
  mutate(
    Factor = factor(as.character(Factor),
                    levels = unique(Factor[order(as.numeric(gsub("\\D", "", as.character(Factor))))])),
    View = factor(View, levels = c("metabolites","metabolites_lagged","spectra","others"))
  )

top_weights <- top_weights %>%
  group_by(Factor, View) %>%
  arrange(Weight) %>%  # ascending or descending depending on your preference
  mutate(feature = factor(feature, levels = feature)) %>%
  ungroup()

top_weights <- top_weights %>%
  mutate(Facet = paste(Factor, View, sep = " - "))
library(ggforce)

ggplot(top_weights, aes(x = reorder(feature, Weight), y = Weight, fill = Weight > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Facet, scales = "free_y") +
  labs(x = "Feature", y = "Weight",
       title = "Top MOFA Factor Weights by Factor × View") +
  theme_minimal() +
  scale_fill_manual(values = c("lightsalmon", "cornflowerblue"), guide = FALSE)


ggplot(top_weights, aes(x = reorder(feature, Weight), y = Weight, fill = Weight > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~Facet, scales = "free_y") +
  labs(x = "Feature", y = "Weight",
       title = "Top MOFA Factor Weights per Factor × View") +
  theme_minimal() +
  scale_fill_manual(values = c("lightsalmon", "cornflowerblue"), guide = FALSE)

library(patchwork)

# Ensure Factor is ordered numerically and View is consistent
top_weights <- top_weights %>%
  mutate(
    Factor = factor(as.character(Factor),
                    levels = paste0("Factor", 1:3)),  # adjust if your factor names are different
    View = factor(View, levels = c("genes_expression", "metabolites", "spectra", "others"))
  )



top_weights <- top_weights %>%
  mutate(
    # Create a character column: "Positive" if Weight > 0, "Negative" otherwise
    color_group = ifelse(Weight > 0, "Positive", "Negative")
  )
# List of factors in reverse order so Factor1 is at the top when stacked
factors <- levels(top_weights$Factor)

# Create plots for each factor with side label
factor_rows <- lapply(factors, function(f) {
  df <- top_weights %>% filter(Factor == f)
  
  # Create plots for each view
  view_plots <- lapply(levels(df$View), function(v) {
    df_view <- df %>% filter(View == v)
    
    # Order features by weight per panel
    df_view <- df_view %>% arrange(desc(abs(Weight))) %>%
      mutate(feature = factor(feature, levels = feature))
    
    ggplot(df_view, 
           aes(x = feature, y = Weight, fill = color_group)) + # Use the new character column
      geom_col() +
      coord_flip() +
      labs(title = v) +
      theme_minimal() +
      # Map colors by the string name. This is unambiguous.
      scale_fill_manual(
        values = c("Negative" = "#dc143c", "Positive" = "cornflowerblue"), # <--- Explicitly define string-to-color mapping
        guide = FALSE
      ) + 
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10)
      )
  })
  
  # Combine the 1x4 views horizontally
  row_plot <- wrap_plots(view_plots, nrow = 1)
  
  # Create a text label for the factor on the side
  factor_label <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = f, angle = 90, size = 6, fontface = "bold")
  
  # Combine factor label and the row horizontally
  wrap_plots(factor_label, row_plot, ncol = 2, widths = c(0.05, 0.95))
})

# Stack all factor rows vertically
final_plot <- wrap_plots(factor_rows, ncol = 1) +
  plot_annotation(title = "Top MOFA Factor Weights by Factor × View")

final_plot

ggsave("MOFA_weights_factors_blocks_fruit_difxscale_genesvine_metabs_modelA_k=5_NoExtraScaling.pdf", final_plot, width=30, height=30, dpi=300, limitsize=FALSE)


#########################
#cols share y axis
# --- 1️⃣ Compute shared y-limits per View ---
ylim_table <- top_weights %>%
  group_by(View) %>%
  summarise(
    ymin = min(Weight, na.rm = TRUE),
    ymax = max(Weight, na.rm = TRUE)
  )

# --- 2️⃣ Create plots per factor, per view ---
factors <- rev(levels(top_weights$Factor))

factor_rows <- lapply(factors, function(f) {
  df <- top_weights %>% filter(Factor == f)
  
  view_plots <- lapply(levels(df$View), function(v) {
    df_view <- df %>% filter(View == v)
    yl <- ylim_table %>% filter(View == v)
    
    # Order features by abs weight
    df_view <- df_view %>%
      arrange(desc(abs(Weight))) %>%
      mutate(feature = factor(feature, levels = feature))
    
    ggplot(df_view, aes(x = feature, y = Weight, fill = Weight > 0)) +
      geom_col() +
      # Enforce common y-axis using the calculated global limits
      coord_flip(ylim = c(yl$ymin, yl$ymax)) +
      labs(title = v) +
      theme_minimal() +
      # FIX: Explicitly map TRUE/FALSE (logical) to color strings
      scale_fill_manual(
        values = c("FALSE" = "#dc143c", "TRUE" = "cornflowerblue"), 
        guide = FALSE
      ) +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10)
      )
  })
  
  # Combine the 4 views horizontally
  row_plot <- wrap_plots(view_plots, nrow = 1)
  
  # Add vertical label for the factor
  factor_label <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = f, angle = 90,
             size = 6, fontface = "bold")
  
  # Combine label + row
  wrap_plots(factor_label, row_plot, ncol = 2, widths = c(0.05, 0.95))
})

# --- 3️⃣ Stack rows vertically (Factor 7 on top, Factor 1 at bottom) ---
final_plot <- wrap_plots(rev(factor_rows), ncol = 1) +
  plot_annotation(title = "Top MOFA Factor Weights by Factor × View (Shared Y-Axis per View)")

final_plot
ggsave("MOFA_weights_factors_blocks_fruit_samexscale_genesvine_metabs_modelA_k=3.pdf", final_plot, width=30, height=30, dpi=300, limitsize=FALSE)


#### 2. Factor Scores / Values
# Extract factor scores (values per sample)
scores <- get_factors(mofa_cov, factors = "all", groups = "all")
# MOFA returns a list (per group); extract first if single group
scores_df <- as.data.frame(scores[[1]])
scores_df$sample <- rownames(scores_df)

# Long format for ggplot
scores_long <- pivot_longer(scores_df, -sample, names_to = "Factor", values_to = "Score")

# Plot sample distributions of scores
ggplot(scores_long, aes(x = Factor, y = Score)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "MOFA Factor Scores across Samples") +
  theme_minimal()

cor(scores_df[, "Factor1"], result_fruit2$Assimilable_Nitrogen__mg_L_, use = "pairwise.complete.obs")
cor(scores_df[, "Factor1"], result_fruit2$`Brix_lab__%_`, use = "pairwise.complete.obs")

scores_numeric <- scores_df[, sapply(scores_df, is.numeric)]
cor_matrix <- cor(scores_numeric, cov_block_select_fruit, use = "pairwise.complete.obs")
pheatmap(
  cor_matrix,
  
  cluster_rows = FALSE,        # cluster factors
  cluster_cols = TRUE,        # cluster genes/metabolites
  display_numbers = TRUE,     # show correlation values
  fontsize_number = 10,
  color = colorRampPalette(c("blue", "white", "red"))(50),  # blue = -1, red = 1
  breaks = seq(-1, 1, length.out = 51),  # force scale from -1 to 1
  main = "Berry: Correlation MOFA Factors vs Covariates(numeric)",
  angle_col = 45
)


colnames(spectra_fruit) <- sub("^wl_", "", colnames(spectra_fruit))
wavelengths <- as.numeric(colnames(spectra_fruit))
spectra_fruit <- as.numeric(spectra_fruit)

bin_size <- 5  # in nm or in number of features per bin
# Define bins (every 5 nm, adjust as needed)
bins <- cut(wavelengths, breaks = seq(min(wavelengths), max(wavelengths), by = bin_size), include.lowest = TRUE)
# Compute mean reflectance per bin
spectra_binned <- sapply(split(seq_along(wavelengths), bins), function(idx) {
  rowMeans(spectra_fruit[, idx, drop = FALSE], na.rm = TRUE)
})
# Convert to a proper matrix (samples × bins)
spectra_binned <- as.matrix(spectra_binned)
# Rename columns by bin midpoint (optional)
colnames(spectra_binned) <- round(tapply(wavelengths, bins, mean), 2)

cor_matrix <- cor(scores_numeric, spectra_fruit, use = "pairwise.complete.obs")
cor_plot <- pheatmap(
  cor_matrix,
  cluster_rows = FALSE,       # keep Factor1..12 order
  cluster_cols = FALSE,       # keep wavelengths in order
  display_numbers = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  angle_col = 45,
  fontsize_col = 8.5,
  fontsize_number = 4,
  main = "Correlation: MOFA Factors vs Spectral Features"
)


ggsave("Berry_covariates_mofaFactorsvsSpectra.pdf", cor_plot, width=20, height=10, dpi=300)

#################################### <3 <3 <3 ##########
#Identify top-weighted features per factor
# assume 'loadings' is a matrix (features × factors)
# and 'scores' is a matrix (samples × factors)
# and 'features_meta' contains info (e.g., gene name, metabolite, wavelength)

comb_loadings <- do.call(rbind, weights)


top_n <- 15   # number of top features per factor

top_features <- lapply(1:ncol(comb_loadings), function(f) {
  ord <- order(abs(comb_loadings[, f]), decreasing = TRUE)
  head(data.frame(
    Feature = rownames(comb_loadings)[ord],
    Weight  = comb_loadings[ord, f],
    stringsAsFactors = FALSE
  ), top_n)
})

names(top_features) <- paste0("Factor", 1:ncol(comb_loadings))
top_features[[1]]  # look at Factor1

#
scores <- get_factors(mofa_cov)
f <- 1  # factor to inspect
features_to_plot <- top_features[[f]]$Feature
#features_to_plot <- setdiff(features_to_plot, c("year2024", "year2022"))
colnames(result_fruit2) <- colnames(result_fruit2) %>%
  gsub("º", "deg", .) %>%       # replace degree symbol
  gsub("\\.", "_", .) %>%       # replace dots with underscore
  gsub("-", "_", .) %>%         # replace dashes with underscore
  gsub(" ", "_", .)              # replace spaces with underscore


colnames(result) <- colnames(result) %>%
  gsub("º", "deg", .) %>%       # replace degree symbol
  gsub("\\.", "_", .) %>%       # replace dots with underscore
  gsub("-", "_", .) %>%         # replace dashes with underscore
  gsub(" ", "_", .)              # replace spaces with underscore

#Leaf
combined <- bind_cols(spectra_pcs[,c("PC1","PC5","PC3","PC2")], result[,c(
                                                                     "Chlorophyll_b__mg_gFM_",
                                                                     "Chlorophyll_a__mg_gFM_",
                                                                     "Carotenoids__mg_gFM_",
                                                                     "Anthocyanins__mg_gFM_",
                                                                     "Starch__mg_gFM_",
                                                                     "ROS_O2__ABS_gFM_",
                                                                     "Phenols__g_GAE_mL_",
                                                                     "yday",
                                                                     "Tmax__degC__roll7d_mean",
                                                                     "Water_P"
                                                                     )], genes_expr[,c("RCA", 
                                                                                       "LBCY")])

combined <- bind_cols(spectra_fruit_pcs[,c("PC6")], result_fruit2[,c("Brix_lab__%_", 
                                                              "Alpha_Amino_Nitrogen__mg_L_",
                                                              "HRmed_____roll7d_mean",
                                                              "Tmax__degC__lag1d",
                                                              "P__mm__roll7d_mean",
                                                              "P__mm__lag1d",
                                                              "Tmax__degC__roll7d_mean",
                                                              "Water_P",
                                                              "Ammoniacal_Nitrogen__mg_L_",
                                                              "Water_P_lag7daytol", 
                                                              "yday","Assimilable_Nitrogen__mg_L_")])
combined <- combined %>% rename(PC6 = `...1`)
plot_data <- data.frame(scores = scores$group1[, f], combined[, features_to_plot], year = result$year  )
plot_data_long <- melt(plot_data, id.vars = c("scores","year"), variable.name = "Feature", value.name = "Value")
# Converta as colunas problemáticas para character
plot_data_long$Feature <- as.character(plot_data_long$Feature)
plot_data_long$year <- as.character(plot_data_long$year)
f_value <- as.character(f)
x_label_final <- paste0("Factor ", f, " scores")
######### IMPORTANT !!!!!! ######################
ggplot(plot_data_long, aes(x = scores, y = Value)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~ Feature, scales = "free") +
  # CORREÇÃO APLICADA AQUI:
  labs(x= x_label_final,#paste0("Factor ", f, " scores"), 
       y = "Feature value",
       title = "Top-weighted features vs Factor 1 scores") +
  theme_minimal()

plot_factor_cor(mofa_cov)
library(corrplot)
factors <- get_factors(mofa_cov)
# The output of get_factors is a list; we want the first element (the matrix)
# and convert it to a data frame.
factors_df <- factors[[1]] %>% 
  as.data.frame()

# 3. Calculate the Correlation Matrix
# The cor() function calculates the pairwise correlation between the columns (factors)
cor_matrix <- cor(factors_df)

# 4. Generate the Correlation Plot (Heatmap)
# Use the corrplot function for a clean visual representation.
corrplot(
  cor_matrix, 
  method = "color",           # Use color for the cells
  type = "upper",             # Display only the upper triangle
  order = "hclust",           # Optional: reorder factors based on clustering
  addCoef.col = "black",      # Add correlation coefficients to the plot
  tl.col = "black",           # Color of the factor labels
  tl.srt = 45,                # Text label rotation
  diag = FALSE,               # Do not display correlation on the diagonal (it's always 1)
  main = "Leaf: Correlation Matrix of MOFA Factors" # Title
)


ggplot(result_fruit2, aes(x = yday, y = Alpha_Amino_Nitrogen__mg_L_, color = Tmax__degC__lag1d)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Day of the year", y = "Alpha_Amino_Nitrogen__mg_L_") + theme_bw()

result_fruit2 <- result_fruit2 %>%
  mutate(Water_P_bin = cut(Water_P, breaks = quantile(Water_P, probs = seq(0, 1, 0.33), na.rm = TRUE),
                           include.lowest = TRUE, labels = c("Low", "Medium", "High")))

ggplot(result_fruit2, aes(x = yday, y = Alpha_Amino_Nitrogen__mg_L_, color = Water_P)) +
  geom_point(aes(color = Water_P), size = 3) +
  #geom_smooth(aes(group = Water_P_bin), method = "lm", se = FALSE, linetype = "dashed") +
  scale_color_gradient(low = "red", high = "blue") +
  labs(x = "Day of Year", y = "Alpha_Amino_Nitrogen__mg_L_", color = "Water_P (MPa)") +
  theme_minimal()


# make evenly spaced daily series (NA allowed), or use aggregated daily mean
tmax_series <- ts(result_fruit2$Tmax__degC__lag1d, frequency = 7)    # or frequency = 365
water_series <- ts(result_fruit2$Water_P, frequency = 7)

ccf_res <- ccf(tmax_series, water_series, lag.max = 21, na.action = na.pass) 
plot(ccf_res)

lm(Water_P ~ Tmax__degC__lag1d, data = result_fruit2)

lm(Water_P ~ Tmax__degC__lag1d + Irrigation + P__mm__roll7d_mean, data = result_fruit2)
cor(result_fruit2$Water_P, result_fruit2$Tmax__degC__lag1d, use="complete.obs")

m1 <- lm(Water_P ~ Tmax__degC__lag1d * Irrigation + P__mm__roll7d_mean, data=result_fruit2)
summary(m1)

ggplot(result_fruit2, aes(x = yday)) +
  geom_point(aes(y = Tmax__degC__lag1d, color = "Tmax 1day lag (°C)")) +
  geom_point(aes(y = Water_P*30, color = "Water_P (scaled)")) +  # scale Water_P for plotting
  geom_col(aes(y = P__mm__roll7d_mean), alpha = 0.3, fill = "skyblue") +
  scale_y_continuous(
    name = "Tmax / Precip",
    sec.axis = sec_axis(~./30, name = "Water_P (MPa)")) +
  labs(color = "") + theme_minimal()

pc6_loadings_df <- as.data.frame(pc6_loadings)
pc6_loadings_df$wavelength <- as.numeric(sub("wl_", "", rownames(loadings)))
pc2_loadings_df <- as.data.frame(pc2_loadings)
pc2_loadings_df$wavelength <- as.numeric(sub("wl_", "", rownames(loadings)))
# Example code
ggplot() +
  geom_line(data = pc6_loadings_df, aes(x = wavelength, y = pc6_loadings, color = "PC6 (positive)"), size=1.2) +
  geom_line(data = pc2_loadings_df, aes(x = wavelength, y = pc2_loadings, color = "PC2 (negative)"), size=1.2, linetype="dashed") +
  labs(x="Wavelength (nm)", y="Absorbance loading (a.u.)",
       color="Spectral PC") +
  theme_minimal()

ggplot(scores, aes(x = scores, y = spectra_fruit_pcs$PC6, color = Irrigation)) +
  geom_point(size=3) + geom_smooth(method="lm", se=FALSE) +
  labs(x="Factor 1", y="Spectral PC6 (absorbance)") +
  theme_minimal()
#### 3. Relationship (Weights × Features = Scores)
# First factor column after 'feature'
factor_name <- setdiff(colnames(weights_df), "feature")[1]

w <- weights_df[, c("feature", factor_name), drop = FALSE]
colnames(w) <- c("feature", "Weight")

# Take feature matrix (scaled values for features)
# Assumes your original input data is available
X <- get_data(MOFA_pcs, views = "metabolites")[[1]]
X_df <- as.data.frame(t(X)) # features × samples

# Compute scores manually
manual_scores <- as.matrix(t(X_df)) %*% as.matrix(w$Weight)

# Compare with MOFA factor scores
cbind(
  MOFA_score = scores_df[[factor_name]],
  Manual_score = manual_scores
) %>% head()

weights_mat <- get_weights(mofa_cov, views = "spectra")[[1]]  # features × factors
scores <- get_factors(mofa_cov)
scores_mat <- as.matrix(scores$group1)
# Example for one view and one factor
factor_name <- "Factor1"
corrs <- apply(weights_mat, 1, function(w) 
  cor(w, scores_mat[, factor_name], use = "pairwise.complete.obs")
)

#Get reconstructed data for a given view
recon <- reconstruct_data(mofa_cov, views = "genes_expression")[[1]]  # samples × features

# Example: correlation of each feature’s values across samples with Factor1
factor_name <- "Factor1"
feature_corrs <- apply(recon, 2, function(x)
  cor(x, scores_mat[, factor_name], use = "pairwise.complete.obs")
)

hist(corrs, main = paste("Correlation of features with", factor_name, "scores"),
     xlab = "Correlation (r)", breaks = 30)

####Factor scores by covariates (boxplot or violin)
scores <- get_factors(mofa_cov, factors = "all", as.data.frame = TRUE)
#scores_annot <- merge(scores, result_fruit2[, c(, "Irrigation")])
scores_annot <- cbind(scores, Irrigation = result_fruit2$Irrigation, 
                      Year = result_fruit2$year, Test_site = result_fruit2$Test_site, 
                      Cultivar = result_fruit2$Cultivar)
library(ggpubr)

ggplot(scores_annot, aes(x = Irrigation, y = value, fill = Irrigation)) +
  geom_boxplot() +
  facet_wrap(~ factor, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "MOFA factors by irrigation treatment", y = "Factor score", x = "")

my_comparisons <- list(
  c("NI", "I")  # Same for all, or customize per facet
)
my_comparisons_year <- list(
  c("2021","2022"),
  c("2021","2023"),
  c("2021","2024"),
  c("2022","2023"),
  c("2022","2024"),
  c("2023","2024")
)

my_comparisons_site <- list(
  c("Quinta dos Aciprestes","Quinta do Seixo"),
  c("Quinta dos Aciprestes", "Quinta de Vale de Cavalos"),
  c("Quinta do Seixo", "Quinta de Vale de Cavalos")
)

my_comparisons_cultivar <- list(
  c("Touriga Nacional","Moscatel Galego Branco"))

ggplot(scores_annot, aes(x = Cultivar, y = value, fill = Cultivar)) +
  geom_boxplot() +
  facet_wrap(~factor, scales = "free_y") +  # <-- precisa do +
  theme_minimal(base_size = 14) +
  labs(title = "Factor scores by Cultivar",
       y = "Factor score", x = "") +
  theme(legend.position = "none") +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = my_comparisons_cultivar
  )

ggsave("Berry_covariates_Factor scores by Cultivar.pdf", width=20, height=12, dpi=300)

#########################
#Feature–score correlation plot
#explore how much each feature correlates with a factor’s scores.
#This is a way to visualize the biological strength of association between features and that factor.
##A skewed histogram (mostly positive or mostly negative) means the factor drives most features in one direction.
## A bimodal shape means opposing feature groups (e.g., genes up vs down in stress).
scores_all <- get_factors(mofa_cov, as.data.frame = TRUE)
# reshape from long to wide: one column per factor
scores_wide <- scores_all %>%
  select(sample, factor, value) %>%
  pivot_wider(names_from = factor, values_from = value)
weights <- get_weights(mofa_cov, views = "spectra", as.data.frame = TRUE)

X <- as.matrix(views[[view_name]])
w_vec <- W[, factor_name]
if (length(w_vec) != ncol(X)) stop("Weights length does not match number of features in view")
# Example for one view and one factor
library(tidyverse)

all_corrs <- tibble()

for (view_name in names(weights_list)) {
  W <- weights_list[[view_name]]            # features × factors
  X <- as.matrix(views[[view_name]])   # samples × features
  cat("Processing view:", view_name, "\n")
  
  # --- Align features by name ---
  common_feats <- intersect(rownames(W), colnames(X))
  
  if (length(common_feats) == 0) {
    warning(paste("No overlapping features for view", view_name))
    next
  }
  
  W <- W[common_feats, , drop = FALSE]
  X <- X[, common_feats, drop = FALSE]
  
  for (factor_name in colnames(W)) {
    w_vec <- W[, factor_name]
    contrib <- t(t(X) * w_vec)   # samples × features
    scores_vec <- scores[[factor_name]]
    
    # Compute correlation across samples
    corrs <- apply(contrib, 2, function(fvec)
      cor(fvec, scores_vec, use = "pairwise.complete.obs"))
    
    temp <- tibble(
      feature = colnames(X),
      view = view_name,
      factor = factor_name,
      correlation = corrs
    )
    all_corrs <- bind_rows(all_corrs, temp)
  }
}

# Optional: inspect results
print(all_corrs)


ggplot(all_corrs, aes(x = correlation)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_grid(view ~ factor, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Feature–Factor Correlations per View",
    x = "Correlation",
    y = "Count"
  )
##########################
#Cross-omics correlation of factor expressions
scores_all <- get_factors(mofa_cov, as.data.frame = TRUE)
# reshape from long to wide: one column per factor
scores_wide <- scores_all %>%
  select(sample, factor, value) %>%
  pivot_wider(names_from = factor, values_from = value)

cor_mat <- cor(scores_wide[ , -1], use = "pairwise.complete.obs")
cor_mat[upper.tri(cor_mat)] <- NA

pheatmap(cor_mat, main = "Berry + covariates: Cross-factor correlations across all views", 
         display_numbers = TRUE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         na_col = "white", 
         breaks = seq(-1, 1, length.out = 100))

### weights x weights cor plots
# 1️⃣ Extract weights for all views
weights <- get_weights(mofa_cov)

# 2️⃣ Compute correlation matrices for each view and store in tidy format
cor_long <- lapply(names(weights), function(view_name) {
  w <- as.data.frame(weights[[view_name]])
  
  # 🔍 Drop any NA or empty columns
  #w <- w[, colSums(!is.na(w)) > 0, drop = FALSE]
  cor_mat <- cor(w, use = "pairwise.complete.obs")
  
  # Turn into long format for ggplot
  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Factor1", "Factor2", "Correlation")
  cor_df$View <- view_name
  return(cor_df)
}) %>% bind_rows()

# Optional: ensure factors are ordered properly
cor_long <- cor_long %>%
  mutate(
    Factor1 = factor(Factor1, levels = paste0("Factor", 1:12)),
    Factor2 = factor(Factor2, levels = paste0("Factor", 1:12)),
    View = factor(View, levels = c("genes_expression", "metabolites", "spectra", "others"))
  )

# 3️⃣ Plot all correlation matrices in one faceted heatmap
ggplot(cor_long, aes(x = Factor1, y = Factor2, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), size = 2.5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1)
  ) +
  facet_wrap(~View, ncol = 2) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pairwise correlations of MOFA factor weights across views",
    fill = "Correlation",
    x = NULL, y = NULL
  )

##Scores vs Scores (Sample scatterplots)
#It shows how samples distribute along two latent dimensions (factors).

# Example: plot Factor1 vs Factor2
ggplot(scores_wide, aes(x = Factor1, y = Factor2, color = Irrigation)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Sample distribution: Factor 1 vs Factor 2",
       x = "Factor 1 score", y = "Factor 2 score") +
  theme(legend.position = "bottom")

# assuming row order of MOFA samples == order in result_fruit2
#scores_all$sample_index <- seq_len(nrow(scores))
result_fruit2$sample <- seq_len(nrow(result_fruit2))
result2_noOut_starch2$sample <- seq_len(nrow(result2_noOut_starch2))

scores_annot <- merge(scores_wide, 
                      result_fruit2[, c("sample", "Irrigation",
                                        "year","Cultivar", "Test_site",
                                        "yday","FW_sample")])

scores_annot <- merge(scores_wide, 
                      result2_noOut_starch2[, c("sample", "Irrigation",
                                        "year","Cultivar", "Test_site",
                                        "yday","FW_sample","Chlorophyll_b__mg_gFM_",
                                        "Chlorophyll_a__mg_gFM_","CHLG (Cq)","LBCY (Cq)","RCA (Cq)")])


ggplot(scores_annot, aes(x = Factor1, y = Factor2, color = Irrigation)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "Sample distribution: Factor 1 vs Factor 2",
       x = "Factor 1 score", y = "Factor 2 score") +
  theme(legend.position = "bottom")

ggplot(scores_annot, aes(x = year, y = Factor2, group = Cultivar, color = Cultivar)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(title = "Temporal trend of Factor 2 across years",
       y = "Factor 2 score")

cor_mat <- cor(scores_annot[, grep("Factor", names(scores_annot))], 
               cov_block_select_fruit, use = "pairwise.complete.obs")

pheatmap::pheatmap(cor_mat,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   breaks = seq(-1, 1, length.out =100),
                   main = "Correlation of factor scores with numeric covariates")

scores_long <- scores_annot %>%
  pivot_longer(
    cols = starts_with("Factor"),
    names_to = "Factor",
    values_to = "Score"
  ) %>%
  # Ensure Factor levels are ordered numerically
  mutate(
    Factor = factor(
      Factor,
      levels = paste0("Factor", 1:9)  # adjust if you have more/less
    )
  )

# Plot all factors at once
ggplot(scores_long, aes(x = year, y = Score, group = Cultivar, color = Cultivar)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Factor, scales = "fixed", ncol = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Temporal trends of MOFA Factors across years and cultivars",
    x = "Year",
    y = "Factor score"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  )
#Within-year dynamics
ggplot(scores_long, aes(x = yday, y = Score, color = Cultivar)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +  # smooth trend lines
  facet_wrap(~ Factor, scales = "fixed", ncol = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Continuous temporal dynamics of MOFA Factors across day of year",
    x = "Day of Year (yday)",
    y = "Factor score"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  )

ggplot(scores_long, aes(x = yday, y = Score,
                        color = Irrigation, shape = Test_site)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(aes(group = Irrigation), method = "loess",
              se = FALSE, span = 0.8, linewidth = 0.9) +
  #scale_shape_manual(values = c(22,24)) +
  scale_shape_manual(values = c(21, 24, 22)) +  # 21 = circle, 24 = triangle (both hollow)
  facet_wrap(~ Factor, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Temporal dynamics of MOFA factors across day of year",
    subtitle = "Rows = Test_site, Columns = Factors, Point shapes = Irrigation",
    x = "Day of Year (yday)",
    y = "Factor Score",
    color = "Test_site",
    shape = "Irrigation"
  ) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold", size = 11),
    strip.text.y = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )
#### keep this one on hold
ggplot(scores_long, aes(x = yday, y = Score,
                        fill = Cultivar, color = Irrigation, linetype = Irrigation)) +
  geom_point(shape = 21, size = 2.8, stroke = 1.2, alpha = 0.8) +
  geom_smooth(aes(group = interaction(Cultivar, Irrigation)),
              method = "loess", se = FALSE, span = 0.8,
              linewidth = 1, na.rm = TRUE,
              method.args = list(surface = "direct")) +
  facet_wrap(~ Factor, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Irrigated" = "#0072B2", "Non-irrigated" = "#D55E00")) +
  scale_fill_brewer(palette = "Set2") +
  scale_linetype_manual(values = c("Irrigated" = "solid", "Non-irrigated" = "dashed")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Temporal dynamics of MOFA factors across day of year",
    subtitle = "Fill = Cultivar; Outline & line = Irrigation",
    x = "Day of Year (yday)",
    y = "Factor Score",
    fill = "Cultivar",
    color = "Irrigation",
    linetype = "Irrigation"
  ) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold", size = 11),
    strip.text.y = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

scores_corr <- scores_annot %>%
  select(starts_with("Factor")) %>%
  cor(method = "pearson", use = "pairwise.complete.obs")

scores_yearly <- scores_annot %>%
  group_by(year) %>%
  summarize(across(starts_with("Factor"), mean, na.rm = TRUE))

scores_corr <- scores_yearly %>%
  select(starts_with("Factor")) %>%
  cor(method = "pearson", use = "pairwise.complete.obs")

library(ggcorrplot)
ggcorrplot(scores_corr, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE,          # <-- add correlation values as text labels
           lab_size = 3)

# Perform hierarchical clustering based on correlation distance
dist_mat <- as.dist(1 - scores_corr)  # convert correlation to distance
hc <- hclust(dist_mat, method = "average")

# Plot the dendrogram
plot(hc, main = "Temporal clustering of MOFA factors (Vitis vinifera) across years",
     xlab = "", sub = "", cex = 0.8)

ggcorrplot(
  scores_corr,
  hc.order = TRUE,
  type = "lower",
  lab = TRUE,
  outline.col = "white",
  ggtheme = ggplot2::theme_minimal(),
  title = "Temporal correlation between MOFA factors (Vitis vinifera) across years",
  colors = c("blue", "white", "red"))
  
clusters <- cutree(hc, k = 4)  # choose 3 clusters, adjust as needed
  
  cluster_df <- data.frame(Factor = names(clusters), Cluster = clusters)
  
  ggplot(cluster_df, aes(x = Factor, y = Cluster, fill = as.factor(Cluster))) +
    geom_tile() +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    ggtitle("Clusters of temporally linked MOFA factors")

# Plot 1: CHLG vs Chlorophyll_a
ggplot(scores_annot, aes(x = `CHLG (Cq)`, y = Chlorophyll_a__mg_gFM_)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "CHLG expression vs Chlorophyll a",
       x = "CHLG expression",
       y = "Chlorophyll a content") +
  theme_minimal()

df_long_leaf <- result2_noOut_starch2 %>%
  pivot_longer(-Test_site, names_to = "Feature", values_to = "Value") %>%
  mutate(OmicType = case_when(
    Feature %in% c("MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                   "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
                   "LBCY (Cq)","CHLG (Cq)") ~ "transcriptome",
    Feature %in% c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                   "Phenols__g_GAE_mL_","ROS__mmol_gFM_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                   "Sugar__mg_gFM_","Anthocyanins__mg_gFM_") ~ "metabolome",
    TRUE ~ "other"
  ))

df_long_leaf <- df_leaf_joined_complete_lagged %>%
  pivot_longer(
    cols = c("RCA (Cq)",
             "LBCY (Cq)","CHLG (Cq)","Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
             "Phenols__g_GAE_mL_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
             "Sugar__mg_gFM_","Anthocyanins__mg_gFM_"),
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  mutate(
    OmicType = case_when(
      Feature %in% c("RCA (Cq)",
                     "LBCY (Cq)","CHLG (Cq)") ~ "transcriptome",
      Feature %in% c("Carotenoids__mg_gFM_","Chlorophyll_a__mg_gFM_","Chlorophyll_b__mg_gFM_",
                     "Phenols__g_GAE_mL_","ROS_O2__ABS_gFM_","Starch__mg_gFM_",
                     "Sugar__mg_gFM_","Anthocyanins__mg_gFM_") ~ "metabolome"
    )
  )

table(df_long_leaf$OmicType, df_long_leaf$Test_site)
table(df_long_leaf$OmicType, df_long_leaf$Code_field)

tab_df_plant <- as.data.frame.matrix(table(df_long_leaf$OmicType, df_long_leaf$Code_field))

tab_long_plant <- as.data.frame(table(df_long_leaf$OmicType, df_long_leaf$Code_field))
colnames(tab_long_plant) <- c("OmicType", "Code_field", "Count")

ggplot(tab_long_plant, aes(x = Code_field, y = Count, fill = OmicType)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  coord_flip()

df_long_leaf %>%
  group_by(Test_site) %>%
  summarise(n_omics = n_distinct(OmicType))

df_long_leaf %>%
  filter(Feature == "CHLG (Cq)") %>%
  group_by(Test_site) %>%
  summarise(
    n_total = n(),
    n_nonNA = sum(!is.na(Value)),
    n_NA = sum(is.na(Value))
  )

df_long_leaf %>%
  dplyr::filter(Feature %in% c("CHLG (Cq)", "Chlorophyll_a__mg_gFM_")) %>%
  dplyr::select(Code_field, Test_site, Feature, Value) %>%
  dplyr::filter(!is.na(Value)) %>%
  pivot_wider(names_from = Feature, values_from = Value) %>%
  mutate(has_both = ifelse(!is.na(`CHLG (Cq)`) & !is.na(`Chlorophyll_a__mg_gFM_`), TRUE, FALSE)) %>%
  count(Test_site, has_both)


df_long_leaf %>%
  filter(Feature %in% c("CHLG (Cq)", "Chlorophyll_a__mg_gFM_")) %>%
  select(Test_site, Date, Code_field, Feature, Value) %>%
  pivot_wider(names_from = Feature, values_from = Value) %>%
  mutate(has_both = ifelse(!is.na(`CHLG (Cq)`) & !is.na(`Chlorophyll_a__mg_gFM_`), TRUE, FALSE)) %>%
  count(Test_site, Date, has_both)

df_clean_leaf <- df_long_leaf %>%
  filter(Feature %in% c("CHLG (Cq)", "Chlorophyll_a__mg_gFM_")) %>%
  group_by(Code_field, Test_site, Date, Feature) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Feature, values_from = Value) %>%
  filter(!is.na(`CHLG (Cq)`), !is.na(`Chlorophyll_a__mg_gFM_`))

ggplot(df_clean_leaf, aes(x = 2^(-`LBCY (Cq)`), y = `Chlorophyll_a__mg_gFM_`, color = Test_site)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "Relationship between CHLG expression and Chlorophyll a content",
    subtitle = "Only overlapping samples (both CHLG and Chlorophyll measured)",
    x = "CHLG expression (Cq)",
    y = "Chlorophyll a (mg/g FM)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggplot(df_long_leaf %>% filter(Feature %in% c("CHLG (Cq)", "Chlorophyll_a__mg_gFM_")),
       aes(x = Date, y = Value, color = Feature)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = FALSE) +
  facet_wrap(~ Test_site, scales = "free_y") +
  theme_minimal() +
  labs(title = "Temporal profiles of CHLG expression and Chlorophyll a",
       subtitle = "Each site over time",
       y = "Value (Cq or mg/g FM)")

df_long_leaf %>%
  filter(Feature %in% c("CHLG (Cq)", "Chlorophyll_a__mg_gFM_")) %>%
  group_by(Test_site, Feature) %>%
  summarise(median_date = median(Date),
            n_timepoints = n_distinct(Date))

library(lubridate)
df_long_leaf %>%
  mutate(Year = year(Date)) %>%
  group_by(Test_site, Year, Feature) %>%
  summarise(
    median_date = median(Date, na.rm = TRUE),
    n_timepoints = n_distinct(Date),
    .groups = "drop"
  ) %>%
  arrange(Test_site, Year, Feature)

df_long_leaf %>%
  # keep only CHLG and chlorophylls
  filter(Feature %in% c("LBCY (Cq)", "RCA (Cq)", "CHLG (Cq)",
                        "Chlorophyll_a__mg_gFM_", 
                        "Carotenoids__mg_gFM_","ROS_O2__ABS_gFM_",
                        "Phenols__g_GAE_mL_", "Starch__mg_gFM_", "Sugar__mg_gFM_")) %>%
  mutate(Year = year(Date)) %>%
  ggplot(aes(x = Date, y = Value, color = Feature)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(se = FALSE, span = 0.8) +
  facet_grid(Test_site ~ Year, scales = "free_x") +
  labs(
    title = "Temporal profiles of CHLG expression and Chlorophyll a/Carotenoids per year",
    x = "Sampling date",
    y = "Value (Cq or mg/g FM)",
    color = "Feature"
  ) +
  theme_minimal(base_size = 13)


df_plot <- df_long_leaf %>%
  filter(Feature %in% c("CHLG (Cq)", "LBCY (Cq)", "RCA (Cq)",
                        "ROS__mmol_gFM_", "ROS_O2__ABS_gFM_",
                        "Phenols__g_GAE_mL_", "Starch__mg_gFM_", "Sugar__mg_gFM_")) %>%
  mutate(
    Year = year(Date),
    # Convert transcript features to expression
    Value = if_else(Feature %in% c("CHLG (Cq)", "LBCY (Cq)", "RCA (Cq)"),
                    2^(-Value), Value)
  )

# Compute automatic scale factor to align both y-axis groups visually
range_transcript <- range(df_plot$Value[df_plot$Feature %in%
                                          c("CHLG (Cq)")],
                          na.rm = TRUE)

range_metabolite <- range(df_plot$Value[df_plot$Feature %in%
                                          c("ROS_O2__ABS_gFM_", "ROS__mmol_gFM_"
                                            )],
                          na.rm = TRUE)

# Scale factor = how much to multiply metabolite values by so that they fit
scale_factor <- diff(range_transcript) / diff(range_metabolite)

ggplot() +
  # Metabolite features (scaled for dual axis)
  geom_point(data = df_plot %>%
               filter(Feature %in% c("ROS__mmol_gFM_", "ROS_O2__ABS_gFM_"
                                     )),
             aes(x = Date, y = Value * scale_factor, color = Feature),
             size = 2, alpha = 0.7) +
  geom_smooth(data = df_plot %>%
                filter(Feature %in% c("ROS__mmol_gFM_", "ROS_O2__ABS_gFM_")),
              aes(x = Date, y = Value * scale_factor, color = Feature),
              se = FALSE, span = 0.8) +
  
  # Transcript features (primary axis)
  geom_point(data = df_plot %>%
               filter(Feature %in% c( "CHLG (Cq)")),
             aes(x = Date, y = Value, color = Feature),
             size = 2, alpha = 0.7) +
  geom_smooth(data = df_plot %>%
                filter(Feature %in% c("CHLG (Cq)")),
              aes(x = Date, y = Value, color = Feature),
              se = FALSE, span = 0.8) +
  
  facet_grid(Test_site ~ Year, scales = "free_x") +
  scale_y_continuous(
    name = "Transcript expression (2^-Cq)",
    sec.axis = sec_axis(~ . / scale_factor,
                        name = "metabolite concentration (ABS/g FM)")
  ) +
  labs(
    title = "Transcript–metabolite dynamics (CHLG vs ROS)",
    x = "Sampling date",
    color = "Feature"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  )

df_long_leaf <- df_long_leaf %>%
  mutate(Year = year(Date))
df_chlg_chl <- df_long_leaf %>%
  filter(Feature %in% c("CHLG (Cq)", "Chlorophyll_a__mg_gFM_", "Chlorophyll_b__mg_gFM_")) %>%
  mutate(Expression = ifelse(Feature == "CHLG (Cq)", 2^(-Value), Value)) %>%
  select(Test_site, Date, Year, Feature, Expression)

# Pivot wide to have columns for each variable
df_wide <- df_chlg_chl %>%
  pivot_wider(names_from = Feature, values_from = Expression)
df_wide_clean <- df_wide %>%
  mutate(
    `CHLG (Cq)` = map_dbl(`CHLG (Cq)`, ~ mean(as.numeric(.x), na.rm = TRUE)),
    Chlorophyll_a__mg_gFM_ = map_dbl(Chlorophyll_a__mg_gFM_, ~ mean(as.numeric(.x), na.rm = TRUE))
  )

df_wide_clean %>%
  group_by(Test_site, Year) %>%
  summarise(
    n_non_na = sum(!is.na(`CHLG (Cq)`) & !is.na(Chlorophyll_a__mg_gFM_))
  )
# Compute correlation and cross-correlation per site-year
library(purrr)

xcorr_results <- df_wide_clean %>%
  group_by(Test_site, Year) %>%
  summarise(
    # keep only complete pairs for both variables
    data = list(na.omit(select(cur_data(), `CHLG (Cq)`, Chlorophyll_a__mg_gFM_))),
    .groups = "drop"
  ) %>%
  mutate(
    cor_now = map_dbl(data, ~ if (nrow(.x) > 1)
      cor(.x$`CHLG (Cq)`, .x$Chlorophyll_a__mg_gFM_)
      else NA_real_),
    
    xcorr = map(data, ~ if (nrow(.x) > 1)
      ccf(.x$`CHLG (Cq)`, .x$Chlorophyll_a__mg_gFM_, plot = FALSE)
      else NA),
    
    lag = map_dbl(xcorr, ~ if (is.list(.x))
      .x$lag[which.max(abs(.x$acf))]
      else NA_real_),
    
    max_corr = map_dbl(xcorr, ~ if (is.list(.x))
      .x$acf[which.max(abs(.x$acf))]
      else NA_real_)
  )


xcorr_results_site <- df_wide_clean %>%
  group_by(Test_site) %>%
  summarise(
    data = list(na.omit(select(cur_data(), `CHLG (Cq)`, Chlorophyll_a__mg_gFM_))),
    cor_now = map_dbl(data, ~ if (nrow(.x) > 1)
      cor(.x$`CHLG (Cq)`, .x$Chlorophyll_a__mg_gFM_) else NA_real_),
    xcorr = map(data, ~ if (nrow(.x) > 1)
      ccf(.x$`CHLG (Cq)`, .x$Chlorophyll_a__mg_gFM_, plot = FALSE) else NA)
  )

xcorr_results_site <- df_wide_clean %>%
  group_by(Test_site) %>%
  summarise(
    data = list(na.omit(select(cur_data(), `CHLG (Cq)`, Chlorophyll_a__mg_gFM_))),
    .groups = "drop"
  ) %>%
  mutate(
    cor_now = map_dbl(data, ~ if (nrow(.x) > 1)
      cor(.x$`CHLG (Cq)`, .x$Chlorophyll_a__mg_gFM_)
      else NA_real_),
    xcorr = map(data, ~ if (nrow(.x) > 1)
      ccf(.x$`CHLG (Cq)`, .x$Chlorophyll_a__mg_gFM_, plot = FALSE)
      else NA),
    lag = map_dbl(xcorr, ~ if (is.list(.x))
      .x$lag[which.max(abs(.x$acf))] else NA_real_),
    max_corr = map_dbl(xcorr, ~ if (is.list(.x))
      .x$acf[which.max(abs(.x$acf))] else NA_real_)
  )

xcorr_results %>% select(Test_site, Year, lag, max_corr, cor_now)

df_wide_clean %>%
  mutate(both_non_na = !is.na(`CHLG (Cq)`) & !is.na(Chlorophyll_a__mg_gFM_)) %>%
  group_by(Test_site, Year) %>%
  summarise(n_pairs = sum(both_non_na)) %>%
  arrange(desc(n_pairs))


ggplot(df_long_leaf %>% 
         filter(Feature %in% c("CHLG (Cq)", "LBCY (Cq)", "RCA (Cq)",
                               "Chlorophyll_a__mg_gFM_", "Chlorophyll_b__mg_gFM_")),
       aes(x = yday, y = Value, color = Feature)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE, span = 0.7) +
  facet_wrap(~ Test_site, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "Feature concentration vs. Day of Year",
       x = "Day of Year", y = "Value")

ggplot(df_long_leaf %>% filter(Feature == "CHLG (Cq)"),
       aes(x = yday, y = Value)) +
  geom_point(aes(color = Tmax..ºC._roll7d_mean )) +
  geom_smooth(se = FALSE) +
  scale_color_viridis_c(name = "Tmax (°C)") +
  facet_wrap(~ Test_site) +
  theme_minimal(base_size = 13) +
  labs(title = "CHLG expression vs. Day of Year, colored by temperature",
       x = "Day of Year", y = "CHLG (Cq)")

genes_df <- result[,c("MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
                      "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
                      "LBCY (Cq)","CHLG (Cq)","Date","Test_site","Water_P")] #"GAPDH (Cq)"
genes_df <- genes_df %>%
  filter(Test_site == "Quinta dos Aciprestes") %>%
  select(
    "MYBA (Cq)","MYB14 (Cq)","MYB15 (Cq)","ABCC1 (Cq)","C4H1 (Cq)","CHS1 (Cq)",
    "DFR (Cq)","PAL1 (Cq)","MATE1 (Cq)","UFGT1 (Cq)","FLS1 (Cq)","RCA (Cq)",
    "LBCY (Cq)","CHLG (Cq)","Date","Water_P"
  ) %>%
  mutate(across(c(-Date,-Water_P), ~ 2^(-.)))
  
metabs_df <- result[,c("Chlorophyll_a__mg_gFM_", "Chlorophyll_b__mg_gFM_",
                       "Carotenoids__mg_gFM_", "ROS_O2__ABS_gFM_",
                       "Phenols__g_GAE_mL_", "Starch__mg_gFM_", 
                       "Anthocyanins__mg_gFM_", "Sugar__mg_gFM_",
                       "Date","Test_site","Water_P")] 

metabs_df <- metabs_df %>%
  filter(Test_site == "Quinta dos Aciprestes") %>%
  select(
    "Chlorophyll_a__mg_gFM_", "Chlorophyll_b__mg_gFM_",
                       "Carotenoids__mg_gFM_", "ROS_O2__ABS_gFM_",
                       "Phenols__g_GAE_mL_", "Starch__mg_gFM_", 
                       "Anthocyanins__mg_gFM_", "Sugar__mg_gFM_",
                       "Date","Test_site","Water_P")


#df_climate_all <- df_climate_all %>%rename(Date = Date)

df <- df_climate_all %>%
  left_join(metabs_df, by = "Date")
# Scale factor for secondary axis
#scale_factor <- max(df$`Tmax (ºC)`, na.rm = TRUE) / max(df$`RCA (Cq)`, na.rm = TRUE)


############### PLOT weather and concentrations #####
library(ggtext)
precip_scale_factor <- 1
html_y_title <- "<span style='color:red;'>Max Temperature (°C)</span> / <span style='color:dodgerblue4;'>Scaled Precipitation (mm)</span>"
# Plot with monthly x-axis
ggplot(df, aes(x = Date)) +
  geom_col(
    aes(y = `P (mm)` / precip_scale_factor, fill = "Precipitation"),
    alpha = 0.6 # Use alpha for transparency so lines are visible
  ) +
  geom_line(aes(y = `Tmax (ºC)`, color = "Temperature"), size = 0.5) +
  geom_point(aes(y = `CHLG (Cq)` * scale_factor, color = "Gene Expression"), size = 2) +
  scale_y_continuous(
    name = " Temperature (°C) /Precipitation (mm)", #html_y_title
    sec.axis = sec_axis(~./scale_factor, name = "CHLG Expression")
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y"
  ) +
  scale_color_manual(
    name = "",
    values = c("Temperature" = "red", "Gene Expression" = "blue")
  ) +
  scale_fill_manual(
    name = "",
    values = c("Precipitation" = "dodgerblue4") 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.left = element_text(color = "black", size = 12),
    axis.title.y.right = element_text(color = "blue", size = 12),
    legend.position = "top",
    #axis.title.y.left = element_markdown(size = 12)
  ) +
  labs(x = "Month")


df <- df %>%
  mutate(
    Date = as.Date(Date), # Assegura que a coluna está no formato de data correto
    Year = factor(lubridate::year(Date)), # Extrai o ano para o faceting
    Month = lubridate::month(Date) # Extrai o número do mês para a filtragem
  )

# 2. FILTRAR DADOS: Apenas Junho (Mês 6) a Setembro (Mês 9)
df_filt <- df %>%
  filter(Month >= 6 & Month <= 9)

# 2. DEFINE THE HTML AXIS TITLE
html_y_title <- "<span style='color:red;'>Max Temperature (°C)</span> / <span style='color:dodgerblue4;'>Scaled Precipitation (mm)</span>"

ggplot(df_filt, aes(x = Date)) +
  # Precipitation Bars
  geom_col(
    aes(y = `P (mm)` / precip_scale_factor, fill = "Precipitation"),
    alpha = 0.6
  ) +
  
  # Existing Temperature and Gene Expression layers
  geom_line(aes(y = `Tmax (ºC)`, color = "Temperature"), size = 1.2) +
  #geom_line(aes(y = `Tmax (ºC)`), colour = "red", size = 1.2) +
  #geom_line(aes(y = `CHLG (Cq)` * scale_factor, color = Irrigation), size = 1.2) +
  geom_point(aes(y = `RCA (Cq)` * scale_factor, color = "Gene Expression"), size = 2) +
  #Gene Expression
  # 3. ADD FACET WRAP!
  # This creates a separate plot for every unique 'Year' value
  facet_wrap(~Year, scales = "free_x", nrow = 2) + # Use 'free_x' so each facet only shows its own dates
  
  # Define the Y-AXES
  scale_y_continuous(
    name = " Temperature (°C) /Precipitation (mm)", #html_y_title,
    sec.axis = sec_axis(~./scale_factor, name = "RCA Expression")
  ) +
  
  # Keep X-axis labels monthly, but now they reset per facet
  scale_x_date(
    date_breaks = "3 months", # Often better for monthly facets
    date_labels = "%b"        # Show only month abbreviation
  ) +
  # 4. MUDANÇA PRINCIPAL: Ajustar o Eixo X para mostrar Dia e Mês
  scale_x_date(
    date_breaks = "7 days", # Quebrar a cada 10 dias
    date_labels = "%d %b"    # Formato: Dia (dd) Mês (Abrev.)
  ) +
  # Existing Scales
  scale_color_manual(name = "", values = c("Temperature" = "red", "Gene Expression" = "purple")) +
  scale_fill_manual(name = "", values = c("Precipitation" = "dodgerblue4")) +
  
  # Remaining theme and guides
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_minimal() +
  theme(
    # Custom ggtext rendering for the left axis title
    #axis.title.y.left = element_markdown(size = 12),
    # Adjust X-axis text angle as needed
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "purple", size = 12),
    axis.title.y.left = element_text(color = "black", size = 12),
    legend.position = "top",
    # Style the facet labels (the year headings)
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(x = "Month")

df_filt_scale <- df_filt %>%
  mutate(
    # RCA * 10
    `RCA (Cq) x 10` = `RCA (Cq)` * 10,
    # LBCY * 100
    `LBCY (Cq)` = `LBCY (Cq)` * 150,
    # CHLG mantido inalterado (ou criar uma coluna para consistência)
    `CHLG (Cq) x 100` = `CHLG (Cq)` * 150
  )

max_tmax <- max(df_filt_scale$`Tmax (ºC)`, na.rm = TRUE)

# 2. Identificar o valor máximo global entre TODOS os genes (Escala B)
# Usamos max() no vetor combinado (c()) dos valores absolutos dos genes
# O abs() garante que considera a magnitude se houver valores negativos.
max_gene_cq <- max(
  c(
    abs(df_filt_scale$`CHLG (Cq) x 100`), 
    abs(df_filt_scale$`RCA (Cq) x 10`), 
    abs(df_filt_scale$`LBCY (Cq)`)
  ), 
  na.rm = TRUE
)

max_metabs <- max(df_filt$ROS_O2__ABS_gFM_, na.rm=TRUE)

# 3. Calcular o Fator de Escala
# Se max_gene_cq for 0, evitamos a divisão por zero.
if (max_metabs == 0) {
  scale_factor <- 1
} else {
  scale_factor <- max_tmax / max_metabs#max_gene_cq
}
ggplot(df_filt, aes(x = Date)) +
  
  # 1. PRECIPITAÇÃO E TEMPERATURA
  geom_col(
    aes(y = `P (mm)` / precip_scale_factor, fill = "Precipitation"),
    alpha = 0.6
  ) +
  geom_line(aes(y = `Tmax (ºC)`, color = "Temperature"), size = 1.2) +
  # --- Linha Horizontal (NOVA) ---
  geom_hline(yintercept = 40, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 35, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_line(aes(y = `HRmed (%)`/1.3, color = "HR (%)"), linewidth = 1.2) +
  
  geom_point(aes(y = `Water_P`*10, color = "Water Base Potential"), size = 2, shape=1, alpha=0.4 )+
  geom_smooth(aes(y = `Water_P`*10 , color = "Water Base Potential",fill = "Water Base Potential"), method = "loess", se = TRUE, linewidth = 1,  alpha = 0.3) + 
  # --- PONTOS E TENDÊNCIA PARA CADA GENE (NOVA ESTRUTURA) ---
  
  # 2. CHLG (Cq)
  #geom_point(aes(y = `CHLG (Cq) x 100` * scale_factor, color = "CHLG (Cq)"), size = 2, shape = 1, alpha = 0.4) +
  #geom_smooth(aes(y = `CHLG (Cq) x 100` * scale_factor, color = "CHLG (Cq)",fill = "CHLG (Cq)"), method = "loess", se = TRUE, linewidth = 1, ,alpha = 0.3) +
  
  # 3. RCA (Cq)
  #geom_point(aes(y = `RCA (Cq) x 10` * scale_factor, color = "RCA (Cq)"), size = 2, shape = 1, alpha = 0.4) +
  #geom_smooth(aes(y = `RCA (Cq) x 10` * scale_factor, color = "RCA (Cq)",fill = "RCA (Cq)",), method = "loess", se = TRUE, linewidth = 1,  alpha = 0.3) +
  
  # 4. LBCY (Cq)
  #geom_point(aes(y = `LBCY (Cq)` * scale_factor, color = "LBCY (Cq)"), size = 2, shape = 1, alpha = 0.4) +
  #geom_smooth(aes(y = `LBCY (Cq)` * scale_factor, color = "LBCY (Cq)",fill = "LBCY (Cq)"), method = "loess", se = TRUE, linewidth = 1, , alpha = 0.3) +
  
geom_point(aes(y = `ROS_O2__ABS_gFM_` * scale_factor, color = "ROS"), size = 2, shape = 1, alpha = 0.4) +
  geom_smooth(aes(y = `ROS_O2__ABS_gFM_` * scale_factor, color = "ROS",fill = "ROS"), method = "loess", se = TRUE, linewidth = 1, alpha = 0.3) +
  # --- AJUSTES FINAIS E ESCALAS ---
  
  # FACET WRAP
  facet_wrap(~Year, scales = "free_x", nrow = 2) + 
  
  # Define the Y-AXES
  scale_y_continuous(
    name = " Temperature (°C) / Precipitation (mm) / Relative Humidity/1.5 (%) / BWPx10 (MPa)",
    breaks = seq(from = -20, to = 70, by = 5),
    # Usa GENE_SCALE_FACTOR para garantir que o rótulo do eixo secundário está correto
    #sec.axis = sec_axis(~. / scale_factor, name = "Gene expression (scaled)")
    sec.axis = sec_axis(~. / scale_factor, breaks = seq(from = 0, to = 8, by = 2),
                        name = "Gene expression (scaled)")
  ) +
  
  # Ajustar o Eixo X para mostrar Dia e Mês
  scale_x_date(
    date_breaks = "7 days", 
    date_labels = "%d %b"
  ) +
  
  # ESCALAS DE COR: Incluir todos os 3 genes e a Temperatura
  scale_color_manual(
    name = "", 
    values = c(
      "Temperature" = "red", 
      "HR (%)" = "grey",
      "Water Base Potential" = "lightblue",
      "ROS" = "orange"
      #"CHLG (Cq)" = "blue", 
      #"RCA (Cq)" = "darkgreen", 
      #"LBCY (Cq)" = "hotpink"
    )
  ) +
  #scale_fill_manual(name = "", values = c("Precipitation" = "dodgerblue4")) +
  scale_fill_manual(
    name = "", 
    values = c(
      "Precipitation" = "dodgerblue4",
      "Water Base Potential" = "lightblue",
      "ROS" = "orange"
     # "CHLG (Cq)" = "blue",          # Cor de preenchimento do CHLG
      #"RCA (Cq)" = "darkgreen",      # Cor de preenchimento do RCA
      #"LBCY (Cq)" = "hotpink"        # Cor de preenchimento do LBCY
    ),
    # Mostra apenas a Precipitação na lenda FILL
    breaks = c("Precipitation") 
  ) +
  # Remaining theme and guides
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
  theme_minimal() +
  coord_cartesian(ylim = c(-20, 70)) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Cor do título do eixo secundário (um neutro)
    axis.title.y.right = element_text(color = "black", size = 12),
    axis.title.y.left = element_text(color = "black", size = 12),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12)
  ) +
  labs(x = "Date")

