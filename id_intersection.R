setwd("/Users/joanaleite/Downloads/data_omicbots")

rdata_files <- list.files(pattern = "\\.RData$")
rdata_files

for (file in rdata_files) {
  load(file)  # loads all objects from the file into the global environment
}

obj_list <- c("Total_Anthocyanins__missing_unit_DATA.RData", "Tartaric_acid__g_L_objs_DATA.RData", "Tannins__mg_gFM__Seed_objs_DATA.RData", "Tannins__g_GAE_mL_1__Seed_objs_DATA.RData"
              , "Sugar__mg_gFM__Leaf_objs_DATA.RData","Starch__mg_gFM_objs_DATA.RData","ROS_DATA.RData","Potential_Alcohol__Vol._Vol.__Fruit_DATA.RData","Phenols__g_GAE_mL_objs_DATA.RData",
              "Malic_acid__g_L_DATA.RData","Hormones_LEAF_DATA.RData","Glucose___Fructose__g_L_objs_DATA.RData","Genes_LEAF_FRUIT_DATA.RData","Chlorophyll_b__mg_gFM_objs_DATA.RData",
              "Chlorophyll_a__mg_gFM_objs_DATA.RData","Carotenoids__mg_gFM_objs_DATA.RData","Brix_lab_objs_DATA.RData","BRIX_DATA.RData","Assimilable_Nitrogen__mg_L__Fruit_objs_DATA.RData",
              "Anthocyanins__mg_gFM_objs_DATA.RData","Ammoniacal_Nitrogen__mg_L_DATA.RData","Alpha_Amino_Nitrogen__mg_L__Fruit_objs_DATA.RData")

objets <- list()
for (f in obj_list) {
  e <- new.env()       # cria um ambiente temporário
  load(f, envir = e)   # carrega o RData nesse ambiente
  objets[[f]] <- e    # adiciona à lista
}

# Extrair vetores de IDs de cada objeto
vetors_ids <- lapply(objets, function(env) {
  obj <- env[[ls(env)[1]]]
  
  if (is.data.frame(obj)) {
    return(obj$ID)
  } else {
    return(obj)
  }
})

vectors_codes <- lapply(objets, function(env) {
  obj <- env[[ls(env)[1]]]
  
  if (is.data.frame(obj)) {
    return(obj$Code_field)
  } else {
    return(obj)
  }
})


common_ids <- Reduce(intersect, vetors_ids)
common_ids

# Juntar todos os IDs em um único vetor com duplicatas
todos_ids <- unlist(vetores_ids)
todos_codes <- unlist(vectors_codes)

vetores_ids_unique <- lapply(vetores_ids, unique)
vetores_codes_unique <- lapply(vectors_codes, unique)
todos_ids_unique <- unlist(vetores_ids_unique)
todos_codes_unique <- unlist(vetores_codes_unique)

# Criar tabela contando em quantos objetos cada ID aparece
id_counts <- table(todos_ids)
codes_counts <- table(todos_codes)

id_counts <- table(todos_ids_unique)
codes_counts <- table(todos_codes_unique)

id_counts_df <- as.data.frame(id_counts)
colnames(id_counts_df) <- c("id", "count")

codes_counts_df <- as.data.frame(codes_counts)
colnames(codes_counts_df) <- c("code_field", "count")

id_counts_df
codes_counts_df
