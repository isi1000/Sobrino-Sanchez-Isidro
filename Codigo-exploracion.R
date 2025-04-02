#IMPORTAR DATOS
dataset= read.csv("DataValues_S013.csv", sep=',')
#head(dataset)
dim(dataset)

#QUEDARSE SOLO CON T0 Y T5
library(dplyr)
dataset_T0 <- dataset %>% select(contains("T0"))
dataset_T5<-dataset%>% select(contains("T5"))
dim(dataset_T0)
dim(dataset_T5)

#localizar donde está la primera columna de metabolitos
index_ile=which(names(dataset_T0) == "Ile_T0")
index_ile
#filtramos y creamos el dataset de muestras y el de metabolitos
##T0
T0_sample=dataset_T0[,1:(index_ile-1)]
dataset_T0=dataset_T0[,index_ile:172]
##T5
T5_sample=dataset_T5[,1:(index_ile-1)]
dataset_T5=dataset_T5[,index_ile:172]

dim(T0_sample)
dim(T5_sample)
dim(dataset_T0)
dim(dataset_T5)
#extraer datos pacientes
sample_dataset=dataset[,3:6]


#unir con los datos temporales
sample_dataset=cbind(sample_dataset,T0_sample,T5_sample)
sample_dataset[,1:4] <- lapply(sample_dataset, as.factor)

dim(sample_dataset)
#head(sample_dataset)

#TRASPONER
matriz_T0 <- t(as.matrix(dataset_T0))
matriz_T5 <- t(as.matrix(dataset_T5))
#cambiamos los nombres de las columnas porque se han perdido
colnames(matriz_T0)=rownames(sample_dataset)
colnames(matriz_T5)=rownames(sample_dataset)
#head(matriz_T0)
#eliminar sufijo
rownames(matriz_T0)=gsub("_T0$", "", rownames(matriz_T0))
rownames(matriz_T5)=gsub("_T5$", "", rownames(matriz_T5))
#juntarlo en una lista
assays_list <- list(T0 = matriz_T0,T5 = matriz_T5)
#construimos la matriz y nos aseguramos de que los nombres coincidan
RowDataframe=data.frame(metabolite=rownames(matriz_T0))
rownames(RowDataframe)=rownames(matriz_T0)
ColDataframe=sample_dataset
rownames(ColDataframe)=colnames(matriz_T0)
library(SummarizedExperiment)

#construir objeto SE
met_exp <- SummarizedExperiment(
  assays = assays_list,
  rowData = RowDataframe,
  colData = ColDataframe)
met_exp
#calcular valores faltantes  en T5
valores_faltantes_T5= colSums(is.na(assay(met_exp, "T5")))
faltantes_T5= data.frame(
  name  = names(valores_faltantes_T5),
  count = valores_faltantes_T5,
  group = "T5"
)

#calcular valores faltantes en T0
valores_faltantes_T0= colSums(is.na(assay(met_exp, "T0")))
faltantes_T0= data.frame(
  name  = names(valores_faltantes_T0),
  count = valores_faltantes_T0,
  group = "T0"
)

#unir ambos 
faltantes= rbind(faltantes_T0, faltantes_T5)


#graficar con ggplot y separar por grupo
library(ggplot2)

ggplot(faltantes, aes(x = name, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  labs(title = "Valores faltantes por muestra (T0 y T5)",
       x = "Muestras",
       y = "Nº de valores faltantes") +
  theme_minimal() 
library(POMA)
library(ggtext)
library(magrittr)
#descartar muestras sin datos
met_exp_filtrado <- met_exp[, -c(7,9,10,11,22,28,29,31,36,37,38,39)]
#met_exp_filtrado

met_T0 <- SummarizedExperiment(
  assays = list(T0 = assay(met_exp_filtrado, "T0")),
  rowData = rowData(met_exp_filtrado),
  colData = colData(met_exp_filtrado)
)
met_T5 <- SummarizedExperiment(
  assays = list(T5 = assay(met_exp_filtrado, "T5")),
  rowData = rowData(met_exp_filtrado),
  colData = colData(met_exp_filtrado)
)
#imputar valores faltantes con POMA, se usa el método k-vecinos mas cercanos, los 0 se consideran NA y se eliminan aquellos metabolitos con más de un 20% de valores faltantes.
imputed_T5 = met_T5 %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 20)
imputed_T0 = met_T0 %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 20)
#boxplot T0
PomaBoxplots(imputed_T5, x = "samples") +
  labs(title = "Distribución de intensidades (T5)")
#boxplot T5
PomaBoxplots(imputed_T0, x = "samples") +
  labs(title = "Distribución de intensidades por muestra (T0)")
#boxplot T0
PomaBoxplots(imputed_T5[1:10,], x = "features") +
  labs(title = "Distribución de intensidades de los primeros 10 metabolitos (T5)")
#boxplot T5
PomaBoxplots(imputed_T0[1:10,], x = "features") +
  labs(title = "Distribución de intensidades de los primeros 10 metabolitos (T0)")

##QUITAR NEGATIVOS
assay(imputed_T5)[assay(imputed_T5)<0]=NA
assay(imputed_T0)[assay(imputed_T0)<0]=NA
no0_T5 = imputed_T5 %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 50)
no0_T0 = imputed_T0 %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 50)

##NORMALIZAR
normalized_T0 <- no0_T5 %>% 
  PomaNorm(method = "log_pareto")

normalized_T5 <- no0_T0 %>% 
  PomaNorm(method = "log_pareto")
#boxplot T0
PomaBoxplots(normalized_T5, x = "samples") +
  labs(title = "Distribución de intensidades (T5)")
#boxplot T5
PomaBoxplots(normalized_T0, x = "samples") +
  labs(title = "Distribución de intensidades por muestra (T0)")
#density T0
PomaDensity(normalized_T5[1:10,], x = "features") +
  labs(title = "Densidad de intensidades de los primeros 10 metabolitos (T5)")
#density T5
PomaDensity(normalized_T0[1:10,], x = "features") +
  labs(title = "Densidad de intensidades de los primeros 10 metabolitos (T0)")
PCA_T0=PomaPCA(
  normalized_T0,
  outcome = "Group",
  center = TRUE,
  scale = TRUE,
  ncomp = 4,
  labels = FALSE,
  ellipse = TRUE,
  load_length = 1
)

PCA_T5=PomaPCA(
  normalized_T5,
  outcome = "Group",
  center = TRUE,
  scale = TRUE,
  ncomp = 4,
  labels = FALSE,
  ellipse = TRUE,
  load_length = 1
)
PCA_T0$factors_plot+
  labs(title = "PCA (T0)")
PCA_T0$loadings_plot+
  labs(title = "Loading plot (T0)")
PCA_T0$eigenvalues_plot+
  labs(title = "Scree Plot (T0)")


PCA_T5$factors_plot+
  labs(title = "PCA (T5)")
PCA_T5$loadings_plot+
  labs(title = "Loading plot (T5)")
PCA_T5$eigenvalues_plot+
  labs(title = "Scree Plot (T5)")

PomaHeatmap(normalized_T0,covs=c(1,3,4))+
  labs(title = "heatmap (T0)")

POMA_T0=PomaClust(
  normalized_T0,
  method = "euclidean",
  k = NA,
  k_max = 8,
  show_clusters = TRUE,
  labels = TRUE
)
POMA_T0$mds_plot
POMA_T0$optimal_clusters_plot

####GUARDAR EL OBJETO
save(normalized_T0, file = "OBJETOSE_T0.Rda")
####GUARDAR EL OBJETO
save(normalized_T5, file = "OBJETOSE_T5.Rda")