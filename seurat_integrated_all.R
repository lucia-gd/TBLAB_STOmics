

# Lucía García Delgado

# Analysis, visualization, and integration of spatial datasets with Seurat
# Anotación usando los 4 niveles de especificidad: D1 (cancer), D2, D3, D4 (all)


# 1. INSTALAR LIBRERÍAS
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(presto)
library(devtools)


# 2. ARGUMENTOS
args = commandArgs(trailingOnly=TRUE)

# Inicializamos las variables con valores predeterminados
sample <- ""
data_dir <- ""
n.markers <- 30
max.dimensions <- 20
clustering.resolution <- 0.5

# Comprobamos si se pasaron los 5 argumentos
if(length(args) == 5){
  sample <- args[1]
  data_dir <- args[2]        
  n.markers <- as.numeric(args[3])  
  max.dimensions <- as.numeric(args[4])  
  clustering.resolution <- as.numeric(args[5])  
} 

# Imprimir los valores de las variables
print(paste("Sample: ", sample))
print(paste("Data path: ", data_dir))
print(paste("Number of Markers: ", n.markers))
print(paste("Max Dimensions: ", max.dimensions))
print(paste("Clustering Resolution: ", clustering.resolution))

list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
data.dir = data_dir


# 3. CARGA DE DATOS DE LA MUESTRA
visium_data <- Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  bin.size = NULL,
  filter.matrix = TRUE)
head(visium_data@meta.data)

# 4. PREPROCESAMIENTO
# Normalizar los datos
visium_data <- NormalizeData(visium_data)

# Buscar las características variables (genes)
visium_data <- FindVariableFeatures(visium_data)

# Escalar los datos
visium_data <- ScaleData(visium_data)

# Realizar PCA
visium_data <- RunPCA(visium_data, features = VariableFeatures(visium_data))

# Realizar UMAP
visium_data <- RunUMAP(visium_data, dims = 1:max.dimensions)

# Clustering 
visium_data <- FindNeighbors(visium_data, dims = 1:max.dimensions)
visium_data <- FindClusters(visium_data, resolution = clustering.resolution)

# Calcular % varianza explicada:
pca_var <- visium_data[["pca"]]@stdev^2
pca_var_percent <- pca_var / sum(pca_var) * 100
pca_var_percent  # Este vector te da el % de varianza por cada dimensión

varianza_tabla <- data.frame(
  PCs = c("PC1-10", "PC1-20", "PC1-30"),
  VarianzaExplicada = c(
    sum(pca_var_percent[1:10]),
    sum(pca_var_percent[1:20]),
    sum(pca_var_percent[1:30])
  )
)


##############################################################################
# 5. MARCADORES DIFERENCIALES
markers <- FindAllMarkers(visium_data)
markers$gene <- toupper(markers$gene)

cluster_markers <- markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  top_n(n.markers, avg_log2FC)


# 6. FUNCION PARA ASIGNAR TIPOS CELULARES
assign_cell_type <- function(cluster_markers, dataset, dataset_name) {
  merged <- merge(cluster_markers, dataset, by = "gene")
  assigned_cell_types <- merged %>%
    group_by(cluster) %>%
    summarize(assigned_cell_type = names(sort(table(cell_type), decreasing = TRUE))[1], .groups = "drop") %>%
    mutate(dataset = dataset_name)
  return(assigned_cell_types)
}


# 7. CARGA DE DATASETS DE REFERENCIA ORDENADOS POR NIVELES DE ESPECIFICIDAD
D1 <- read.csv("/home/luciagd/twoblab/luciagd/seurat_analysis/D1.csv")
D2 <- read.csv("/home/luciagd/twoblab/luciagd/seurat_analysis/D2.csv")
D3 <- read.csv("/home/luciagd/twoblab/luciagd/seurat_analysis/D3.csv")
D4 <- read.csv("/home/luciagd/twoblab/luciagd/seurat_analysis/D4.csv")

# Dataset cell_type - symbol
cell_type_dataset <- function(D) {
  data.frame(
    gene = toupper(D$Symbol),
    cell_type = D$cell_name,
    stringsAsFactors = FALSE
  ) %>% filter(gene != "" & !is.na(gene))
}

dataset_D1 <- cell_type_dataset(D1)
dataset_D2 <- cell_type_dataset(D2)
dataset_D3 <- cell_type_dataset(D3)
dataset_D4 <- cell_type_dataset(D4)


# 8. ASIGNACION INICIAL: inicialización del data.frame "cell_type_assignment"
todos_clusters <- unique(cluster_markers$cluster)

cell_type_assignment <- data.frame(
  cluster = todos_clusters,
  assigned_cell_type = "Unknown",
  dataset = NA,
  stringsAsFactors = FALSE
)

# 9. ANOTACIÓN POR ORDEN DE ESPECIFICIDAD
# Asignación con D1
res_D1 <- assign_cell_type(cluster_markers, dataset_D1, "D1") #anota: 0,2,4,5,6
cell_type_assignment$assigned_cell_type[match(res_D1$cluster, cell_type_assignment$cluster)] <- res_D1$assigned_cell_type
cell_type_assignment$dataset[match(res_D1$cluster, cell_type_assignment$cluster)] <- res_D1$dataset

# Asignación con D2
faltan <- cell_type_assignment$cluster[cell_type_assignment$assigned_cell_type == "Unknown"]
if (length(faltan) > 0) {
  res_D2 <- assign_cell_type(cluster_markers %>% filter(cluster %in% faltan), dataset_D2, "D2")
  cell_type_assignment$assigned_cell_type[match(res_D2$cluster, cell_type_assignment$cluster)] <- res_D2$assigned_cell_type
  cell_type_assignment$dataset[match(res_D2$cluster, cell_type_assignment$cluster)] <- res_D2$dataset
}

# Asignación con D3
faltan <- cell_type_assignment$cluster[cell_type_assignment$assigned_cell_type == "Unknown"]
if (length(faltan) > 0) {
  res_D3 <- assign_cell_type(cluster_markers %>% filter(cluster %in% faltan), dataset_D3, "D3")
  cell_type_assignment$assigned_cell_type[match(res_D3$cluster, cell_type_assignment$cluster)] <- res_D3$assigned_cell_type
  cell_type_assignment$dataset[match(res_D3$cluster, cell_type_assignment$cluster)] <- res_D3$dataset
}

# Asignación con D4
faltan <- cell_type_assignment$cluster[cell_type_assignment$assigned_cell_type == "Unknown"]
if (length(faltan) > 0) {
  res_D4 <- assign_cell_type(cluster_markers %>% filter(cluster %in% faltan), dataset_D4, "D4")
  cell_type_assignment$assigned_cell_type[match(res_D4$cluster, cell_type_assignment$cluster)] <- res_D4$assigned_cell_type
  cell_type_assignment$dataset[match(res_D4$cluster, cell_type_assignment$cluster)] <- res_D4$dataset
}

# ASIGNACIÓN DE LOS TIPOS CELULARES AL OBJETO SEURAT
cell_type_vector <- setNames(cell_type_assignment$assigned_cell_type, cell_type_assignment$cluster)

visium_data$cell_type <- sapply(visium_data$seurat_clusters, function(x) {
  key <- as.character(x)
  if (key %in% names(cell_type_vector)) {
    cell_type_vector[[key]]
  } else {
    "Unknown"
  }
})


#################################################################################################3

# 9. OBTENCIÓN DE ARCHIVOS INTERMEDIOS

# Filtrado: x genes más importantes por cluster:
top_per_cluster <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n.markers, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

top_per_cluster <- top_per_cluster %>%
  left_join(cell_type_assignment, by = "cluster")

# Tabla con el tipo celular de cada cluster + su tissue_type, cancer_type y tissue_class
cell_type_markers_extended <- function(D) {
  data.frame(
    gene = toupper(D$Symbol),
    cell_type = D$cell_name,
    tissue_type = D$tissue_type,
    cancer_type = D$cancer_type,
    tissue_class = D$tissue_class,
    stringsAsFactors = FALSE
  ) %>%
    filter(gene != "" & !is.na(gene))
}

dataset_D1_ext <- cell_type_markers_extended(D1)
dataset_D2_ext <- cell_type_markers_extended(D2)
dataset_D3_ext <- cell_type_markers_extended(D3)
dataset_D4_ext <- cell_type_markers_extended(D4) 

# FUnción de asignación
assign_cell_type_extended <- function(cluster_markers, dataset_ext, dataset_name) {
  merged <- merge(cluster_markers, dataset_ext, by = "gene")
  if (nrow(merged) == 0) {
    return(data.frame())  # Si no hay coincidencias
  }
  assigned_info <- merged %>%
    group_by(cluster) %>%
    summarize(
      assigned_cell_type = names(sort(table(cell_type), decreasing = TRUE))[1],
      assigned_tissue_type = names(sort(table(tissue_type), decreasing = TRUE))[1],
      assigned_cancer_type = names(sort(table(cancer_type), decreasing = TRUE))[1],
      assigned_tissue_class = names(sort(table(tissue_class), decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    mutate(dataset = dataset_name)
  
  return(assigned_info)
}

# Inicialización del data.frame extendido
cell_type_assignment_extended <- data.frame(
  cluster = sort(unique(cluster_markers$cluster)),
  assigned_cell_type = "Unknown",
  assigned_tissue_type = NA,
  assigned_cancer_type = NA,
  assigned_tissue_class = NA,
  dataset = NA,
  stringsAsFactors = FALSE
)


# Asignar con D1
res_D1 <- assign_cell_type_extended(cluster_markers, dataset_D1_ext, "D1")
cell_type_assignment_extended$assigned_cell_type[match(res_D1$cluster, cell_type_assignment_extended$cluster)] <- res_D1$assigned_cell_type
cell_type_assignment_extended$assigned_tissue_type[match(res_D1$cluster, cell_type_assignment_extended$cluster)] <- res_D1$assigned_tissue_type
cell_type_assignment_extended$assigned_cancer_type[match(res_D1$cluster, cell_type_assignment_extended$cluster)] <- res_D1$assigned_cancer_type
cell_type_assignment_extended$assigned_tissue_class[match(res_D1$cluster, cell_type_assignment_extended$cluster)] <- res_D1$assigned_tissue_class
cell_type_assignment_extended$dataset[match(res_D1$cluster, cell_type_assignment_extended$cluster)] <- res_D1$dataset

# Asignar con D2
faltan <- cell_type_assignment_extended$cluster[cell_type_assignment_extended$assigned_cell_type == "Unknown"]
if (length(faltan) > 0) {
  res_D2 <- assign_cell_type_extended(cluster_markers %>% filter(cluster %in% faltan), dataset_D2_ext, "D2")
  cell_type_assignment_extended$assigned_cell_type[match(res_D2$cluster, cell_type_assignment_extended$cluster)] <- res_D2$assigned_cell_type
  cell_type_assignment_extended$assigned_tissue_type[match(res_D2$cluster, cell_type_assignment_extended$cluster)] <- res_D2$assigned_tissue_type
  cell_type_assignment_extended$assigned_cancer_type[match(res_D2$cluster, cell_type_assignment_extended$cluster)] <- res_D2$assigned_cancer_type
  cell_type_assignment_extended$assigned_tissue_class[match(res_D2$cluster, cell_type_assignment_extended$cluster)] <- res_D2$assigned_tissue_class
  cell_type_assignment_extended$dataset[match(res_D2$cluster, cell_type_assignment_extended$cluster)] <- res_D2$dataset
}

# Asignar con D3
faltan <- cell_type_assignment_extended$cluster[cell_type_assignment_extended$assigned_cell_type == "Unknown"]
if (length(faltan) > 0) {
  res_D3 <- assign_cell_type_extended(cluster_markers %>% filter(cluster %in% faltan), dataset_D3_ext, "D3")
  cell_type_assignment_extended$assigned_cell_type[match(res_D3$cluster, cell_type_assignment_extended$cluster)] <- res_D3$assigned_cell_type
  cell_type_assignment_extended$assigned_tissue_type[match(res_D3$cluster, cell_type_assignment_extended$cluster)] <- res_D3$assigned_tissue_type
  cell_type_assignment_extended$assigned_cancer_type[match(res_D3$cluster, cell_type_assignment_extended$cluster)] <- res_D3$assigned_cancer_type
  cell_type_assignment_extended$assigned_tissue_class[match(res_D3$cluster, cell_type_assignment_extended$cluster)] <- res_D3$assigned_tissue_class
  cell_type_assignment_extended$dataset[match(res_D3$cluster, cell_type_assignment_extended$cluster)] <- res_D3$dataset
}

# Asignar con D4
faltan <- cell_type_assignment_extended$cluster[cell_type_assignment_extended$assigned_cell_type == "Unknown"]
if (length(faltan) > 0) {
  res_D4 <- assign_cell_type_extended(cluster_markers %>% filter(cluster %in% faltan), dataset_D4_ext, "D4")
  cell_type_assignment_extended$assigned_cell_type[match(res_D4$cluster, cell_type_assignment_extended$cluster)] <- res_D4$assigned_cell_type
  cell_type_assignment_extended$assigned_tissue_type[match(res_D4$cluster, cell_type_assignment_extended$cluster)] <- res_D4$assigned_tissue_type
  cell_type_assignment_extended$assigned_cancer_type[match(res_D4$cluster, cell_type_assignment_extended$cluster)] <- res_D4$assigned_cancer_type
  cell_type_assignment_extended$assigned_tissue_class[match(res_D4$cluster, cell_type_assignment_extended$cluster)] <- res_D4$assigned_tissue_class
  cell_type_assignment_extended$dataset[match(res_D4$cluster, cell_type_assignment_extended$cluster)] <- res_D4$dataset
}


# Cuantificación del número de spots por cluster:
# data.frame con número de spots por cluster
spot_counts <- as.data.frame(table(visium_data$seurat_clusters))
colnames(spot_counts) <- c("cluster", "n_spots")

# Añadimos la columna con el tipo celular asignado: cluster, n_spots, assigned_cell_type
table_nspots_cell_type <- merge(spot_counts, cell_type_assignment, by = "cluster")


################################################################################################

# 10. CARPETA DE SALIDA PARA GUARDAR LOS RESULTADOS
output_base <- "/home/luciagd/twoblab/luciagd/seurat_analysis/results_integrated"
output_dir <- file.path(output_base, sample)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 11. GUARDADO DE GRÁFICAS

# CLustering plots with no annotation
p1 <- DimPlot(visium_data, reduction = "umap", label = TRUE, label.box=TRUE)
p2 <- SpatialDimPlot(visium_data, label = TRUE, label.size = 4, pt.size.factor = 5.5) + labs(fill = "") +
  theme(legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(override.aes = list(size = 4))  
  )
clustering_images <- p1 + p2
clustering_images <- clustering_images + ggtitle("UMAP plot and Spatial Plot - Clusters")
ggsave(filename = file.path(output_dir, "cluster_images.png"), plot = clustering_images, width = 12, height = 6, dpi = 300)


# UMAP Plot colored by cell types
umap_plot_annotated <- DimPlot(visium_data, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE, label.box=TRUE) +
  theme(legend.position = "none") +
  ggtitle("UMAP Plot -") # HE ELIMINADO EL TÍTULO
ggsave(filename = file.path(output_dir, "umap_plot_annotated.png"), plot = umap_plot_annotated, width = 6, height = 6, dpi = 300)

# Spatial Plot colored by cell types
spatial_plot_annotated <- SpatialDimPlot(visium_data, group.by = "cell_type", label = TRUE, repel = TRUE, label.box=TRUE, label.size = 6, pt.size.factor = 5.5) + 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 13)) + 
  labs(fill = "Tipo celular") + 
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  ggtitle(paste("Spatial Plot -", sample, 
                "\nn.markers:", n.markers, 
                " | max.dimensions:", max.dimensions, 
                " | clustering.resolution:", clustering.resolution))
ggsave(filename = file.path(output_dir, "spatial_plot_annotated.png"), plot = spatial_plot_annotated, width = 10, height = 8, dpi = 300)

# Combined plots
combined_plot <- umap_plot_annotated / spatial_plot_annotated  # Vertical layout
ggsave(filename = file.path(output_dir, "combined_plot.png"), plot = combined_plot, width = 10, height = 18, dpi = 300) 


anotated_plots <- umap_plot_annotated + spatial_plot_annotated
ggsave(filename = file.path(output_dir, "anotated_plots.png"), plot = combined_plot, width = 12, height = 6, dpi = 300) 

# 12. GUARDADO DE TABLAS
write.csv(top_per_cluster, file = file.path(output_dir, "top_markers_with_celltype.csv"), row.names = FALSE)
write.csv(cell_type_assignment, file = file.path(output_dir, "cell_type_assignment.csv"), row.names = FALSE)
write.csv(cell_type_assignment_extended, file = file.path(output_dir, "cell_type_assignment_extended.csv"), row.names = FALSE)
write.csv(table_nspots_cell_type, file = file.path(output_dir, "cluster_nspots_celltype.csv"), row.names = FALSE)
write.csv(spot_counts, file = file.path(output_dir, "table_spots_counts.csv"), row.names = FALSE)
write.csv(pca_var_percent, file = file.path(output_dir, "pca_var_percent.csv"), row.names = FALSE)
write.csv(varianza_tabla,file = file.path("varianza_explicada.csv"), row.names = FALSE)



# Ejecución desde la terminal:
# Rscript seurat_integrated_all.R sample data_dir n.markers max.dimensions clustering.resolution
# Rscript seurat_integrated_all.R A694Tumor '/home/luciagd/twoblab/luciagd/LUCIA_SpatialTranscriptomics/results/A694Tumor/outs' 10 10 0.5


