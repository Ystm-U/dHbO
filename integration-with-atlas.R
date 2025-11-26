library(Seurat)
library(tidyverse)
library(Matrix)
options(future.globals.maxSize = 16000*1024^2)

setwd(" ... /dHBO_project/merge/")

d30 <- readRDS("../ ... /d30-merged.rds")
hdev <- readRDS("./Data/human_dev_subset.rds") # 100k cells from Braun et al., 2023, Science 

#### processing for developing human brain atlas ####
# extract count matrix and metadata
counts <- GetAssayData(hdev, assay = "RNA", slot = "counts")
meta_data <- hdev@meta.data
feature_data <- hdev[["RNA"]]@meta.features

# map Ensembl IDs to gene symbols
ensembl_ids <- rownames(counts)
gene_symbols <- feature_data$Gene[match(ensembl_ids, feature_data$Accession)]

# replace missing gene names with Ensembl IDs
gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- ensembl_ids[is.na(gene_symbols) | gene_symbols == ""]

# confirm uniqueness
rownames(counts) <- make.unique(gene_symbols)

# create a new Seurat object
hdev <- CreateSeuratObject(counts = counts)

# restore metadata
meta_data <- meta_data[colnames(hdev), , drop = F]
hdev <- AddMetaData(hdev, metadata = meta_data)

#### preprocessing and integration ####
d30[["SCT"]]@misc <- list()
seurat_list <- list(hdev, d30)

seurat_list <- lapply(seurat_list, function(x) {
  SCTransform(x, verbose = F)
})

features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

merged <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

#### dimention reduction ####
merged <- RunPCA(merged, verbose = F)
merged <- RunUMAP(merged, dims = 1:30,
                  reduction="integrated.cca", reduction.name = "umap_cca")

#### batch name correlction ###
table(merged$Subregion, useNA = "ifany")  # before replacement
merged$Subregion[is.na(merged$Subregion)] <- "X_dHbO"
table(merged$Subregion, useNA = "ifany")  # after replacement

table(merged$batch, useNA = "ifany")  # before replacement
merged$batch[is.na(merged$batch)] <- "human_dev"
table(merged$batch, useNA = "ifany")  # after replacement

table(merged$seurat_clusters, useNA = "ifany")  # before replacement
merged$seurat_clusters[is.na(merged$seurat_clusters)] <- "human_dev"
table(merged$seurat_clusters, useNA = "ifany")  # after replacement

#### UMAP plots ###
# clusters from brain atlas
clusters <- levels(factor(merged$Subregion))
highlight <- c("Midbrain","Midbrain dorsal","Midbrain ventral","Medulla","Pons","Hindbrain")
order <- c(setdiff(clusters, highlight), highlight)

p <- DimPlot(merged, group.by = "Subregion",
             cols = c(rep("#cccccc",13), 
                      rep("#99dddd",3), "#6666ff", "#ff7777","#339933"),
             order = rev(order),
             raster = F
) + NoAxes() +
  NoLegend() +
  theme(
    plot.title = element_blank(),
    panel.background = element_blank(),
    plot.margin=unit(c(0,0,0,0), "cm")
  )
ggsave(paste0("hdev_merged_dimplot_Subregion.tiff"), plot = p, device = "tiff",
       path = "Output/", width = 6, height = 6, units = "in", dpi = 600)

# clusters from d30 dHbOs
clusters <- levels(factor(merged@active.ident))
highlight <- rev(c("1", "5", "3", "0", "6", "7"))
order <- c(setdiff(clusters, highlight), highlight)

p <- DimPlot(merged, #group.by = "active.ident",
             cols = c(rep("#cccccc",66), 
                      rev(c("#329767","#67c0a2","#f07047","#faad61","#4293c3","#2168ac"))),
             order = rev(order),
             raster = F
) + NoAxes() +
  NoLegend() +
  theme(
    plot.title = element_blank(),
    panel.background = element_blank(),
    plot.margin=unit(c(0,0,0,0), "cm")
  )
ggsave(paste0("hdev_merged_dimplot_clusters.tiff"), plot = p, device = "tiff",
       path = "Output/", width = 6, height = 6, units = "in", dpi = 600)

# Feature plots
genes <- c("SOX2", "LHX9", "HOXB2", "CRABP1", "MAPT", "ELAVL3")
for (i in genes){
FeaturePlot(merged, i, pt.size = 0.1, max.cutoff = 2,
            cols = c("lightgrey", "#001F8F"),
            order = T, raster = F) + NoLegend() + NoAxes() +
    theme(
      plot.title = element_text(size = 28, face = "italic"),
      panel.background = element_blank(),
      plot.margin=unit(c(0,0,0,0), "cm")
    ) -> p
ggsave(paste0(i, ".tiff"), plot = p, device = "tiff",
       path = "Output/", width = 4, height = 4, units = "in", dpi = 600)
}

saveRDS(merged, "./Data/merged.rds")
