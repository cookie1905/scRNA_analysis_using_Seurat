
library(dplyr) # 5.3.0
library(Seurat) # 1.1.4
library(patchwork) # 1.3.2
library(scales)


setwd("C:\\scRNA_analysis_using_Seurat")

donor1_counts <- Read10X_h5("data/Donor1_filtered_feature_bc_matrix.h5")
donor2_counts <- Read10X_h5("data/Donor2_filtered_feature_bc_matrix.h5")
donor3_counts <- Read10X_h5("data/Donor3_filtered_feature_bc_matrix.h5")

donor1 <- CreateSeuratObject(
  counts = donor1_counts,
  min.cells = 3,   
  min.features = 200   
)
donor1


Cells(donor1) # get cell barcodes
Features(donor1) # get feature names
donor1[["RNA"]]$counts # Get counts/expression matrix


donor2 <- CreateSeuratObject(
  counts = donor2_counts,
  min.cells = 3,   
  min.features = 200   
)

donor3 <- CreateSeuratObject(
  counts = donor3_counts,
  min.cells = 3,   
  min.features = 200   
)

donor1$sample <- "donor1"
donor1[[]] 

donor2$sample <- "donor2"
donor3$sample <- "donor3"

combined <- merge(donor1, y = list(donor2, donor3), add.cell.ids = c("donor1", "donor2", "donor3"))
dim(combined)
head(Cells(combined))
head(combined[[]])

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
head(combined[[]])

VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
dim(combined)

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
head(combined[["RNA"]]$data) # Normalized counts is in layer data

combined <-  FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined


top10 <- head(VariableFeatures(combined), 10)
top10


# all.genes <- rownames(combined)
# combined <- ScaleData(combined, features = all.genes)

combined <- ScaleData(combined, vars.to.regress = "percent.mt")

combined <- RunPCA(combined)
ElbowPlot(combined)

combined <- FindNeighbors(combined, dims = 1:10)
combined <- FindClusters(combined, resolution = 0.25)

combined <- RunUMAP(combined, dims = 1:10) 
DimPlot(combined, reduction = "umap", label = TRUE)

saveRDS(combined, file = "data/combined.rds")
