library(hdf5r)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load outputs of the spaceranger pipeline (this is an example for sample10) into Seurat using the Load10X_Spatial() function and returns a Seurat object that contains both the spot-level expression data along with the associated image of the tissue slice
sample10 <- Load10X_Spatial(
  "Sample10",
  filename = "Sample10_NSD_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "tissue_hires_image",
  filter.matrix = TRUE,
  to.upper = FALSE
)

# SCTransform normalizes the data, detects high-variance features, and stores the data in the SCT assay
sample10 <- SCTransform(sample10, assay = "Spatial", verbose = FALSE)
# Dimensionality reduction and clustering
sample10 <- RunPCA(sample10, assay = "SCT", verbose = FALSE)
sample10 <- FindNeighbors(sample10, reduction = "pca", dims = 1:30)
sample10 <- FindClusters(sample10, verbose = FALSE)
sample10 <- RunUMAP(sample10, reduction = "pca", dims = 1:30)

#Create seurat object from Allen Brain Atlas table of cell metadata and gene expression matrix (https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x)
# load metadata once
meta = fread("metadata.csv")
subclass = unique(meta$subclass_label)
simplify = c("sample_name", "full_genotype_label", "neighborhood_label", 
             "subclass_label", "class_label", "donor_sex_label", "region_label")
meta[, names(meta)[!(names(meta) %in% simplify)] := NULL]
print(paste("metadata size:", format(object.size(meta), units="auto")))

# load gene counts
counts = fread("matrix.csv")
gn = ncol(counts)

# for each matrix: read in, make Seurat and remove from memory
sc.list = lapply(subclass, function(i) {
  # inner join to ensure all cells have metadata
  dt = merge(counts, meta[subclass_label == i], by="sample_name")
  # transpose counts (g x sn), get rownames and meta.data from dt3
  x = CreateSeuratObject(counts=transpose(dt[, 1:gn], make.names="sample_name"), 
                         row.names=names(dt)[2:gn], 
                         project="AIBS_mouse_cort-hipp_10x",
                         meta.data=data.frame(dt[, names(meta), with=F], row.names=1))
  # do you wish to randomly downsample? (https://github.com/satijalab/seurat/issues/3108)
  x = x[, sample(colnames(x), size=min(20000, nrow(dt)), replace=F)]
  remove(dt)
  # update and exit
  print(paste0("I finished ", i))
  return(x)
})

# merge and save
print(paste("sc.list size:", format(object.size(sc.list), units="auto")))
# quick&dirty: takes huge amounts of memory (something about merge)
allen = Reduce(merge, sc.list) 
print(paste("allen size:", format(object.size(allen), units="auto")))
saveRDS(allen, file="allen.rds")
remove(sc.list)


# Integration with single-cell data
# Load the pre-process scRNA-seq reference
allen_reference <- readRDS("allen.rds")
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
+     RunPCA(verbose = FALSE) %>%
+     RunUMAP(dims = 1:30)

# Perform label transfer
anchorsS10 <- FindTransferAnchors(reference = allen_reference, query = sample10, normalization.method = "SCT")

# Output, for each spot, a probabilistic classification for each of the scRNA-seq derived classes
predictions.assayS10 <- TransferData(anchorset = anchorsS10, refdata = allen_reference$subclass, prediction.assay = TRUE, weight.reduction = sample10[["pca"]], dims = 1:30)

# Add these predictions as a new assay in the Seurat object
sample10[["predictions"]] <- predictions.assayS10
DefaultAssay(sample10) <- "predictions"

# Distinguish prediction score between distinct sequential layers of the neocortex
SpatialFeaturePlot(sample10, features = c("L2/3 IT", "L4", "L5 IT", "L6 CT", "CA1", "DG"), pt.size.factor = 1.6, crop = TRUE, alpha = c(0.1, 1))
