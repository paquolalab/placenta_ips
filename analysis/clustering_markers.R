

### Here we perform Clustering Analysis and find Markers of each cluster by using Seurat V.3.0.1 R package 
### Download the R object with counts matrix here: '/Rdata/counts.rds'

library(Seurat)
options(stringsAsFactors = FALSE)

## Read Counts Matrix
counts <- readRDS("counts.rds")

## Filtering and Normalization
# Here we want to create a Seurat object for the counts matrix
# We'll filter out cells with less than 1000 genes, and remove genes expressed in less than 3 cells
sc <- CreateSeuratObject(counts, min.features = 1000, min.cells = 3)

# Here we want to calculate the percentage of MT RNA for each cell, and remove cells with more than 20% of MT RNA.
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
sc <- subset(sc, subset = percent.mt < 20 & nFeature_RNA > 1000)

# Perform SCT normalization, regressing for library size and MT percentage
sc <- SCTransform(sc, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)

# Run PCA and UMAP (based on 30 first PCs) 
sc <- RunPCA(sc, verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:30, verbose = FALSE)

# Create a vector with the colors we want to use for each sample in the UMAP plot
colors <- levels(sc)
names(colors) <- c("#4292C6", "#08306B", "#9E9AC8", "#6A51A3","#3F007D","#006D2C", "#F16913", "#D94801", "#7F2704", "#CB181D") 

# Plot the UMAP to observe the cell diversity over the samples and days

DimPlot(sc,
        reduction = "umap",
        cols = names(colors),
        label = TRUE,
        pt.size = 0.3,
        label.size = 0.1)

## Clustering
# Find clusters of cells
sc <- FindNeighbors(sc, dims = 1:30, verbose = FALSE)
sc <- FindClusters(sc, verbose = FALSE)

# Plot the UMAP to observe the clusters
DimPlot(sc,
        reduction = "umap",
        label = TRUE,
        pt.size = 0.3)
dev.off()


# Rename clusters - from cluster 0 to cluster 18
new.cluster.ids <- c("BR_4", "TS_2", "BR_3", "SC", "iPS_2","TS_3", "BR", "BR_2", "TS_1", "iPS_1",
                     "TS_5", "SC_3", "TS_6", "BR_1", "TS","SC_1", "iPS_3", "SC_2", "TS_4")
names(new.cluster.ids) <- levels(sc)
sc <- RenameIdents(sc, new.cluster.ids)
sc <- AddMetaData(sc, Idents(sc), col.name="new_cluster_ids")

# Plot UMAP with renamed clusters 
cluster.colors <- c("#54278F", "#A63603", "#807DBA", "#08306B", "#238B45", "#7F2704", "#3F007D", "#6A51A3", "#FD8D3C", "#00441B",
                    "#F16913", "#08519C", "#A50F15", "#9E9AC8","#67000D","#4292C6", "#74C476", "#2171B5", "#D94801")
DimPlot(sc,
        reduction = "umap",
        cols = cluster.colors, 
        label = TRUE,
        pt.size = 0.3,
        label.size = 0.1)


# Build a table with numbers of cells per cluster and sample
cluster.sample <- lapply(new.cluster.ids, function(x){ return(sc@meta.data[sc@meta.data$new_cluster_ids==x,]) })
cluster.sample <- sapply(cluster.sample, function(x) table(x$orig.ident))
colnames(cluster.sample) <- levels(sc)
write.table(cluster.sample,"cells_sample_cluster.tsv", sep='\t')


## Find Markers of each cluster
# The expression of each cluster is tested against the union of all remaining clusters
sc.markers <- FindAllMarkers(sc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")

sc.markers <- sc.markers[sc.markers$p_val_adj < 0.01,]
# Save table with markers
write.table(sc.fmarkers, "markers_clusters.tsv", sep='\t', row.names=FALSE)

## Plots 







