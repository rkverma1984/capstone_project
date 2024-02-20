library(Matrix) #Matrix is a package that provides classes for dense and sparse matrices and operations on them using BLAS and LAPACK libraries.
library(Signac) #Signac is a package for the analysis of single-cell chromatin data.
library(Seurat) #Seurat is a package designed for QC, analysis, and exploration of single-cell RNA-seq data.
library(SeuratDisk) #SeuratDisk is a package that provides a disk-backed storage backend for Seurat objects.
library(SeuratObject) #SeuratObject is a package that provides a class for storing
library(tidyverse) #tidyverse is a collection of R packages designed for data science.
library(janitor) #janitor is a package that provides a clean API for cleaning data.
library(conflicted) #conflicted is a package that provides a conflict resolution system for R.
library(scDblFinder) #scDblFinder is a package for identifying doublets in single-cell RNA-seq data.
library(ggplot2) #ggplot2 is a system for declaratively creating graphics, based on The Grammar of Graphics.
library(EnsDb.Hsapiens.v86) #EnsDb.Hsapiens.v86 is a Bioconductor package that provides the full genome annotation
library(BSgenome.Hsapiens.UCSC.hg38) #BSgenome.Hsapiens.UCSC.hg38 is a Bioconductor package that provides the full genome
library(BiocManager) #BiocManager is a package that provides a convenient way to install and update Bioconductor packages.
library(BiocGenerics) #BiocGenerics provides S4 generic functions used in Bioconductor.
library(ggforce) #ggforce is a package that extends ggplot2 to provide a wider range of functionality for visualizing data.
library(TFBSTools) #TFBSTools is an R package for the analysis of transcription factor binding site (TFBS) sequences.
library(JASPAR2020) #JASPAR is the largest open-access database of curated and non-redundant transcription factor (TF) binding profiles from six different taxonomic groups.
library(ArchR) #ArchR is a package for the analysis of single-cell chromatin data.
library(ggseqlogo) #ggseqlogo is a package that provides a 'ggplot2' extension for drawing publication-ready sequence logos.

## Use dplyr::filter instead of stats::filter
conflict_prefer("filter", "dplyr")    # dplyr::filter() over stats::filter()


#### ----------------------------------------------------------------------------------------

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) #GetGRangesFromEnsDb is a function that extracts genomic ranges from an EnsDb object.
seqlevels(annotations) <- paste0('chr', seqlevels(annotations)) #seqlevels is a function that returns the sequence levels of a GRanges object.
seqlevelsStyle(annotations) <- 'UCSC' #seqlevelsStyle is a function that sets the style of the sequence levels of a GRanges object.
genome(annotations) <- "hg38" #genome is a function that sets the genome of a GRanges object.


#### ----------------------------------------------------------------------------------------

#### set directory paths #####
path_to_base_dir <- 'base path' #modify the base path
base_dir <- paste0(path_to_base_dir, "/capstone_project")
w_dir <- file.path(base_dir, "scripts")
data_dir <- file.path(base_dir, "data")


#### change current directory to working directory ####
{setwd(w_dir)} #setwd is a function that sets the current working directory to the specified path.
print(paste("WORKING DIR", w_dir))


#### ----------------------------------------------------------------------------------------

#### read data ####

counts <- Read10X_h5(filename = file.path(data_dir, "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5"))


#Create a seurat object containing the RNA data and metadata
lymph <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
) #CreateSeuratObject is a function that creates a Seurat object from a matrix of gene expression counts.


lymph #print the seurat object


head(lymph@meta.data) #print the metadata of the seurat object


print(paste("lymph dimensions", dim(lymph))) #print the dimensions of the seurat object


atac_fragment_file <-  file.path(data_dir, "lymph_node_lymphoma_14k_atac_fragments.tsv.gz") #file.path is a function that creates a path to a file or directory.

atac_fragment_file #print the path to the atac fragment file


# Add the ATAC data as another assay in the Seurat object. CreateChromatinAssay is a function that creates a chromatin assay from a matrix of peak counts.
chrom_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = atac_fragment_file,
  annotation = annotations,
)


lymph[["ATAC"]] <- chrom_assay #add the ATAC data as another assay in the Seurat object

head(lymph)

rm(chrom_assay, counts) #remove the chromatin assay and counts


#### ----------------------------------------------------------------------------------------

#### Perform QC ####

DefaultAssay(lymph) <- "ATAC" #DefaultAssay is a function that sets the default assay for a Seurat object.

lymph <- NucleosomeSignal(lymph) #NucleosomeSignal is a function that computes the nucleosome signal for a Seurat object.

head(lymph)

lymph <- TSSEnrichment(lymph, fast=TRUE) #TSSEnrichment is a function that computes the TSS enrichment for a Seurat object.

head(lymph)

options( scipen = 999 ) #options is a function that sets or retrieves a number of global options for the current R session.
p <- VlnPlot(
  object = lymph,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0,
  adjust = 1
)

svg(filename = "../results/violin_plots_rna_atac.svg",  width = 12, height = 5)
p


saveRDS(lymph, file="../results/lymph_final.rds")
fs::file_size("../results/lymph_final.rds")


#### filter out low quality cells using cutoffs for RNA and ATAC data and TSS enrichment
filtered_lymph <- subset(
  x = lymph,
  subset =nCount_RNA < 10000 &
    
    nCount_ATAC < 60000 &
    
    TSS.enrichment < 7 &
    TSS.enrichment > 2 &
    
    nucleosome_signal < 1.25 &
    nucleosome_signal > 0.4 
)

print(lymph)

print(filtered_lymph)


#### plot filtered data ####

p2 <- VlnPlot(
  object = filtered_lymph,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0,
)

svg(filename = "../results/violin_plots_rna_atac_filtered.svg",  width = 12, height = 5)
p2


#### ----------------------------------------------------------------------------------------

# call peaks using macs3. CallPeaks is a function that calls peaks using MACS3. MACS3 is a software for identifying transcription factor binding sites and histone modification sites from ChIP-seq data. Note that the path to the MACS2 executable must be specified.
peaks <- CallPeaks(filtered_lymph, macs2.path = '/home/ravikumarv/mambaforge/envs/py38/bin/macs3')

saveRDS(peaks, file="../results/peaks_final.rds")
fs::file_size("../results/peaks_final.rds")

rm(lymph) # remove lymph object

head(peaks)

#### remove peaks on nonstandard chromosomes and in genomic blacklist regions. 
#### keepStandardChromosomes is a function that removes peaks on nonstandard chromosomes from a GRanges object. 
#### subsetByOverlaps is a function that subsets a GRanges object based on overlaps with another GRanges object.

blacklist_hg38_unified <- ArchR:::.getBlacklist("hg38") #getBlacklist is a function that retrieves the genomic blacklist for a specified genome.
peaks_filtered <- keepStandardChromosomes(peaks, pruning.mode = "coarse") #remove peaks on nonstandard chromosomes from a GRanges object.
peaks_filtered <- subsetByOverlaps(x = peaks_filtered, ranges = blacklist_hg38_unified, invert = TRUE) #remove not well defined (blacklisted) peaks from a GRanges object.

peaks_filtered

#### quantify counts in each peak
macs3_counts <- FeatureMatrix(
  fragments = Fragments(filtered_lymph),
  features = peaks_filtered,
  cells = colnames(filtered_lymph)
)

filtered_lymph

head(filtered_lymph)

# add the peaks to the Seurat object
filtered_lymph[["peaks"]] <- CreateChromatinAssay(
  counts = macs3_counts,
  fragments = atac_fragment_file,
  annotation = annotations
)

head(filtered_lymph)

filtered_lymph

saveRDS(filtered_lymph, file="../results/filtered_lymph_final.rds")


saveRDS(filtered_lymph, file="../results/filtered_lymph_after_peaks_final.rds")
fs::file_size("../results/filtered_lymph_after_peaks_final.rds")

filtered_lymph <-readRDS("../results/filtered_lymph_after_peaks_final.rds")

#### run TSS enrichment on peaks and save the object
filtered_lymph <- TSSEnrichment(filtered_lymph, fast = FALSE, assay='peaks')

saveRDS(filtered_lymph, file="../results/filtered_lymph_after_peaks_TSS_final.rds")
fs::file_size("../results/filtered_lymph_after_peaks_TSS_final.rds")


#### ----------------------------------------------------------------------------------------


filtered_lymph <-readRDS("../results/filtered_lymph_after_peaks_TSS_final.rds")


rm(filtered_lymph_sctransform)

######## plot DensityScatter for filtered data ########
#### plot nCount_RNA vs TSS.enrichment ####
d1 <- DensityScatter(filtered_lymph, x = 'nCount_RNA', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
svg(filename = "../results/nCount_RNA_vs_TSS_enrichment.svg",  width = 8, height = 6)
d1


#### plot nCount_peaks vs TSS.enrichment ####
d2 <- DensityScatter(filtered_lymph, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
svg(filename = "../results/nCount_peaks_vs_TSS_enrichment.svg",  width = 8, height = 6)
d2

filtered_lymph$high.tss <- ifelse(filtered_lymph$TSS.enrichment > 2.5, 'High', 'Low') #create a new column in the metadata of the Seurat object

DefaultAssay(filtered_lymph) <- "peaks"
t1 <- TSSPlot(filtered_lymph, group.by = 'high.tss') + NoLegend() # TSSPlot is a function that plots the TSS enrichment for a Seurat object.
svg(filename = "../results/TSS_plot.svg",  width = 10, height = 5)
t1

DefaultAssay(filtered_lymph) <- "peaks"
filtered_lymph$nucleosome_group <- ifelse(filtered_lymph$nucleosome_signal > 0.6, 'NS > 0.6', 'NS < 0.6') # create a new column in the metadata of the Seurat object
t2 <- FragmentHistogram(object = filtered_lymph, group.by = 'nucleosome_group') # plot a histogram of the fragment length distribution for a Seurat object.
svg(filename = "../results/fragment_histogram.svg",  width = 10, height = 5)
t2

##### plot violin plots for filtered data #####
p3 <- VlnPlot(
  object = filtered_lymph,
  features = c("nCount_RNA", "TSS.enrichment", "nucleosome_signal", "nCount_peaks"),
  ncol = 4,
  pt.size = 0,
)+ NoLegend()

svg(filename = "../results/violin_plots_rna_peaks.svg",  width = 12, height = 5)
p3


#### additional filtering of data to filter out low quality cells
filtered_lymph <- subset(
  x = filtered_lymph,
  subset = nCount_peaks < 20000 &
     nCount_peaks > 408 &
     
    nCount_RNA < 10000 &
    nCount_RNA > 330 &
    
#    nCount_ATAC < 50000 &
#    nCount_ATAC > 783 & 
#    
#    TSS.enrichment < 7 &
#    TSS.enrichment > 2 &
#    
#    nucleosome_signal < 1.25 &
#    nucleosome_signal > 0.4    
)


#### plot violin plots for filtered data ####
p3 <- VlnPlot(
  object = filtered_lymph,
  features = c("nCount_RNA", "TSS.enrichment", "nucleosome_signal", "nCount_peaks"),
  ncol = 4,
  pt.size = 0,
)+ NoLegend()

svg(filename = "../results/violin_plots_rna_peaks_filtered.svg",  width = 12, height = 5)
p3


#### plot DensityScatter for filtered data ####
d3 <- DensityScatter(filtered_lymph, x = 'nCount_RNA', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
svg(filename = "../results/nCount_RNA_vs_TSS_enrichment_filtered.svg",  width = 8, height = 6)
d3

d4 <- DensityScatter(filtered_lymph, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
svg(filename = "../results/nCount_peaks_vs_TSS_enrichment_filtered.svg",  width = 8, height = 6)
d4

saveRDS(filtered_lymph, file="../results/filtered_lymph_after_cleaned_peaks_final.rds") # save the filtered Seurat object
fs::file_size("../results/filtered_lymph_after_cleaned_peaks_final.rds") 


#### ----------------------------------------------------------------------------------------

filtered_lymph <-readRDS("../results/filtered_lymph_after_cleaned_peaks_final.rds") # read the filtered Seurat object


#### scale and run PCA on RNA data
DefaultAssay(filtered_lymph) <- "RNA"
filtered_lymph_sctransform <- SCTransform(filtered_lymph) # SCTransform performs variance stabilizing transformation and regresses out unwanted sources of variation for a Seurat object.
filtered_lymph_sctransform <- RunPCA(filtered_lymph_sctransform) # RunPCA is a function that runs PCA on a Seurat object.

filtered_lymph_sctransform <- RunUMAP(filtered_lymph_sctransform, reduction='pca', dims = 1:50, reduction.name = 'umap.rna_50', reduction.key = 'rnaUMAP_50_') # RunUMAP and append the results to the Seurat object

dim1 <- DimPlot(filtered_lymph_sctransform , reduction = 'umap.rna_50', label = FALSE, repel = TRUE, label.size = 3) + NoLegend() # plot dimensionality reduction data for RNA.
svg(filename = "../results/umap_rna_50.svg",  width = 8, height = 6)
dim1

#### find top features
DefaultAssay(filtered_lymph) <- "peaks"
filtered_lymph_sctransform <- FindTopFeatures(filtered_lymph_sctransform, min.cutoff = 5)
cat('\n')
print(filtered_lymph_sctransform)


filtered_lymph_sctransform <- RunTFIDF(filtered_lymph_sctransform) # compute the TF-IDF transformation for a Seurat object.
filtered_lymph_sctransform <- RunSVD(filtered_lymph_sctransform) # dimention reduction using SVD.
filtered_lymph_sctransform <- RunUMAP(filtered_lymph_sctransform, reduction='lsi', dims = 2:50, reduction.name = 'umap.atac_50', reduction.key = 'atacUMAP_50_') 

dim2 <- DimPlot(filtered_lymph_sctransform , reduction = 'umap.atac_50', label = FALSE, repel = TRUE, label.size = 3) + NoLegend() # plot dimensionality reduction data for ATAC.
svg(filename = "../results/umap_atac_50.svg",  width = 8, height = 6)
dim2


#### ----------------------------------------------------------------------------------------

#### build a joint neighbor graph using both assays ####
filtered_lymph_sctransform <- FindMultiModalNeighbors(
  object = filtered_lymph_sctransform,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  verbose = TRUE
) 

filtered_lymph_sctransform <- RunUMAP(
  object = filtered_lymph_sctransform,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE,
  reduction.name = "wnn.umap_50", 
  reduction.key = "wnnUMAP_50_", 
) # RunUMAP and append the results to the Seurat object

dim3 <- DimPlot(filtered_lymph_sctransform , reduction = 'wnn.umap_50', label = FALSE, repel = TRUE, label.size = 3) + NoLegend() # plot dimensionality reduction data for combined RNA and ATAC.
svg(filename = "../results/wnn_umap_50.svg",  width = 8, height = 6)
dim3


#### ----------------------------------------------------------------------------------------

#### Add cell annotations to the Seurat object

ref_file <- file.path(data_dir, "pbmc_multimodal.h5seurat")
reference <- LoadH5Seurat(ref_file, assays = list("SCT" = "counts"), reductions = 'spca') # load annotations

reference <- UpdateSeuratObject(reference) # UpdateSeuratObject is a function that updates a Seurat object to the current version of Seurat.

# transfer cell-type labels from reference to query. FindTransferAnchors is a function that finds anchors between two Seurat objects.
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = filtered_lymph_sctransform,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)


#### transfer data from a reference to a query using a set of anchors. ####
predictions <- TransferData(
  anchorset = transfer_anchors,
  reference = reference,
  refdata =  "celltype.l1",
  weight.reduction = filtered_lymph_sctransform[['pca']],
  dims = 1:50
)

#### AddMetaData is a function that adds metadata to a Seurat object. ####
filtered_lymph_sctransform <- AddMetaData(
  object = filtered_lymph_sctransform,
  metadata = predictions
) 

#### set the cell identities to the cell type predictions ####
Idents(filtered_lymph_sctransform) <- "predicted.id"

head(filtered_lymph_sctransform@meta.data)


#### ----------------------------------------------------------------------------------------

#### Clustering using the combined RNA and ATAC data ####
FindClusters(filtered_lymph_sctransform, graph.name = "wsnn", algorithm = 3, resolution = 10, verbose = FALSE)

saveRDS(filtered_lymph_sctransform, file="../results/filtered_lymph_sctransform_after_annotation_final.rds")
fs::file_size("../results/filtered_lymph_sctransform_after_annotation_final_final.rds")

rm(reference, filtered_lymph)

filtered_lymph_sctransform <- readRDS('../results/filtered_lymph_sctransform_after_annotation_final.rds')

#### print the predicted cell types identities
table(Idents(filtered_lymph_sctransform)) 

# Identify top 10 markers for each cluster
filtered_lymph_sctransform.markers.seurat <- FindAllMarkers(filtered_lymph_sctransform, only.pos = TRUE, logfc.threshold = 0.25, verbose = F) 

head(filtered_lymph_sctransform.markers.seurat )

top10 <- filtered_lymph_sctransform.markers.seurat %>% group_by(cluster) %>% top_n(n = 10 , wt = avg_log2FC)

DoHeatmap(filtered_lymph_sctransform, features = as.character(top10$gene)) #  plot a heatmap of the expression of a set of features for a Seurat object.


#### ----------------------------------------------------------------------------------------

#### plot known B cell markeers on umaps ####
to_plot <- c("MS4A1", "BANK1", "PAX5", 'CD19', 'CD22', 'NFATC1', 'TCF4', 'IRF8') # list of B cell markers

feat1 <- FeaturePlot(filtered_lymph_sctransform, to_plot, reduction='umap.rna_50', label.size = 3) + NoLegend() # plot b-cell markers on RNA umap
svg(filename = "../results/Bcell_markers_on_umap_rna_50.svg",  width = 12, height = 8)
feat1

feat2 <- FeaturePlot(filtered_lymph_sctransform, to_plot, reduction='umap.atac_50', label.size = 3) + NoLegend() # plot b-cell markers on ATAC umap
svg(filename = "../results/Bcell_markers_on_umap_atac_50.svg",  width = 12, height = 8)
feat2

feat3 <- FeaturePlot(filtered_lymph_sctransform, to_plot, reduction='wnn.umap_50', label.size = 3) + NoLegend() # plot b-cell markers on combined RNA and ATAC umap
svg(filename = "../results/Bcell_markers_on_combined_rna_atac_50.svg",  width = 12, height = 8)
feat3

p5 <- VlnPlot(filtered_lymph_sctransform, features = to_plot)
svg(filename = "../results/Bcell_markers_violin_plot.svg",  width = 12, height = 8)


#### plot umaps with annotations ####
dim4 <- DimPlot(filtered_lymph_sctransform, reduction = "umap.rna_50", group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 3) # plot annotations on RNA umap
svg(filename = "../results/umap_rna_50_with_annotations.svg",  width = 8, height = 6)
dim4

dim5 <- DimPlot(filtered_lymph_sctransform, reduction = "umap.atac_50", group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 3) # plot annotations on ATAC umap
svg(filename = "../results/umap_atac_50_with_annotations.svg",  width = 8, height = 6)
dim5

dim6 <- DimPlot(filtered_lymph_sctransform, reduction = "wnn.umap_50", group.by = "predicted.id", label = TRUE, repel = TRUE, label.size = 3) # plot annotations on combined RNA and ATAC umap
svg(filename = "../results/wnn_umap_50_with_annotations.svg",  width = 8, height = 6)
dim6


#### ----------------------------------------------------------------------------------------

#### Link peaks to genes. try to find set of peaks that may regulate the expression of a gene of interest. 
#### This is done by computing the correlation between gene expression and accessibility at nearby peaks, 
#### and correcting for bias due to GC content, overall accessibility, and peak size.
DefaultAssay(filtered_lymph_sctransform) <- "peaks"

#### compute the GC content for each peak
filtered_lymph_sctransform <- RegionStats(filtered_lymph_sctransform, genome = BSgenome.Hsapiens.UCSC.hg38) # RegionStats computes statistics for a set of genomic regions.

#### link peaks to genes
filtered_lymph_sctransform <- LinkPeaks(
  object = filtered_lymph_sctransform,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = to_plot
)

idents.plot <- c("B", "NK", "CD4 T", "CD8 T", "other T", "Mono", 'DC', 'other')
p1 <- CoveragePlot(
  object = filtered_lymph_sctransform,
  region = "MS4A1",
  features = "MS4A1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p3 <- CoveragePlot(
  object = filtered_lymph_sctransform,
  region = "PAX5",
  features = "PAX5",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

pat1 <- patchwork::wrap_plots(p1, p3, ncol = 1)
svg(filename = "../results/coverage_plots.svg",  width = 8, height = 12)
pat1


#### ----------------------------------------------------------------------------------------

######## Add motif information to seurat object ########

#### get the JASPAR2020 motif matrix set ####
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# AddMotifs is a function that adds motif information to a Seurat object. 
filtered_lymph_sctransform_with_motif <- AddMotifs(
  object = filtered_lymph_sctransform,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

rm(filtered_lymph_sctransform)

saveRDS(filtered_lymph_sctransform_with_motif, file="../results/filtered_lymph_sctransform_with_motif_final.rds")

filtered_lymph_sctransform_with_motif <- readRDS("../results/filtered_lymph_sctransform_with_motif_final.rds")


#### FindMarkers is a function that finds differentially accessible peaks between two groups of cells in a Seurat object. ####
da_peaks <- FindMarkers(
  object = filtered_lymph_sctransform_with_motif,
  ident.1 = 'B',  # B cells
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
) 

#### get top differentially accessible peaks ####
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

#### find motifs enriched in the sequences of a set of genomic regions using FindMotifs.
enriched.motifs <- FindMotifs(
  object = filtered_lymph_sctransform_with_motif,
  features = top.da.peak
) 

head(enriched.motifs) # print the enriched motifs ####

#### plots the enrichment of motifs in a set of genomic regions using motiplot.
motif1 <-MotifPlot(
  object = filtered_lymph_sctransform_with_motif,
  motifs = head(rownames(enriched.motifs))
)
svg(filename = "../results/motif_plot.svg",  width = 8, height = 6)
motif1
saveRDS(filtered_lymph_sctransform_with_motif, file="../results/final_capstone_project_data.rds")
