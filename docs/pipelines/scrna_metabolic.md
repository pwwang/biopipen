# scrna_metabolic

Metabolic landscape analysis for single-cell RNA-seq data

An abstract from https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape

## Reference

> Xiao, Zhengtao, Ziwei Dai, and Jason W. Locasale. "Metabolic landscape of the tumor microenvironment at single cell resolution." Nature communications 10.1 (2019): 1-12.

## Run the pipeline

### Run from CLI:

```shell
pipen run scrna_metabolic_landscape ScrnaMetabolicLandscape [options]
```

### Serve as part of a pipeline:

```python
from biopipen.ns.scrna_metabolic_landscape import ScrnaMetabolicLandscape

pipeline = ScrnaMetabolicLandscape(<options>)

# You can specify dependencies so that the whole metabolic landscape pipeline
# works as a part of another pipeline
```

## Inputs

- `metafile`: Either a metafile or an rds file of a Seurat object.
    If it is a metafile, it should have two columns: `Sample` and `RNAData`.
    `Sample` should be the first column with unique identifiers for the samples
    and `RNAData` indicates where the barcodes, genes, expression matrices are.
    The data will be loaded and an unsupervised clustering will be done.
    Currently only 10X data is supported.
    If it is an rds file, the seurat object will be used directly.

- `gmtfile`: The GMT file with the metabolic pathways. The gene names should
    match the gene names in the gene list in RNAData or the Seurat object.
    You can also provide a URL to the GMT file.

- `group_by`: Group the data by the given column in the metadata. For example, `cluster`.

- `subset_by`: (Optional) Subset the data by the given column in the metadata.
    For example, `Response`. `NA` values will be removed in this column.
    If None, the data will not be subsetted.

- `mutaters`: (Optional) Add new columns to the metadata for grouping/subsetting.
    They are passed to `sobj@meta.data |> mutate(...)`. For example,
    `{"timepoint": "if_else(treatment == 'control', 'pre', 'post')"}`
    will add a new column `timepoint` to the metadata with values of
    `pre` and `post` based on the `treatment` column.

- `cases`: (Optional) Multiple cases for the analysis. If you have multiple different
    grouping/subsetting scenarios, you can specify them here. Each case can have
    its own `subset_by`, `group_by`, and analysis parameters.

## Advanced Configuration

### Multiple Cases

You can define multiple analysis cases with different grouping/subsetting strategies:

```toml
[ScrnaMetabolicLandscape]
metafile = "test_data/scrna_metabolic/seurat_obj.rds"
gmtfile = "test_data/scrna_metabolic/KEGG_metabolism.gmt"

[ScrnaMetabolicLandscape.MetabolicPathwayActivity.envs.cases]
"By Treatment" = { group_by = "seurat_clusters", subset_by = "treatment" }
"By Response" = { group_by = "seurat_clusters", subset_by = "response" }
```

### Custom Plotting

Each process supports customizable plots:

```toml
[ScrnaMetabolicLandscape.MetabolicPathwayActivity.envs.plots]
"Custom Heatmap" = {
    plot_type = "heatmap",
    show_row_names = true,
    devpars = { width = 1200, height = 800, res = 150 }
}
"Custom Violin" = {
    plot_type = "violin",
    add_box = true,
    devpars = { width = 1000, height = 600 }
}

[ScrnaMetabolicLandscape.MetabolicFeatures.envs.plots]
"Top 5 Summary" = {
    plot_type = "summary",
    top_term = 5,
    level = "subset"
}
```

## A step-by-step example

### Prepare the seurat object

Using the data from:

> Yost KE, Satpathy AT, Wells DK, Qi Y et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. Nat Med 2019 Aug;25(8):1251-1259. PMID: 31359002

```r
library(Seurat)

# Download data (tcell rds and metadata) from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813

count_file <- "test_data/scrna_metabolic/GSE123813_bcc_scRNA_counts.txt.gz"
meta_file <- "test_data/scrna_metabolic/GSE123813_bcc_all_metadata.txt.gz"

counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = F)
metadata <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = F)
```

```r
# Subset 1000 cells for just demo purpose
counts = counts[, sample(1:ncol(counts), 1000)]
metadata = metadata[colnames(counts),]
```

```r
# Create seurat object
seurat_obj = CreateSeuratObject(counts=counts)
seurat_obj@meta.data = cbind(
    seurat_obj@meta.data,
    metadata[rownames(seurat_obj@meta.data),]
)
seurat_obj = NormalizeData(seurat_obj)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- FindVariableFeatures(object = seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Output exceeds the size limit. Open the full output data in a text editor

# Warning message:
# "Feature names cannot have underscores ('_'), replacing with dashes ('-')"
# Centering and scaling data matrix

# PC_ 1
# Positive:  CALD1, EMP2, SPARC, CAV1, EMP1, ACTN1, KRT17, BGN, DSP, RAB13
# 	   RND3, S100A16, KRT14, S100A14, S100A2, CD9, FHL2, IER3, GJB2, TRIM29
# 	   TSC22D1, KRT15, SFN, FSTL1, DDR1, IGFBP7, JUP, PTRF, KRT5, FOSL1
# Negative:  CD74, CRIP1, RARRES3, RGS1, LAT, CD69, NKG7, HLA-DPA1, TNFRSF4, CD27
# 	   GZMA, HLA-DRB1, CXCR6, CTSW, GPR183, LDLRAD4, ICOS, HIST1H4C, GZMK, AC092580.4
# 	   SLC9A3R1, CCR7, S100A4, HLA-DPB1, HMGB2, GBP5, SELL, CXCR3, LAG3, HLA-DRB5
# PC_ 2
# Positive:  DSP, DSC3, SFN, SERPINB5, KRT17, S100A14, KRT15, KRT16, IRF6, TRIM29
# 	   KRT14, LGALS7B, JUP, PERP, TACSTD2, GJB2, KRT5, KRT6B, DDR1, S100A2
# 	   PKP1, CDH3, KRT6A, FXYD3, GJB3, MPZL2, CXADR, DSC2, DSG3, GRHL3
# Negative:  COL1A2, COL3A1, LUM, COL6A1, MXRA8, FN1, CTSK, COL1A1, COL6A3, DCN
# 	   MMP2, COL6A2, PRRX1, FKBP10, TNFAIP6, FAP, PCOLCE, PDGFRB, NNMT, AEBP1
# 	   C1S, CCDC80, SFRP2, RCN3, PDPN, SERPINF1, COL5A2, CTHRC1, COL5A1, COL12A1
# PC_ 3
# Positive:  LYZ, TYROBP, SPI1, FCER1G, KYNU, C15orf48, BCL2A1, HLA-DRA, CD68, AIF1
# 	   CST3, CTSZ, SERPINA1, CSF2RA, HLA-DRB5, SLC7A11, LST1, MS4A7, ALDH2, FAM49A
# 	   IFI30, HLA-DRB1, GPR157, PLEK, HLA-DPB1, IL1B, CD86, CCDC88A, HLA-DPA1, RNF144B
# Negative:  MT2A, MT1X, BGN, MT1E, COL6A2, COL1A2, COL6A1, COL5A2, MXRA8, DCN
# 	   COL12A1, COL3A1, COL1A1, LUM, MFAP2, PCOLCE, MT1F, MMP2, COL5A1, COL6A3
# 	   PRRX1, C1S, AEBP1, CTSK, FAP, LAT, PDGFRB, RARRES3, THY1, EFEMP2
# 	   GJB3, LYPD3, GRHL3, PVRL4, KRT16, TACSTD2, GJB5, SERPINB13, MPZL2, KRT23
# Negative:  UBE2C, GTSE1, BIRC5, RRM2, CCNA2, TYMS, TOP2A, TK1, DLGAP5, MKI67
# 	   PKMYT1, CENPA, KIFC1, CDCA5, UHRF1, ASF1B, AURKB, FAM111B, TROAP, CKAP2L
# 	   HJURP, ESCO2, FOXM1, CDK1, ZWINT, E2F2, CLSPN, HIST1H1B, CDT1, MCM10
```

```r
# By default, the pipeline will do the clustering using the SeuratClustering process
# If you want to do the clustering yourself, you can set `is_seurat = true`
# when running the pipeline, which will skip the clustering step.

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
head(Idents(seurat_obj))
# Computing nearest neighbor graph

# Computing SNN

# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# Number of nodes: 1000
# Number of edges: 30631

# Running Louvain algorithm...
# Maximum modularity in 10 random starts: 0.8795
# Number of communities: 9
# Elapsed time: 0 seconds

# bcc.su009.post.tcell_CCATTCGCAATCTACG 0
# bcc.su004.pre.tcell_ACAGCTACACTTCTGC  0
# bcc.su006.pre.tcell_CGTGTAAAGTGACTCT  1
# bcc.su001.post.tcell_CCTTTCTGTACCGTTA 2
# bcc.su007.pre.tcell_ATTTCTGAGAAGGACA  1
# bcc.su009.pre.tcell_AGACGTTTCCTGCAGG8 8
# Levels:
# '0''1''2''3''4''5''6''7''8'
```

```r
# Check the meta.data
head(seurat_obj@meta.data)
# 	orig.ident	nCount_RNA	nFeature_RNA	patient	treatment	sort	cluster	UMAP1	UMAP2	RNA_snn_res.0.5	seurat_clusters
# <fct>	<dbl>	<int>	<chr>	<chr>	<chr>	<chr>	<dbl>	<dbl>	<fct>	<fct>
# bcc.su009.post.tcell_CCATTCGCAATCTACG	bcc.su009.post.tcell	4242	1726	su009	post	CD45+ CD3+	CD8_mem_T_cells	-9.1816435	0.7484789	0	0
# bcc.su004.pre.tcell_ACAGCTACACTTCTGC	bcc.su004.pre.tcell	3159	1466	su004	pre	CD45+ CD3+	CD8_ex_T_cells	4.1067562	3.3754938	0	0
# bcc.su006.pre.tcell_CGTGTAAAGTGACTCT	bcc.su006.pre.tcell	3279	1401	su006	pre	CD45+ CD3+	Tregs	0.4816457	11.9388428	1	1
# bcc.su001.post.tcell_CCTTTCTGTACCGTTA	bcc.su001.post.tcell	5057	2218	su001	post	CD45+ CD3+	CD8_ex_T_cells	2.3045983	7.4856248	2	2
# bcc.su007.pre.tcell_ATTTCTGAGAAGGACA	bcc.su007.pre.tcell	3701	1413	su007	pre	CD45+ CD3+	CD8_mem_T_cells	-1.9617293	5.5546365	1	1
# bcc.su009.pre.tcell_AGACGTTTCCTGCAGG	bcc.su009.pre.tcell	3891	1593	su009	pre	CD45+ CD3+	Tcell_prolif	5.2243228	-1.0460106	8	8
```

```r
# save seurat object
saveRDS(seurat_obj, "test_data/scrna_metabolic/seurat_obj.rds")
```

### Prepare the pathway file

A set of collected metabolic pathways can be found here:

https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape/blob/master/Data/KEGG_metabolism.gmt

Download and save it to `test_data/scrna_metabolic/KEGG_metabolism.gmt`

### Prepare the configuration file

Save at `test_data/scrna_metabolic/config.toml`:

```toml
# Pipeline configuration
[ScrnaMetabolicLandscape]
metafile = "test_data/scrna_metabolic/seurat_obj.rds"
gmtfile = "test_data/scrna_metabolic/KEGG_metabolism.gmt"
group_by = "seurat_clusters"
subset_by = "timepoint"

[ScrnaMetabolicLandscape.mutaters]
timepoint = "if_else(patient != 'su001', NA_character_, treatment)"

# Optional: Configure individual processes
[ScrnaMetabolicLandscape.MetabolicPathwayActivity.envs]
ntimes = 1000

[ScrnaMetabolicLandscape.MetabolicFeatures.envs.plots]
"Summary Plot" = { plot_type = "summary", top_term = 5 }
```

### Run the pipeline

```shell
pipen run scrna_metabolic_landscape ScrnaMetabolicLandscape --config test_data/scrna_metabolic/config.toml
```

### Check out the results/reports

The results can be found at `./ScrnaMetabolicLandscape_results/`, and reports can be found at `./ScrnaMetabolicLandscape_results/REPORTS`. To check out the reports, open `./ScrnaMetabolicLandscape_results/REPORTS/index.html` in your browser.

There are 3 parts of the results:

- `MetabolicPathwayActivity`:

    The pathway activities for groups (defined by `group_by` in the configuration) for each subset (defined by `subset_by`)

- `MetabolicPathwayHeterogeneity`:

    Pathway heterogeneity for groups for each subset

- `MetabolicFeatures`:

    The pathway enrichment analysis in detail for each group against the rest of the groups in each subset.
