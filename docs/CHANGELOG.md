## 0.1.9

- 🐛 Load `all_config_annotations.csv` if `filtered_contig_annotations.csv` doesn't exist for `tcr.ImmunarchLoad`
- 🐛 Calculate diversity for all clones only if filtering by clone sizes failed for `tcr.ImmunarchAdvanced`
- 🚑 Fix seurat object creating when expressions are named "Gene Expression" for scrna.SeuratPreparing
- ✨ Add `tcr.TCRClustering`
- ✨ Add `raw` to immdata for `tcr.immunarchLoading`
- ✨ Add `on_raw` env to `tcr.TCRClustering`
- ✨ Add `bam.ControlFREEC`

## 0.1.8

- ✨ Add tcr.Attach2Seurat

## 0.1.7

- ➕ Add datar dep for scrna_metabolic pipeline
- 🚑 Fix scrna_metabolic.MetabolicPathwayActivity
- ✨ Add bcftools.BcftoolsFilter
- 👽️ Don't wrap job report in `report_jobs` report macro (to adopt pipen-report 0.2)
- ✨ Add more options for scrna.DimPlots

## 0.1.6

- ✨ Convert CNVpytor results to gff and bed
- 🚑 Make scrna_metabolic pipeline work standalone
- ➕ Add datar dep for scrna_metabolic pipeline
- 🚑 Fix scrna_metabolic.MetabolicPathwayActivity
- ✨ Add bcftools.BcftoolsFilter

## 0.1.5

- ✨ Add features and fix issues for immunopipe 0.0.4
- ✨ Add some vcf processes

## 0.1.4

- 🐛 Fix bam.CNVpytor when snpfile is not provided
- ✨ Add metabolic pathway analysis for single-cell RNA-seq data

## 0.1.3

- Add gsea.GSEA and scrna.SCImpute
- Add gene name conversions
- Add gsea.FGSEA
- Add venn plots and refactor ImmunarchFilter
- Add plot.Heatmap
- Reuse plot.Heatmap for scrna.GeneExpressionInvestigation
- Attach metadata to seurat object in scrna.SeuratPreparing
- Add envs.group_subset for scrna.GeneExpressionInvestigation
- Fix typo for scrna.GeneExpressionInvestigation
- Add docs


## 0.1.2

- ✨ Add envs.qc for scrna.SeuratPreparing

## 0.1.1

- Finish processes for immunopipe

## 0.1.0

- Adopt pipen 0.2+
