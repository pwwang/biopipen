## 0.4.6

- ✨  [vcf.TruvariBench] Allow `multimatch` to be passed
- ✨  [vcf.TruvariConsistency] Add report

## 0.4.5

- ✨ [bam.CNVpytor] Generate and fix VCF file as result
- 📝 [vcf.TruvariBench] Update docs to show other arguments for `truvari bench`
- ✨ [vcf.TruvariBench] Allow `sizemax` to be passed
- ✨ [bed.BedConsensus] Add process and tests
- ✨ [core] Add `ref.genome` to configurations
- ⚡️ [bed.BedConsensus] Parallelize and speed up
- 💚 [test] Add bedtools to env `bio`
- 💚 [test] Add chromsome sizes to reference
- 💚 [test] Add r-gsea_r to env `r`
- 💚 [scrna.ScFGSEA] Fix tests⏎

## 0.4.4

- 🐛 [scrna.SeuratPreparing] Fix after tidyseurat being used
- 🐛 [scrna.SeuratPreparing] Fix object `Sample` not found
- 📝 [Housekeeping] Fix API docs
- 📝 [Housekeeping] Make apis show neater docs

## 0.4.3

- ✨ [scrna] Add `filter` for cases in CellsDistribution, MarkersFinder and ScFGSEA
- ✨ [utils] Allow gg object for ggs in plot.R
- 🐛 [scrna_metabolic] Fix reports
- 🐛 [scrna_metabolic] Fix multiple cases
- 🐛 [scrna_metabolic] Fix rmagic for normalization
- ⚡️ [scrna.SeuratClusterStats] Add common gene list
- ⚡️ [scrna.MarkersFinder] Add `filter2` to filter after mutaters
- 🐛 [tcr.Immunarch] Fix missing library tibble in script
- ⚡️ [scrna.ScFGSEA] Make ident hierarchical

## 0.4.2

- 💚 [Housekeeping] Fix CI deploy
- ⚡️ [processes] Use faster do_call() instead of do.call()
- 📝 [tcr] Fix some docstrings with `{{` and `}}`
- ✅ [vcf.TruvariBench] Add ref for test
- 🩹 [tcr.TCRClustering] FIx VGeneScores.txt being generated in current dir
- 📝 [scrna.SeuratPreparing] Update docstring and refactor script
- ✨ [scrna.SeuratClustering] Allow dims to be expanded in arguments
- 📝 [scrna.MarkersFinder] Adopt reduced case configuration level

## 0.4.1

### General
- 👷 [Housekeeping] Add deploy in CI
- 🚚 [Housekeeping] Move tests/test_tcr/TCRClustering to tests/test_tcr/TCRClusteringStats
- 🔧 [Tests] Add r-tidyseurat to env_r.toml

### Processes
- 🩹 [scrna.CellsDistribution] Reduce envs.cases levels
- 🩹 [scrna.CellsDistribution] Allow acurate sizes to be used in orderby
- 🩹 [scrna.ScFGSEA] Reduce envs.cases levels
- ✨ [scrna.ScFGSEA] Allow `{ident}` or `{cluster}` as placeholder in cases
- ✨ [scrna.SeuratClusterStats] Add dimplots
- 🚑 [scrna.SeuratClusterStats] Limit 20 genes by default
- 🐛 [tcr.ImmunarchLoading] Fix multiple "Source" columns in data
- 🩹 [tcr.TCRClustering] Make clusterfile as a meta file that can be used by SeuratMetadataMutater
- ✨ [tcr.TCRClusteringStats] Add shared clusters by grouping
- 📝 [tcr.TCRClusteringStats] Don't show shared TCR clusters between groups if not configured
- 📝 [gsea.FGSEA] Limit pagesize to 10 in report
- ✨ [vcf.TruvariBenchSummary] Add process and test
- ✨ [vcf.TruvariBenchSummary] Add default `input_data`
- ✏️ [bed.Bed2Vcf] Fix typos in doc
- ✨ [bed.Bed2Vcf] Allow records to be skipped
- ✅ [vcf.TruvariBench] Add ref for test

## 0.4.0

- ✨ [scrna.CellsDistribution] Add process and test
- 🗑️ Remove `namespaces` (use `ns` instead)

## 0.3.2

- ✅ Allow tests to run locally only
- 💚 Add pipen-args for tests
- ✅ [plot.Heatmap] Fix test
- ✅ [pipeline.scrna_metabolic] Add ARGS in run.env
- ✅ [scrna.ScFGSEA] Add test
- ✨ [tcr.TCRClusteringStats] Add process
- ✅ [tcr.TCRClustering] Use env r for testing
- ✅ [tcr.TCRClustering] Add test
- ✅ [pipeline.scrna_metabolic] Add test
- ✅ [gsea.GSEA] Add tests
- ✅ [gsea.FGSEA] Add tests
- ✅ [plot.Heatmap] Add tests
- ✅ [gene.GeneNameConversion] Add tests
- ✅ [utils.gene] Add tests
- 💚 [bed.Bed2Vcf] Fix test
- ✅ [vcf.VcfFix] Add test
- ✅ [misc.File2Proc] Use base container for test
- ✅ [misc.File2Proc] Fix test
- 🩹 [scrna.ExprImpute] Use if-statement for requirements
- ✨ [scrna.SeuratClusterStats] Add process and test

## 0.3.1

- 🗑️ Deprecate `biopipen.namespaces`, use `biopipen.ns` instead
- ✨ [bed.Bed2Vcf] Add bed.Bed2Vcf
- ✨ [vcf.VcfFix] Add vcf.VcfFix
- 🐛 [vcf.vcfFix] Fix when a flag in INFO
- ✨ [vcf.TruvariBench] Add vcf.TruvariBench
- ✨ [vcf.TruvariConsistency] Add vcf.TruvariConsistency
- 🐛 [utils.reference] Fix typo in tabix_index
- 🐛 [vcf.VcfIndex] Fix vcf.VcfIndex
- ✨ [bed.Bed2Vcf] Allow to ignore non-existing contigs and index the output file
- ✨ [misc.Shell] Add misc.Shell to run a shell command

## 0.3.0

- ♻️ Refactor some processes for immunopipe
- 🩹 [scrna.SeuratPreparing] Remove tmp datadir for scrna.SeuratPreparing if exsits
- 🩹 [scrna.SeuratPreparing] Add a TODO comment in scrna.SeuratPreparing (#26)
- ✨ [scrna.Subset10X] Add `scrna.Subset10X`
- 💥 [tcr.Immunarch] Merge `tcr.ImmunarchBasic` and `tcr.ImmunarchAdvanced` into `tcr.Immunarch`
- 🩹 [tcr.VJUsage] Fix R script being generated at current direct for `tcr.VJUsage`
- ✨ [scrna.SeuratMetadataMutater] Add `scrna.SeuratMetadataMutater`
- 🐛 [tcr.Immunarch] Fix clonotype tracking not selecting top clones by given top N
- ♻️ [pipeline.scrna_metabolic] Refactor scrna_metabolic
- 📝 [pipeline.scrna_metabolic] Update docs for scrna_metabolic pipeline
- ✨ [pipeline.scrna_metabolic] Allow scrna_metabolic pipeline to handle multiple cases
- 🚑 [scrna.ExprImpute] Fix reticulate not using right python
- 🚑 [scrna.SeuratMetadataMutater] Fix error when input mutaters in None
- 🚑 [scrna_metabolic.MetabolicInputs] Fix diot not imported in script

## 0.2.1

- User rtoml over toml

## 0.2.0

- 📌 Pin deps for docs
- Don't link non-existing files for misc.Glob2Dir
- Upgrade datar to 0.8
- ⬆️ Upgrade pipen to v0.3
- ⚡️ Load 10X TCR and RNA-seq data files more robustly for scrna.SeuratPreparing and tcr.ImmunarchLoading


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
