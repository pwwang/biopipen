# Change Log

## 0.34.4

- chore(scrna): update error flag in MarkersFinder and PseudoBulkDEG classes
- chore(deps): bump pipen to 0.17.11

## 0.34.3

- BREAKING(scrna.MetabolicFeatures): use `:` instead of `,` to separate groups for comparisons
- feat(tcr.ScRepCombiningExpression): add TCR presence indicator to combined expression object
- fix(scrna.CellTypeAnnotation): fix identification of default ident if rds or qs is given as input for celltypist
- fix(scrna.CellTypeAnnotation): ensure n_neighbors is set in neighbors params
- fix(scrna.SeuratClusterStats): correct typo in group_by variable name
- fix(tcr.ScRepLoading): fix an issue where barcode is not the first column when loaded for 10X data
- fix(tcr.TCRClustering): fix input collection for clustering
- fix(tcr.TESSA): refine TCR input preparation
- fix(scrna.SeuratClusterStats): exclude NA entities in groupings
- chore(scrna.CellTypeAnnotation): fix dead links for model files for celltypist
- chore: remove SCP-plot.R
- chore(ci): comment out cache deletion condition
- chore(scrna_metabolic_landscape): enlarge font size for list items in introduction of report templates
- chore(scrna.SeuratClusterStats): adjust plot settings for better visualization
- chore(tcr.TCRClustering): swap importing BLOSUM62 matrix from newer and older version of biopython for GIANA
- chore(scrna_metabolic_landscape): adjust plot dimensions
- chore(scrna.SeuratClusterStats): enhance descriptions of plots
- chore(deps): update versions for copier, pipen-report, and pipen-verbose
- docs(scrna): update mutater documentation to include clone selectors
- ci: update cache number for conda environments to force updates

## 0.34.2

- feat(scrna): add PseudoBulkDEG process for differential gene expression analysis
- fix(test.Seurat): update pipeline function to set starts for PrepareSeurat
- chore(scrna.ScFGSEA): rename allpathway_plots to alleach_plots
- chore(scrna/tcr): update parameter naming for consistency
- chore(scrna.TopExpressingGenes): use common report template
- chore(test): update pipeline function to remove report enabling
- test(scrna): merge tests for Seurat processes to avoid repeatedly load pbmc3k dataset
- test(scrna): move map2ref tests out of tests/scrna/Seurat for being tested locally only
- test(scrna): add set.seed to PrepareSeurat script
- ci(docker): update Dockerfiles to use dynamic REF_NAME argument
- ci: update build condition and add caching for test intermediates

## 0.34.1

- fix(tcr.ClonalStats): update envs assignment to handle todot parameter
- chore(deps): bump up pipen to 0.17.8
- fix(SeuratClusterStats): fix when features are given as a dict (used in heatmap)
- feat(MarkersFinder): add enrichment plot across all subsets by each or all ident.1 in group.by
- feat(ScFGSEA): add support for all pathway plots for all subsets by each
- chore(TopExpressingGenes): adjust plot height for bar plots
- docs(SeuratPreparing): update docs for support loading loom files

## 0.34.0

### New Features
- feat(scrna): add ScVelo analysis for RNA velocity and Slingshot for trajectory inference
- feat(tcr): add ScRepCombiningExpression for combining TCR/BCR and expression data, and ScRepLoading for multiple TCR/BCR data formats
- feat(scrna): add ScRepLoading with support for multiple formats and improved logging
- feat(plot): add Plot class and associated R script for data visualization
- feat(bam): add SamtoolsView for BAM file processing
- feat(utils): implement Reporter class for generating JSON reports for processes
- feat(scrna): support qs2 format for input and output in various processes
- feat(scrna.SeuratPreparing): add mutaters parameter for metadata mutation
- feat(scrna.MarkersFinder): enhance parameter handling and marker processing
- feat: add common Svelte report template for job reporting using JSON reports

### Enhancements
- enh(bam.CNAClinic): change `envs.binsize` to bp instead of kbp
- enh(cnv): replace ggplot2 with plotthis for improved plotting in AneuploidyScoreSummary and TMADScoreSummary
- enh(scrna): adopt biopipen.utils.R v0.1.0 across multiple processes
- enh(scrna_metabolic_landscape): improve flexibility and stability of metabolic landscape analysis
- enh(scrna.ExprImputation): improve error handling and threshold handling for cell imputation
- enh(tcr.ClonalStats): support qs2 format for output

### Bug Fixes
- fix(cnv.TMADScore): correct output filename
- fix(bam): fix report template and argument handling in CNVpytor, CNAClinic, and ControlFREEC
- fix(cnv.AneuploidyScore): replace ggplot with plotthis and add error handling for chromosome detection
- fix(cnv.AneuploidyScoreSummary): rename 'rows' to 'rows_by' for clarity in heatmap function
- fix(delim.SampleInfo): fix plot functions being registered twice to gglogger and ensure reporter saves to correct directory
- fix(scrna): replace readRDS/saveRDS with biopipen.utils functions for consistency across processes
- fix(scrna.SeuratClusterStats): improve data handling and plot saving functionality

### Refactoring
- refactor(scrna.CellCellCommunicationPlots): use scplotter::CCCPlot
- refactor(tcr.CDR3AAPhyschem): adopt input from ScRepCombiningExpression
- refactor: remove utility R scripts and use biopipen.utils.R package

### Development & Infrastructure
- chore(deps): update dependencies and specify versions for bioconductor packages
- chore(docker): refactor Dockerfiles to streamline base image usage and dependency installation
- chore: use filter 'r' for R input/output paths instead of 'quote' in scripts
- test: add comprehensive tests for bam, cnv, and scrna processes (local tests only for some)
- test: add docker-test job and Dockerfile for test image building
- ci: update environment cache and reference data cache

## 0.33.1

- fix(delim.SampleInfo): fix when plot_type ending with "plot"
- feat(scrna.LoomTo10X): add LoomTo10X to convert loom format of scRNA-seq data to 10X format
- ci: update conditions for build and deploy jobs based on event type

## 0.33.0

- chore(dependencies): update pipen (v0.17) and related package versions in pyproject.toml
- fix: update all template filters in script to adopt pipen 0.17, which passes in.file etc as a MountedPath
- ci: add caching for conda environments to improve workflow efficiency
- test: improve test output grouping for better readability
- chore: add descriptive summaries for fgsea and enrichr results (#158)
- chore(snp.PlinkFromVcf): enhance type annotations and set default for keep_allele_order
- feat(snp.Plink2GTMat): enhance genotype coding options and improve documentation
- feat(stats.ChowTest): separate groups in output and add pvalues for the coefficient for the subregressions
- fix(utils/misc.py): enhance error messages in command execution for better debugging
- fix(web.Download): enhance output filename generation by adding URL decoding and improved slugification
- fix(cellranger.CellRangerCount): fix inconsistency between in.id and in.fastqs
- chore(snp.PlinkFilter): remove unnecessary docstring from PlinkFilter script
- fix(cellranger.CellRangerSummary): use plotthis and biopipen.utils.R for plotting, logging and report content generation
- test(cellranger): add tests for CellRangerCount and CellRangerSummary with data download setup
- chore: comment out dev-dependencies section in pyproject.toml
- feat(scrna.CellCellCommunication): add subset and split_by options for CellCellCommunication and update conversion logic

## 0.32.3

- chore: add descriptive summaries for fgsea and enrichr results

## 0.32.2

- chore: update dependencies to latest versions
- feat: add PDF output option for SampleInfo plots
- feat: add PDF output options for violin and scatter plots in Seurat preparation scripts
- feat: add PDF output options for volcano, dotplot, venn, and upset plots; update filters for type hints
- feat: add PDF output option for Enrichr plots in TopExpressingGenes script
- feat: add PDF output options for UMAP plots in SeuratMap2Ref script; update image handling in misc.liq
- feat: add PDF output options for cluster size distribution, shared clusters, and sample diversity plots; update plotting functions to handle multiple output formats
- feat: add PDF output options for various Immunarch scripts; enhance reporting with downloadable PDF files
- feat: add PDF output options for cluster size distribution, dimension plots, and feature plots; enhance reporting with downloadable PDF files
- feat: add PDF output options for radar and bar plots; enhance reporting with downloadable PDF files
- feat: add PDF output options for CloneResidency script; enhance reporting with downloadable PDF files
- feat: add PDF output options for GSEA table and enrichment plots; enhance reporting with downloadable PDF files
- feat: add PDF output options for pie charts, heatmaps, Venn plots, and UpSet plots; enhance reporting with downloadable PDF files
- feat: add PDF output options for Enrichr plots; enhance reporting with downloadable PDF files
- feat: add PDF output options for estimated coefficients and distribution plots; enhance reporting with downloadable PDF files
- chore: add gcc to cnvkit pipeline docker deps

## 0.32.1

- fix(scrna.ScFGSEA): fix case gmtfile not working
- fix(TopExpressingGenes): add `InlineNotification` component to TopExpressingGenes.svelte
- feat(scrna.SeuratPreparing): add `envs.species` so that `percent.mt`, `percent.ribo`, `percent.hb` and `percent.plat` can be correctly calculated for mouse
- fix(scrna.SeuratClusterStats): fix `kind` not being added to the figure file name for plots of features

## 0.32.0

- deps: update pipen-runinfo dependency to version 0.8.1
- feat(scrna): add CellCellCommunication and CellCellCommunicationPlots
- fix(scrna.SeuratMap2Ref): fix report and add stats to report
- fix(utils.single_cell.R): fix categorical data when converting seurat to anndata format
- refactor(scrna.Seurat2AnnData): abstract seurat_to_anndata() for reuse
- enh(tcr.TCRClustering): make GIANA compatible with latest BioPython (v1.84)
- fix(tcr.TCRClstering): fix clusTCR error due to scipy update (v1.14)

## 0.31.7

- deps: bump pipen-args to 0.16
- chore: update pyright configuration to include `biopipen/**/*.py`
- feat(bam): add BamSubsetByBed process for subsetting bam file by regions in a bed file
- feat(bed): add BedtoolsMakeWindows process for generating windows from a BED file or genome size file

## 0.31.6

- deps: pin the channels of conda dependencies for tests
- feat(vcf): adopt truvari v4+ for related processes
- feat(regulatory): add VariantMotifPlot to plot motif and surrounding sequences with mutations
- refactor(regulatory.MotifAffinityTest): optimize code base
- ci: add verbosal output for tests

## 0.31.5

- deps: update pipen to version 0.15.3 and xqute to version 0.5.2
- feat(bam): add BamSampling process for sampling a fraction of reads from a bam file
- feat(protein): add the protein module and Prodigy and ProdigySummary to calculate the binding affinity of a complex structure
- ci: do not print verbose logs for tests
- chore(bam.BamMerge): use logger instead of print for logging

## 0.31.4

- deps: bump pipen-report to 0.20.1 (pipen to 0.15.2)
- fix(plot.VennDiagram): update default devpars and fix issues with computed data
- fix(scrna.SeuratMap2Ref): fix identifying the normalization method of reference

## 0.31.3

- test: fix test not failing when tests failed
- test: fix gene name conversion tests due to external API change
- fix(tcr.CDR3AAPhyschem): fix when chain is not available
- fix(tcr.TCRClustering): fix when chain is not available

## 0.31.2

- fix(tcr.CDR3AAPhyschem): use sequence from TRB chain only
- fix(tcr.TCRClustering): fix for multi-chain TCRs, use TRB only if on_multi is false

## 0.31.1

- enh(scrna.SeuratMap2Ref): check if reference has SCTModel if SCTransform'ed (likely prepared by old Seurat)

## 0.31.0

- deps: bump pipen to 0.15.0

## 0.30.0

### scrna/tcr

- BREAKING(scrna): move clustree plots from `SeuratClustering/SeuratSubClustering` to `SeuratClusterStats`
- feat(scrna.CellTypeAnnotation): allow to merge/not to merge (envs.merge) the clusters with the same labels predicted
- feat(scrna.SeuratPreparing): add scDblFinder to detect doublets
- feat(scrna.SeuratMap2Ref): add `envs.skip_if_normalized` option to skip normalization if query is already normalized using the same method as the reference
- refactor(tcr.Immunarch): source the files for Immunarch scripts for better debugging
- refactor(scnra.SeuratClustering): refactor the script for better debugging
- refactor(scnra.SeuratPreparing): refactor the script for better debugging
- fix(scrna): fix resolution expansion for `SeuratClustering` and `SeuratSubClustering`
- fix(cellranger.CellRangerCount): fix falsy `envs.create_bam` not working for cellranger v7
- fix(scrna): Fix generating PrepSCTFindMarkers  command when no previous commands present
- tests(scrna.ScFGSEA): fix unavailable urls to GMT files
- chore(scrna.SeuratMap2Ref): optimize memory usage
- chore(scrna.MetaMarkers): remove `plugin_opts.poplog_max`
- chore(tcr.CloneResidency): improve logging when handling subjects

### other

- fix(stats.Mediation): fix when NAs in the data
- feat(plot): add Scatter for scatter plots
- tests: use single conda env for tests
- ci: fix CI due to conda env changes
- docs(web): update docs of `envs.tool` for Download/DownloadList
- feat(web): add GCloudStorageDownloadFile and GCloudStorageDownloadBucket to download files from GCP
- chore(regulatory.MotifAffinityTest): use template filter source_r to source R files
- tests(regulatory.MotifAffinityTest): rename regulation to regulatory
- chore: use template filter source_r to source R files
- fix(stats): handle case when p-value is 0 for MetaPvalue and MetaPvalue1

## 0.29.2

- chore(stats.Mediation): make better logging strategy for various number of cases
- chore(scrna.SeuratClusterStats): use ident label length to adjust default height for feature plots
- fix(scrna.MetaMarkers): fix seurat object not updated when expanding cases and run `PrepSCTFindMarkers` when necessary before calling meta-markers
- fix(scrna.MarkersFinder): fix fetching command when composing the PrepSCTFindMarkers command
- fix(scrna_metabolic_landscape): handle null values in for loop in MetabolicFeatures and MetabolicFeaturesIntraSubset for report generation

## 0.29.1

- BREAKING: rename namespace `regulation` to `regulatory`
- choir(plot.Manhattan): default `envs.title` to `None` (don't add title to the plot by default)
- enh(plot.Manhattan): give warnings instead of errors about zooming chromosomes not existing
- fix(plot.Manhattan): fix `envs.ylabel` not working
- feat(stats): add `Mediation` for mediation analysis
- feat(plot.QQPlot): add support for custom theoratical values
- tests(plot.QQPlot): add tests
- chore(snp.MatrixEQTL): allow pvalue cutoffs to be greater than 1 (but 1 will be used anyway)
- fix(snp.PlinkIBD): add `--keep-allele-order` to keep the allele order
- fix(delim.SampleInfo): fix numbers not split up when `each` is specified.
- enh(delim.SampleInfo): make sizes of pie charts proportional to number of samples when `each` is specified
- enh(scrna.MarkersFinder): run `PrepSCTFindMarkers` when necessary before calling FindMarkers
- feat(scrna.SeuratPreparing): add option to cache Seurat object at different steps
- feat(scrna.SeuratPreparing): allow `doubletfinder` to run with a different number of cores
- chore(scrna.SeuratClustering): record `PrepSCTFindMarkers` command in `sobj@commands`
- tests(scrna.SeuratClusterStats): use less stringent p-value cutoff for DEG/MarkersFinder
- tests(scrna.SeuratPreparing): add `doubletfinder` in tests

## 0.29.0

### Depedencies

- deps: update `pipen-filters` to version 0.13
- deps: add `meme` to `env_bio.yml` for tests
- deps: bump pipen-board to 0.15.2
- ci: update GitHub Actions versions and dependencies

### BREAKING

- BREAKING: merge namespace bcftools to vcf

### New Features/Processes

- feat(snp): add `PlinkFromVcf` to convert VCF files to PLINK format
- feat(snp): add `Plink2GTMat` for converting PLINK files to genotype matrix
- feat(snp): add `PlinkIBD` analysis for identity by descent
- feat(snp): add `PlinkHWE` for Hardy-Weinberg Equilibrium test
- feat(snp): add `PlinkHet` for calculating sample heterozygosity
- feat(snp): add `PlinkCallRate` for calculating call rate of samples and variants
- feat(snp): add `PlinkFilter` process for filtering PLINK files
- feat(snp): add `PlinkFreq` process for calculating and filtering by allele frequencies
- feat(snp): add `PlinkUpdateName` proc to update variant names in PLINK files
- feat(gene): add `GenePromoters` for retrieving gene promoter regions and tests for it
- feat(bed): add `BedtoolsIntersect` process for finding the intersection of two BED files
- feat(regulation): add the namespace and `MotifScan` to use fimo for motif scanning
- feat(regulation): add `MotifAffinityTest` to test the affinity of motifs to the sequences and the affinity change due the mutations.
- feat(cnv.AneuploidyScore): allow BED and VCF files as in.segfile
- feat(cnv.TMADScore): allow BED and VCF files as in.segfile
- feat(cnvkit): allow user home directory (~) to be used in `envs.ref` in mulitple processes
- feat(plot): add `ManhattanPlot` to for support for plotting Manhattan plots
- feat(plot): add `QQPlot` proc for generating QQ-plot or PP-plot
- feat(vcf): add `BcftoolsView` process for viewing, subsetting, and filtering VCF files
- feat(vcf): add `run_bcftools` function for running bcftools with given arguments
- feat(vcf.BcftoolsSort): allow sorting contigs based on a chrom size file
- feat(vcf.BcftoolsFilter): allow indexing output file
- feat(vcf.BcftoolsAnnotate): allow providing annotation file as input file and allow indexing output file
- feat(stats): add `MetaPvalue1` to combine pvalues from the same file
- feat(stats.MetaPvalue): add `envs.keep_single` flag to keep the single p-values
- feat(utils.misc.R): add `run_command` function for R
- feat(utils.reference): allow `tabix_index` to index infile directly
- feat(snp.MatrixEQTL): add envs.match_samples flag to subset snp, expr and cov data with common samples
- feat(snp.MatrixEQTL): fix cov data being wrongly transposed
- feat(snp.MatrixEQTL): allow extra columns when snp and gene position file is BED
- feat(tests): add lazy loading for reference data download and --local flag for downloading more references locally

### Refactors/Improvements

- refactor(utils.gene): redesign `gene_name_conversion` functions for both python and R
- refactor(gene.GeneNameConversion): use R for implementation
- refactor(misc.Shell): save envs.cmd to a file and run it to fix the escaping issues of the command
- enh(snp.MatrixEQTL): use rtracklayer to read the position files

### Minor

- choir(misc.Str2File): add default for `in.name`
- choir(gene.GeneNameConversion): allow `envs.notfound` to be `ignore` or `skip` when `envs.output` is `append` or `replace`
- choir(utils.misc.py): flush output for command printing in `run_command` function
- choir(utils.misc.py): print command with a new line for `run_command` function
- choir(snp.PlinkFromVcf): indicate sex is to be handled
- choir(cellranger_pipeline): remove unused import
- choir(bam.CNVpytor): implement cnvnator2vcf directly instead of using cnvnator2vcf.pl

### Tests

- tests(utils.gene): update tests for `gene_name_conversion`
- tests(gene.GeneNameConversion): use right environment for tests
- tests(snp): add tests for plink related processes
- tests(snp): disable report generation for plink related tests
- tests(regulation): specify envs.genome for MotifAffinityTest
- tests: add python package mygene to conda environment `biopipen-r`
- tests: add tests for bcftools processes
- tests: do not download reference data for hg38 at CI
- tests: update `bioconductor-ggmanh` dependency to version 1.9.6
- tests: add bcftools to conda environment dependencies for tests

### Docs

- docs(MatrixEQTL): fix choice items of `envs.model`
- docs(cellranger_pipeline): fix types of some items in docs, which should be 'list', instead of 'type=list'

### Fixes

- fix(utils.reference): avoid index file to be created again for the same infile for `tabix_index` function
- fix(utils.reference): pass `-f` to bgzip or gunzip to overwrite the output if exists
- fix(vcf): fix passing vcffile as a string in fix_vcffile in VcfFix_utils.py
- fix(cnvkit_pipeline): fix sex in process channels
- fix(cnv.AneuploidyScoreSummary): fix when Sample column is already in metafile
- fix(cnv.TMADScoreSummary): fix when Sample column is already in metafile

### Immunopipe-related

- fix(tcr.TCRClusterStats): fix `envs.shared_clusters.heatmap_meta` being broken by `envs.shared_clusters.sample_order`
- choir(scrna.SeuratMap2Ref): present better error message when `envs.use` or values of `envs.MapQuery.refdata` not in reference
- fix(scrna.MarkersFinder): run PrepSCTFindMarkers when needed
- choir(scrna.SeuratClustering): use `FindClusters` to run for multiple resolutions
- choir(scrna.SeuratSubClustering): use `FindClusters` to run for multiple resolutions
- feat(scrna.SeuratClustering): add clustree plot
- feat(scrna.SeuratSubClustering): add clustree plot
- tests(scrna.SeuratClusterStats): add assertion for clustree plot generation

## 0.28.1

- fix(scrna.CellsDistribution): fix `devpars` and `hm_devpars` not working

## 0.28.0

- tests(scrna.CellTypeAnnotation): add tests for `CellTypeAnnotation` using `scCATCH`
- feat(cellranger_pipeline): add docker image building for cellranger pipeline
- chore(cellranger.CellRangerCount): add `envs.create_bam` to control whether create bams (supporting cellranger v8)
- chore(cellranger): add `in.id` to `CellRangerCount` and `CellRangerVdj` to specify sample ids
- chore(cellranger.CellRangerSummary): set the default value of `report_paging` to 8

## 0.27.9

- feat(tcr.TCRClusterStats): add `sample_order` to set sample order on heatmap and `cluster_rows` to switch row clustering on/off

## 0.27.8

- fix(scrna.SeuratClusterStats): fix selected columns not unique for stats
- feat(scrna.SeuratMap2Ref): allow non-SCTransform'ed reference
- feat(scrna.SeuratMap2Ref): allow splitting query object for mapping (pwwang/immunopipe#61)
- deps: update pipen-board to version 0.15.1

## 0.27.7

- fix(utils.gsea): fix gsea table not being printed for `runFGSEA`
- fix(core.filters): fix slugified pathway plot file name in report
- fix(scrna_metabolic_landscape): fix mutaters not working
- fix(scrna_metabolic_landscape.MetabolicFeatures/MetabolicFeaturesIntraSubset): skip groups with less than 5 cells in `do_one_group` and save a warning file under the case
- fix(utils.gsea): do not switch 1st and 2nd columns when 2nd column is numeric for `localizeGmtfile`
- chore: fix typo in class name ExprImpution to ExprImputation
- choir(tests): remove KEGG_metabolism.gmt for `prep_reference.py`
- tests(scrna_metabolic_landscape): fix tests

## 0.27.6

- fix(scrna_metabolic_landscape.MetabolicFeatures): fix return value of groups with less than 5 cells in `do_one_group`
- choir(utils.gsea): avoid printing NULL for `runFGSEA`
- tests: use the return value of `pipen.run()` to test the success

## 0.27.5

- fix(scrna.Subset10X/SeuratTo10X): correct the paths to the scripts
- feat(testing): allow to enable report for testing pipelines
- feat(scrna.SeuratPreparing): add `envs.cell_qc_per_sample` to filter cells before merging instead after
- test: add tests to scrna.SeuratTo10X and scrna.SeuratPreparing
- fix(scrna.SeuratClusterStats): fix color palette for ridge plots

## 0.27.4

- feat: add plot.ROC
- choir(delim.SampleInfo): add alpha to the colors of the plots using biopipen color pallete
- feat: add snp.MatrixEQTL
- docs(tcr/scrna/scrna_metabolic_landscape): update links of images in docs

## 0.27.3

- deps: temporary fix copier breaks with pyyaml-include v2 (copier-org/copier#1568)
- deps: bump pipen-poplog to 0.1.2 (quick fix for populating logs when job fails)
- choir(scrna.ScFGSEA): Skip cases when no cells found (pwwang/immunopipe#50)
- choir(scrna.MarkersFinder): Skip cases when no cells found (pwwang/immunopipe#50)
- choir(scrna.MetaMarkers): Skip cases when no cells found (pwwang/immunopipe#50)
- feat(scrna.SeuratPreparing): support `DoubletFinder`

## 0.27.2

- fix(utils.misc.py): inherit envs when `env` passed for `run_command()`
- fix(scrna.RadarPlots): fix mutaters not working
- feat(tcr.CloneResidency): support `envs.upset_ymax` to set the max value of y axis in upset bar plot.
- feat(tcr.TCRDock): add process
- choir(utils.misc.py): update level to DEBUG for python logger (leaving the filtering to pipen-poplog)
- choir(stats.DiffCoexpr): change log_warn to debug for some verbosal logging messages
- refactor(snp.PlinkSimulation): make the configuration files as input so multiple simulations could run in parallel easily.

## 0.27.1

- BREAKING(scrna.SeuratMap2Ref): rename `envs.name` to `envs.ident` so `envs.MapQuery.refdata` is not required anymore. It will be inferred from `envs.ident` and `envs.use`.

## 0.27.0

- deps: bump pipen to 0.14.5
- deps: bump datar to 0.15.6
- depr(scrna.MarkersFinder): remove `use_presto` as it's used by Seurat v5 by default
- enh(tcr.CloneResidency): support log scale for y axis of upset bar plots
- enh(scrna.SeuratClusterStats): allow to rotate labels in circos plot (pwwang/immunopipe#48)
- enh(scrna.SeuratClusterStats): use pal_biopipen for ident colors in circos plot
- fix(scrna.CellsDistribution): fix the row order of the heatmaps
- fix(scrna.SeuratClusterStats): fix when split-by is specified
- feat(scrna.CellsDistribution): support `prefix_each`
- feat(scrna.MarkersFinder): allow set max number of genes to plot in dotplots
- feat(scrna.MarkersFinder): support setting detailed arguments for overlapping plots
- feat(scrna.MarkersFinder): support `prefix_group`
- feat(scrna.ScFGSEA): support `prefix_each`
- feat(scrna.RadarPlots): support `prefix_each` and `subset`
- choir(scrna.SeuratClusterStats): use logger instead of print for log messages
- choir(tcr.TCRClustering): print session info for clustcr script
- choir(scrna.MarkersFinder): flatten toc when no `section` and no `ident-1` specified
- choir: extract case expansion pattern (scrna.CellsDistribution, scrna.MarkersFinder, scrna.MetaMarkers, scrna.RadarPlots, scrna.ScFGSEA, scrna.TopExpressingGenes)
- docs: add more detailed docs for `envs.section`
- tests: add assertion for success of the pipelines
- tests: add tests for utils.misc.R
- tests: add r-presto to env_r.yml
- tests: simplify tests for r functions
- tests: pin scipy to 1.8.0 for clustcr in env_r.yml
- tests: refactor tests for core.filters using unittest
- tests: refactor tests for utils.common_docstrs using unittest
- tests: refactor tests for utils.gene using unittest
- tests: refactor tests for utils.mutate_helpers using unittest
- tests: refactor tests for utils.single_cell using unittest

## 0.26.2

- deps: bump datar-pandas to 0.5.5 to dismiss deprecated warnings

## 0.26.1

- deps: bump pipen to 0.14.3
- deps: pin ggplot2 to 3.4 due to breaking changes of 3.5 for test
- fix(utils.misc.R): replace latin and greek characters with closest ascii chars for `slugify()`
- feat(scrna.TopExpressingGenes): support `subset`
- fix(scrna.CellsDistribution): fix the row order of the heatmaps.
- enh(tcr.CloneResidency): add legend for multiplets in upset plots.
- feat(scrna.SeuratClusterStats): add circos plot for cell composition stats (pwwang/immunopipe#46).

## 0.26.0

- deps: bump pipen to 0.14.1
- deps: bump pipen-report to 0.18.1
- fix(scrna.CellsDistribution): fix multiple `cells_by` columns and speed up plotting
- choir(tcr.CloneResidency): mark singletons in Venn diagrams more clear
- fix(scrna.RadarPlots): fix the order of groups on radar plots
- choir(scrna.RadarPlots): transpose the count/percentage table to save to files
- fix(scrna.MarkersFinder): fix generating report json file when no significant genes found
- choir(scrna.MarkersFinder): Plot maximum 20 genes in dotplots
- choir(scrna.MarkersFinder): Do not convert dashs in case names to dots
- feat(utils.misc): add logger for python (which allows pipen-poplog to populate logs to running log)
- feat: add rnaseq.UnitConversion and tests
- feat: add rnaseq.Simulation to simulate RNAseq data
- feat: add snp module and snp.PlinkSimulation process
- feat: add stats module
- feat: add stats.ChowTest
- feat: add stats.LiquidAssoc
- feat: add stats.DiffCoexpr
- feat: add stats.MetaPvalue

## 0.25.4

- deps: bump datar to  0.15.4 (support pandas 2.2)
- fix(utils.single_cell.R): fix `immdata_from_expanded` missing other data columns
- fix(tcr.Immunarch): fix `mutaters` not working when no `subset` is set
- fix(scrna.CellsDistribution): fix `hm_devpars` not working

## 0.25.3

- fix(scrna.CellTypeAnnotation): keep factor meta data when input and output are RDS for celltypist

## 0.25.2

- fix(scrna_metabolic_landscape.MetabolicPathwayHeterogeneity): fix output directory path is not slugified
- chore(tcr.Immunarch): change case filling log to debug level

## 0.25.1

- scrna.CellTypeAnnotation: leave the meta data as is in celltypist wrapper

## 0.25.0

- deps: bump pipen to 0.13.2
- feat: add scrna.AnnData2Seurat and scrna.Seurat2AnnData
- scrna.MarkersFinder: allow to cache `FindAllMarkers` results
- scrna.CellTypeAnnotation: support celltypist (#111)
- scrna.SeuratSubClustering: add `envs_depth = 1` to replace whole `envs.cases` when new case assigned
- test: add tests for celltypist of `CellTypeAnnotation`, `AnnData2Seurat` and `Seurat2AnnData`

## 0.24.2

- deps: bump pipen-report to 0.17.3
- chore: use internal `slugify` instead of `slugify` library
- scrna.SeuratPreparing: fix displaying filters in report
- scrna.SeuratPreparing: fix logging Seurat procedure arguments
- cellranger: add CellRangerSummary
- cell_ranger: use `Iframe` in report to have loading indicators
- cellranger_pipeline: add `CellRangerCountPipeline` and `CellRangerVdjPipeline`

## 0.24.1

- tcr.Immunarch: update spectratyping output file extension to png

## 0.24.0

- deps: bump up deps by pipen 0.13
- deps: add pipen-poplog to populate job logs to running log
- deps: bump pipen-poplog to 0.0.2
- feat: add utils.caching.R
- cellranger: fix inferring sample name when fastqs from mulitple lanes
- scrna.SeuratClustering/SeuratSubClustering: cache Seurat procedures step by step
- scrna.MetaMarkers: limit log messages to be populated to 15
- scrna.SeuratPreparing: log procedure arguments at debug level
- scrna_metabolic_landscape.MetabolicFeaturesIntraSubset: use logger to log so can be poplutated to running log
- scrna_metabolic_landscape.MetabolicPathwayActivity: use logger to log so can be poplutated to running log
- tcr.ImmunarchLoading: add logs for steps
- tcr.TCRClustering: use logger to log so can be poplutated to running log
- tcr.TESSA: log command at debug level
- tcr.Immunarch: add plot_type for divs to support boxplots
- tcr.TCRClustering: fix log_info not found
- tcr.Immunarch: make poplog_max 999 to pop all job logs to running log
- scrna_metabolic_landscape.MetabolicFeaturesIntraSubset: change log level for groups from warning to info

## 0.23.8

- scrna.SeuratPreparing: log `Seurat` procedure arguments
- scrna.ScFGSEA: add `subset` to filter cells

## 0.23.7

- scrna.SeuratPreparing: update log message for transformation/scaling step
- scrna_metabolic_landscape.MetabolicPathwayHeterogeneity: add utils.gsea script source to support localizeGmtfile

## 0.23.6

- feat: support url for gmtfile wherever GSEA is performed (#113)
- utils.gsea.R: fix file path in gsea.R
- tcr.Immunarch: add error message for empty filtered/subset data in diversity
- scrna.SeuratPreparing: correct description of default assay in docstr
- scrna.SeuratPreparing: run also the normal normalization procedures when `SCTransform` is used (useful for visualization purposes on RNA assay)
- scrna.SeuratClustering: add related issue link to `PrepSCTFindMarkers`
- scrna.ModuleScoreCalculator: document the names added by cell cycle score (pwwang/immunopipe#34)
- scrna.SeuratPreparing: support sample names as `reference` for `IntegrateLayers`

## 0.23.5

- scrna.SeuratClusterStats: fix when `frac` or `frac_ofall` is true and no `group-by` nor `split-by` is specified for `stats`
- core.filters: fix when no enriched items found for report component `enrichr`
- scrna.MarkersFinder: fix when no enriched items found
- scrna.MetaMarkers: fix when no enriched items found
- scrna.TopExpressingGenes: fix when no enriched items found
- utils.gsea.R: fix when no enriched items found for `runEnrichr`
- scrna_metabolic_landscript: fix adding report when ncores > 1

## 0.23.4

- scrna.TopExpressingGenes: fix colnames while pulling average expression
- scrna.CellsDistribution: fix when cells_by has multiple column names
- scrna.CellTypeAnnotation: fix the order of the clusters for direct method
- scrna.SeuratClusterStats: add `position` options for bar plots for stats
- scrna.RadarPlots: add `colors` to set the colors of the loops in radar and bar plots
- tcr.Immunarch: add `split_by` and `split_order` to put subplots together in one single plots

## 0.23.3

- tcr.ImmunarchLoading: change mode from `single` to `paired` by default

## 0.23.2

- scrna.RadarPlots: fix test error when not enough observations
- scrna.RadarPlots: add `n` and `mean` to test table

## 0.23.1

- scrna.RadarPlots: fix error when generating report for tests when breakdown is not provided

## 0.23.0

- deps: bump pipen to 0.12.5
- deps: bump pipen-report to 0.16.3
- deps: Update seurat to 5.0.1 in test env file
- chore: Add `/tmp` to .gitignore
- scrna.MarkersFinder: Add `envs.use_presto` to use presto to speed up finding markers
- scrna.MarkersFinder: Fix a bug when subsetting cells
- scrna.MarkersFinder: Set `envs.dbs` to `KEGG_2021_Human` and `MSigDB_Hallmark_2020` by default
- scrna.MarkersFinder: Fix FindAllMarkers/FindMarkers for SCTransform'ed data
- scrna.SeuratPreparing: Fix handling of empty path in `RNAData`
- scrna.SeuratPreparing: Set `envs.gene_qc.min_cells` to 0 by default (instead of 3)
- scrna.SeuratPreparing: Add sample integration procedures
- scrna.SeuratPreparing: Allow to filter genes directly
- scrna.SeuratClustering: Add options to limit string and numeric output length to have more exact caching signature
- scrna.SeuratClustering: Set default `random.seed` to `8525` for `FindClusters`
- scrna.SeuratClustering: Allow multiple resolutions for `FindClusters`
- scrna.SeuratClustering: Print table of idents in log for found clusters
- scrna.SeuratClustering: Move integration procedues to `SeuratPreparing` and do only clustering
- scrna.SeuratClustering: Update tests
- scrna.SeuratClustering: Make the cluster labels start with "c1" instead of "0"
- scrna.SeuratClustering: Default reduction of `RunUMAP` and `FindNeighbors` to pca
- scrna.SeuratClustering: Fix test
- scrna.SeuratClustering: Print less verbosal log
- scrna.SeuratClusterStats: Add `ngenes` to plot the number of genes expressed
- scrna.SeuratClusterStats: Add barplot for `features` and allow aggregation of features
- scrna.SeuratClusterStats: Fix matching kind for plots of features
- scrna.SeuratClusterStats: Use new umap for plotting feature and dimplots for sub-clustering
- scrna.SeuratClusterStats: Use default assay for plotting of number of genes expressed
- scrna.SeuratClusterStats: Add `envs.mutaters` to mutate meta data
- scrna.SeuratClusterStats: Add histograms to plot number of cells against another variable
- scrna.SeuratClusterStats: Fix reduction for subclustering for dimplots
- scrna.SeuratClusterStats: Subset seurat object for featureplots when ident is subclusters
- scrna.SeuratClusterStats: Fix argument layer not excluded for heatmaps in features
- scrna.SeuratClusterStats: Add `frac_ofall` and `transpose` for `stats` to calculate fraction within group or against all cells, and transpose ident and group, respectively
- scrna.ModuleScoreCalculator: Fix features not being passed to `AddModuleScore` as a list
- scrna.ModuleScoreCalculator: Support calculating diffusion map components
- scrna.SeuratMap2Ref: Rename `envs.alias` to `envs.name
- scrna.SeuratMap2Ref: Set default value of `envs.MappingScore.ndim` to 30
- scrna.SeuratMap2Ref: Add `envs.ncores` for parallelization
- scrna.SeuratMap2Ref: Remove preset MapQuery arguments
- scrna.SeuratMap2Ref: Raise an error when envs.MapQuery.refdata is not provided
- scrna.SeuratMap2Ref: Default `envs.use` to the key of `envs.MapQuery.refdata` with single key
- scrna.SeuratMap2Ref: Use layer instead of slot in docstring (Seurat v5)
- scrna.SeuratMap2Ref: Make sure the column of cluster labels is a factor
- scrna.ScFGSEA: Allow to ignore small group when fgsea fails due to all NAs for pre-ranks
- scrna.ScFGSEA: Use default assay and use layer instead of slot (Seurat v5)
- scrna.TopExpressingGenes: Use default assay of Seurat object and fix column names of average expression (Seurat v5)
- scrna.TopExpressingGenes: Change default enrichment gene sets to `KEGG_2021_Human` and `MSigDB_Hallmark_2020`
- scrna.MetaMarkers: Change default enrichment gene sets to `KEGG_2021_Human` and `MSigDB_Hallmark_2020`
- scrna.MetaMarkers: Give better message when tests for genes fail
- scrna.MetaMarkers: Give error message when not enough cells in case
- scrna.CellsDistribution: Allow to order clusters by `envs.cluster_orderby`
- scrna.CellsDistribution: Add heatmaps
- scrna.SeuratSubClustering: Add process
- scrna_metabolic_landscape: Add `InlineNotification` component to imports for report
- scrna_metabolic_landscape.MetabolicFeatures: Fix when default assay is SCT
- scrna_metabolic_landscape.MetabolicFeaturesIntraSubset: Fix when default assay is SCT
- scrna_metabolic_landscape.MetabolicPathwayActivity: Fix when default assay is SCT
- scrna_metabolic_landscape.MetabolicPathwayActivity: Use default assay of Seurat object
- scrna_metabolic_landscape.MetabolicPathwayHeterogenetiy: Fix when default assay is SCT
- scrna.CellTypeAnnotation: Use layer instead of slot of Seurat object (Seurat v5) for sctype
- tcr.ImmunarchLoading: Allow empty path in TCRData column in input file
- tcr.ImmunarchLoading: Do not hide `envs.mode` anymore in docs
- tcr.CloneResidency: Fix stringifying the subject in case it is a factor
- tcr.CloneResidency: Make `section` works in report
- tcr.Immunarch: Support paired chain data for VJ conjuction plots
- tcr.TESSA: Change `envs.assay` to None to use default assay of Seurat object
- scrna_basic: remove scrna_basic pipeline, use immunopipe instead
- scrna.GeneExpressionInvestigation: Remove deprecated code
- scrna.Write10X: Use layer instead of slot (Seurat v5)
- scrna.ExprImputation: Use default assay of seurat object
- scrna.SeuratTo10X: Rename `Write10X` to `SeuratTo10X`
- scrna.SeuratSubClustering: Fix original reduction being poluted by subclustering
- scrna.SeuratClusterStats: Add `avgheatmap` to plot more elegant heatmap for average gene expressions
- scrna.SeuratClusterStats: Fix ident not working for dimplots
- scrna.SeuratClusterStats: Fix for hists when x is a factor/character vector
- scrna.SeuratClusterStats: Add cluster_orderby to order clusters for features
- scrna.SeuratClusterStats: Add na_group to keep NA values in group-by
- scrna.SeuratClusterStats: Allow avgheatmap to plot features other than gene expressions
- scrna.SeuratClusterStats: Add mutate_helpers.R source file
- scrna.SeuratClusterStats: Fix data binding for avgheatmap in features
- utils.mutate_helpers: Change arguments id_col and compare_col to id and compare, respectively
- utils.mutate_helpers: Fix that subset can't be an expression for expanded family
- utils.mutate_helpers: Add top to select top entities (e.g clones)
- scrna.RadarPlots: Add `breakdown` and `test` to break down the cell distribution and run statistic test on the fractions

## 0.22.8

- scrna_metabolic_landscape.MetabolicPathwayActivity: Fix `useNames = NA` being deprecated in matrixStats v1.2 (more locations)
- scrna_metabolic_landscape.MetabolicPathwayActivity: Fix heatmap `column_split`
- scrna_metabolic_landscape.MetabolicFeaturesIntraSubset: Sort groups when being processed
- utils.gsea: Fix `useNames = NA` in rowSds for matrixStats v1.2
- utils.mutate_helpers: Fix tests

## 0.22.7

- scrna_metabolic_landscape.MetabolicPathwayActivity: Fix `useNames = NA` being deprecated in matrixStats v1.2

## 0.22.6

- deps: Bump pipen-board to 0.13.10 (pipen-report to 0.16.2)

## 0.22.5

- docs: Bump pipen-board to 0.13.9 (pipen-report to 0.16.1)
- cellranger.CellRangerCount: Update iframe height in report
- cellranger.CellRangerVdj: Update iframe height in report

## 0.22.4

- utils.mutate_helpers: Update docs

## 0.22.3

- utils.mutate_helpers: Return ids only when subset is true and group is not NA for  `uniq = TRUE` in expanded, collapsed, emerged and vanished

## 0.22.2

- docs: Update logo and favicon
- docs: Update logo height in README.md
- core.filters: Add `exclude` argument to dict_to_cli_args filter
- cellranger: Add CellRangerCount and CellRangerVdj
- scrna.CellTypeAnnotation: Allow using NA to exclude clusters from output Seurat object
- scrna.SeuratClusterStats: Fix path of expression table file
- scrna.MarkersFinder: Use `FindAllMarkers` if `ident.1` is not specified
- scrna.CellsDistribution: Don't add rownames to the output table file
- utils.mutate_helpers: Add `debug` and `each` to expanded, collapsed, vanished and emerged

## 0.22.1

- scrna.CellsDistribution: Export table with distinct columns
- scrna.SeuratMetadataMutater: Warn about existing columns in seurat object
- tcr.ImmunarchLoading: Change `metacols` to `extracols` so essential columns get exported
- tcr.Attach2Seurat: Detach prefix from template in code
- tcr.CDR3AAPhyschem: Detach prefix from template in code
- tcr.Immunarch: Use `immdata$prefix` as prefix by default
- tcr.TCRClustering: Use `immdata$prefix` as prefix by default
- tcr.TESSA: Allow `in.immdata` to be either an RDS file of immunarch object or a text file of cell-level expanded data

## 0.22.0

- Bump pipen-board to 0.13.8 (pipen-report to 0.16)
- Use `render_job` filter to generate report
- utils: Add biopipen palette
- scrna.SeuratClusterStats: Add subset for dimplots to
- scrna.CellsDistribution: Add descr for cases in report
- scrna.CellsDistribution: Save the table only with  necessary columns
- scrna.MarkersFinder: Add dot plot
- scrna.MetaMarkers: Use logger to log messages
- scrna.SeuratClustering: Use logger to log messages
- scrna.SeuratClustering: Add cache option to cache the clustering results if nothing changed except ncores
- delim.SampleInfo: Fix handling of null `exclude_cols`

## 0.21.2

- tcr.Immunarch: Add V-J junction circos plots
- tcr.Immunarch: Refactor logging statements using `r-logger`

## 0.21.1

- deps: Update pipen-board and pygments versions
- docs: Adopt mkdocs-rtd 0.0.10
- docs: Fix internal reference in API docs
- delim.SampleInfo: Refactor data subset logic in SampleInfo class

## 0.21.0

- tcr.Immunarch: Fix empty groups in diversity plot after filtering
- tcr.Immunarch: Add `in.metafile` to allow other meta info (i.e. seurat clusters) for future subsetting
- tcr.Immunarch: Change `envs.mutaters` now on expanded (cell-level) data
- tcr.Immunarch: Add `subset` for cases to do analysis on a subset of data
- tcr.Immunarch: Add `separate_by` also works on other diversity plots
- tcr.Immunarch: Add `ymin` and `ymax` to align diversity plots by `separate_by`
- tcr.Immunarch: Add `ncol` to specify # columns in the combined plots
- scrna.RadarPlots: Fix `envs.order` not working
- scrna.MarkersFinder: Add `overlap` to find overlapping markers between cases (pwwang/immunopipe#24)
- scrna.MarkersFinder: Add `subset` for each case to subset cells
- scrna.MarkersFinder: Add dot plots for cases
- scrna.CellsDistribution: Allow multiple columns for `cells_by`
- scrna.CellsDistribution: Add `subset` for cases to subset cells
- cnv.AneuploidyScoreSummary: Ignore `.call` suffix to get sample name by default
- cnv.AneuploidyScoreSummary: Fix image path in report while `envs.group_cols` is a string (not an array)
- utils.single_cell.R: Add functions to expand, filter and restore immunarch objects
- utils.common_docstrs: Extract some common docstrings for procs
- utils.misc.R: Use r-logger for logging for R scripts
- utils.mutate_helpers.R: Add include_emerged for `expanded()` and include_vanished for `collapsed()`
- utils.mutate_helpers.R: Fix tests
- tests: Add r-logger to test dependencies

## 0.20.7

- (delim.SampleInfo) Add `distinct` to case to perform stats on distinct records
- (scrna_basic) Fix docker image building

## 0.20.6

- ⬆️ Bump pipen-board to 0.13.4
- ✨ [scrna.MarkersFinder] Allow to set assay for `Seurat::FindMarkers()`
- ✨ [scrna.CellsDistribution] Add venn/upset plot for overlapping cell groups in different cases

## 0.20.5

- ⬆️ Bump pipen-board to 0.13.3
- 🏗️ [tcr.CloneResidency] Rename `envs.sample_groups` to `envs.section` to be consistent with other processes
- 🏗️ [tcr.CloneResidency] Allow additional metadata from input for more flexible case definition (i.e. analysis for specific seurat clusters)
- 📝 [scrna.ScFGSEA] Remove the link in the summary of the docstring (since they are not transformed in the report)
- 🎨 [tcr.CDR3AAPhyschem] Give better error message when wrong group items are given

## 0.20.4

- 🐛 [scrna.SeuratClusterStats] Fix toc not saved correct (causing report not including the right sections)

## 0.20.3

- 🐛 [scrna.SeuratPreparing] Fix when cell_qc is None
- 🎨 [scrna.MarkersFinder] Add margins to volcano plot
- 🐛 [scrna.SeuratClusterStats] Fix `ident` in cases of `envs.dimplots` not working

## 0.20.2

- 🚑 [scrna.SeuratPreparing] Fix % in docstring to crash the pipeline

## 0.20.1

- 📝 [scrna_basic/scrna_metabolic_landscape/scrna/tcr] Update docstring
- 🎨 [scrna.MarkersFinder] Try include more genes in volcano plot (pwwang/immunopipe#17)
- 🎨 [scrna.CellsDistribution] Give better error message in CellsDistribution if group value not found (pwwang/immunopipe#16)
- 🚚 [tcr.TCRClusterStats] Rename TCRClusteringStats to TCRClusterStats (pwwang/immunopipe#15)

## 0.20.0

- ⬆️ Bump pipen to 0.12

## 0.19.0

- ⬆️ Bump pipen-report 0.13.1 (pwwang/immunopipe#9, 2)
- ⬆️ Bump pipen-board to 0.12.5
- 💄 [docs] Hide unnecessary items in nav bar
- 💄 [docs] Get docs, especially API docs, formatted better
- 🐛 [delim.SampleInfo] Fix order in pie charts
- 🎨 [delim.SampleInfo] Add stricter checker for input file (pwwang/immunopipe#13)
- 🎨 [scrna.SeuratPreparing] Improve QC plots
- 📝 [scrna.SeuratPreparing] Fix type annotation for `envs.features_defaults.ncol` in docstring
- 🐛 [scrna.CellsDistribution] Fix the cluster order in pie charts
- 🐛 [scrna.SeuratClusterStats] Fix the cluster order in pie charts
- 🎨 [scrna.SeuratClusterStats] Indicate the case name in logs when pie is enable for group-by
- ✨ [scrna.SeuratClusterStats] Allow mutiple columns in the file for `envs.features_defaults.features`
- ✨ [scrna.SeuratClustering] Add number of clusters at the end of log
- 🩹 [scrna.ModuleScoreCalculator] Set default assay to RNA in case module scores only caculated using integrated features
- 📝 [tcr.Immunarch] Fix docstring for `envs.div.args`
- 🎨 [tcr.CloneResidency] Allow order to be optional
- 🎨 [tcr.Immunarch] Allow to skip overlap and gene usage analyses by setting method to `none` (pwwang/immunopipe#11, pwwang/immunopipe#12)
- 🐛 [tcr.TCRClusteringStats] Don't cluster on heatmap when there are only 2 samples
- 🐛 [scrna_metabolic_landscape.MetabolicFeatures] Import Seurat explictly to avoid satijalab/seurat#2853
- 🐛 [scrna_metabolic_landscape.MetabolicPathwayActivity] Fix when NA values in data for heatmap
- 🐛 [scrna_metabolic_landscape.MetabolicPathwayHeterogeneity] Fix error when no significant pathways selected

## 0.18.3

- 🐛 [scrna.MarkersFinder] Fix when either ident is empty

## 0.18.2

- 🐛 [tcr.CDR3AAphyschem.R] Fix a bug when the min length of CDR3 seqs > 12

## 0.18.1

- ⬆️ Bump datar to 0.15.3
- 🎨 [scrna.MetaMarkers/ScFGSEA/SeuratClusterStats] Remove `tidyseurat::` prefix for `filter`
- ✨ [tcr.TESSA] Allow the results to be saved to seurat object
- 📝 [tcr.TESSA] Fix docs about envs.assay

## 0.18.0

- 🔧 Update .gitignore
- ⬆️ Bump pipen to 0.11
- ⬆️ Bump datar to 0.15.2
- 🚨 Make line length to 88 for linting
- ✨ [core.filters] Add `skip` argument to `r()`
- 🚑 [tcr.TESSA] Fix type annotation for envs.max_iter
- 🐛 [delim.SampleInfo] Allow `unique:` prefix for `on` in stats cases;  fix sample order in plots
- ♻️ [scrna.SeuratClusterStats] Redesign envs
- ✨ [scrna.MarkersFinder] Add volcano plot
- ✨ [tcr.TESSA] Add `envs.assay` for seurat object input
- 🐛 [tcr.TESSA] Fix when a V-gene/J-gene is missing
- ✅ [gsea.FGSEA] Fix tests
- 🚸 [scrna.SeuratClustering] Add clear message when `k.weight` is too large for `IntegrateData`⏎

## 0.17.7

- ✅ [tests] Allow pass FORCE=true to run local-only tests
- ✅ [tests] Fix receiving VERBOSE and FORCE in test script
- 🚑 [tcr.ImmunarchLoading] Fix when `Sample` is the only column in meta
- ✨ [tcr.TESSA] Add process and test

## 0.17.6

- 👷 Fix CI for publishing the package
- ⬆️ Bump pipen-board to 0.11.5
- 🚑 [scrna.SeuratClusterStats] Adjust default width and height for plots
- 🚑 [scrna.CellTypeAnnotation] Keep order of clusters after hitype annotation

## 0.17.5

- 👷 Do not run CI build for publish job
- 🎨 [tcr.TCRClustering] Add `TCR_Cluster_Size1` in addition to `TCR_Cluster_Size` in `out.clusterfile` to represent #cells and #CDR3 seqs
- ⬆️ Bump up dependencies

## 0.17.4

- ✨ [tcr.TCRClustering] Add `TCR_Cluster_Size` in `out.clusterfile`
- 💥 [scrna.SeuratClusterStats] Rename `envs.exprs` to `envs.features`

## 0.17.3

- ⬆️ Bump pipen-report to 0.12.8
- 📝 [delim.SampleInfo] Show h1 in report only when stats specified
- 📝 [delim.SampleInfo] Fix parsing excluded_cols in report

## 0.17.2

- 📝 [delim.SampleInfo] Add report template

## 0.17.1

- 🎨 [scrna.CellTypeAnnotation] Change `seurat_clusters_old` to `seurat_clusters_id` to save old seurat_clusters
- 💥 [csv] Rename to `delim`
- 🚚 [csv.BindRows] Rename to `delim.RowsBinder`
- ✨ [utils.mutate_helpers.R] Add `paired()` to identify paired records
- ✨ [delim.SampleInfo] Add process

## 0.17.0

- ⬆️ Bump pipen-board to 0.11.4
- 📝 [docs] Update logo
- 📝 [docs] Add css due to mkdocs-rtd change
- 💥 [core.filters] Default `sortkeys` to `False` for filter `r`
- 🐛 [scrna.ModuleScoreCalculator] Fix aggregation values of programs
- 🐛 [scrna.SeuratClusterStats] Fix typo for default stats
- 🐛 [scrna.ModuleScoreCalculator] Fix name for cell cycle scores
- 🐛 [scrna.CellsDistribution] Fix when `cells_by` or `group_by` is not an identifier
- 🚑 [utils.mutate_helpers.R] Allow accessing metadata using `.`
- ✨ [scrna.ModuleScoreCalculator] Add proc

## 0.16.7

- 🔥 [scrna.SeuratMetadataMutater] Remove unnecessary in.mutaters
- 📝 [docs] Use kmdocs-rtd for documentation
- 📝 [scrna_basic] Fix docs
- 📝 [docs] Fix CI when files in docs/ changes
- 📝 [docs] Fix CI when CI config file changes
- 📝 [scrna_basic] Update docs for processes
- 🔧 [scrna_basic] Update example config file
- 📝 [docs] Add logo and favicon
- 📝 [docs] Fix font-sizes in APIs
- 📝 [docs] Fix logo size in README

## 0.16.6

- 🚑 [scrna] Hotfix for docstring when parsed by argparse help

## 0.16.5

- 💥 [scrna.SeuratMetadataMutater] Move mutaters from in to envs
- 🔥 [scrna.CellsDistribution] Remove unnecessary in.casefile
- 🚑 [scrna.CellTypeAnnotation] Hotfix when envs.hitype_db as a file starts with "hitypedb_"

## 0.16.4

- 🚑 [scrna.CellTypeAnnotation] Hotfix passing `envs.newcol`
- ⬆️ Bump pipen-report to 0.12.7

## 0.16.3

- 📝 [scrna_metabolic_landscape] Update docstring
- ✨ [tcr.CDR3AAPhyschem] Allow envs.subset_cols to be separated by comma
- ✨ [scrna.CellTypeAnnotation] Add `envs.newcol` to keep original idents

## 0.16.2

- 🚨 Add .lintr for R lintr
- ⬆️ Bump pipen-board to 0.11.1
- 💄 [report] Separate enrichr_report
- 💄 [scrna.CellsDistribution] Fix reports
- 💄 [scrna.CellsDistribution] Reorganize report
- 💄 [scrna.MarkersFinder] Reorganize report
- 💄 [scrna.ScFGSEA] Reorganize report
- 💄 [scrna.TopExpressingGenes] Reorganize report
- 🚨 [scrna.TopExpressingGenes] Fix linting issues in script
- 🔧 [scrna.MarkersFinder] Set envs.prefix_each to True by default
- 🔧 [scrna.TopExpressingGenes] Set envs.prefix_each to True by default
- ✨ [scrna.MetaMarkers] Add proc⏎

## 0.16.1

- 🚨 Fix some linting issues
- ⬆️ Bump pipen-board to 0.11
- 🎨 [scrna.CellTypeAnnotation] Rename `seurat_clusters.old` to `seurat_clusters_old` to save the old clusters for sctype
- 🐛 [scrna.CellTypeAnnotation] Fix saving annotated cell type to text file for sccatch
- 🎨 [scrna.CellTypeAnnotation] Save old clustering to `seurat_clusters_old` for sccatch
- 🎨 [scrna.CellTypeAnnotation] Save old clustering to `seurat_clusters_old` for direct method
- 📝 [scrna.CellTypeAnnotation] Fix links in docs for sccatch
- ✨ [scrna.SeuratClusterStats] Allow `envs.exprs.genes` to be genes directly (separated by ",")
- 💄 [docs] Update API doc styles for dark mode
- ✨ [tcr.TCRClustering] Save the souce code of GIANA with this package
- ✨ [tcr.TCRClusteringStats] Allow multiple cases
- 📝 [tcr.ImmunarchLoading] Update docstring
- ✨ [utils] Add mutate_helpers to identify expanded, collapsed, emerged and vanished clones
- 🐛 [utils/misc.R] Fix list_setdefault and list_update when value is NULL
- 🐛 [scrna.TopExpressionGenes] Fix expanding cases
- ✨ [scrna.SeuratClustering] Allow envs.FindIntegrationAnchors.reference to be a string separated by comma
- ✨ [scrna.ScFGSEA] Allow multiple cases
- ✨ [scrna.MarkersFinder] Allow to use mutate_helpers in envs.mutaters
- 🎨 [scrna.CellsDistribution] Redesign envs to support multiple cases
- 💄 [tcr.Immunarch] Fix report generation for rarefraction analysis
- 🔧 [tcr.Immunarch] Change envs to be less error prone
- 💄 [scrna.CellsDistribution] Fix reports
- 💄 [scrna.ScFGSEA] Fix reports
- ✅ [tests] Fix tests

## 0.16.0

- ⬆️ Bump pipen-board to 0.10
- 💄 [docs] Update docs styles
- 🚨 [core/testing] Remove unused importings
- 🎨 [scrna] Rename RNADir to RNAData for input data
- 🐛 [gsea.GSEA] Replace `doc.string` with `doc_string` to avoid over parsing by pipen-args
- 🎨 [tcr.Immunarch] Refactor and split into modules
- 🎨 [scrna.CellTypeAnnotation] Rename CellTypeAnnotate to CellTypeAnnotation and add hitype
- 🎨 [tcr.ImmunarchLoading] Make it compatible with immunarch 0.9
- 🎨 [scrna.MakersFinder] Support multiple cases
- 🎨 [scrna.TopExpressionGenes] Support multiple cases
- 🐛 [scrna.RadarPlots] Fix section and devpars not passed to script
- 🐛 [scrna.SeuratClustering] Fix PCA on each sample
- 🎨 [scrna.ExprImpution] Rename from ExprImpute to ExprImputation
- 👷 [scrna.CellTypeAnnotation] Add r-hitype to env_r.yml for testing
- 🐛 [scrna.CellTypeAnnotation] Fix typos for hitype script
- 🐛 [scrna.CellTypeAnnotation] Fix startsWith in hitype script
- 🎨 [scrna_basic] Rename `ScrnaBasicAnnotate` to `ScrnaBasicAnnotation`
- 📝 [scrna_basic] Update docs
- 🐛 [cnvkit_pipeline] Fix docker image building
- 📝 [cnvkit_pipeline] Fix docs

## 0.15.2

- ⬆️ Bump pipen-board to 0.9.1
- ✨ [scrna.RadarPlots] Add process
- 🎨 [tcr.Immunarch] Separate diversity in script into a different file
- ✨ [scrna.TopExpressingGenes] Add process
- 🎨 [scrna.CellsDistribution] Use a different color palette
- 🎨 [scrna.SeuratClusterStats] Warn about heatmap without downsampling

## 0.15.1

- ⬆️ Bump pipen-board to 0.8.0
- ⬆️ Bump pipen-report to 0.12.5 (to fix the pydantic error)
- 🎨 [tcr.CloneResidency] Add indicators during running
- 🎨 [tcr.CloneResidency] Allow multiple cases add mutaters for metadata
- 🐛 [misc.File2Proc] Check if input file exists
- 🎨 [tcr.Immunarch] Allow cases for trackings and add mutaters for metadata

## 0.15.0

- ⬆️ Bump pipen to 0.10.6
- ⬆️ Bump pipen-board to 0.7.8
- ➖ Retire cmdy at all places (#54)
- ✅ [core.filters] Add run.env to test
- ✅ [core.filters] Add test for `dashify=True`
- 🎨 [scrna.MarkersFinder] Make envs.sigmarkers case wise for scrna.MarkersFinder (#53)

## 0.14.3

- ⬆️ Bump pipen to 0.10.5
- 🔧 [scrna_metabolic_landscape] Make proc group options for process readonly
- 🎨 [scrna_metabolic_landscape.MetabolicFeatures] Add indicators during computation

## 0.14.2

- ⬆️ Bump pipen-board to 0.7.4
- ⬆️ Bump pipen-report to 0.12.3
- ⚡️ Replace `do.call` with `do_call` in R scripts to improve performance
- 🐛 [scrna.CellTypeAnnotate] Fix when no cell types is given for direct annotation
- 🐛 [cnv.AneuploidyScore] Fix when `envs.cn_tranform` is a list of thresholds

## 0.14.1

- ⬆️ Bump pipen-board to 0.7.3
- ⬆️ Bump other dependencies
- 🎨 [scrna] Add type=int for envs.ncores in docstrings
- 🚑 [tcr.CloneResidency] Dismiss warnings from pivot_wider

## 0.14.0

- ⬆️ Bump pipen-board to 0.6.3
- 🔧 Fix make-examples.sh for docker images for pipelines
- 🚑 [scrna_basic] Fix "Issued certificate has expired" in making examples for docker
- ✨ [tcr.CDR3AAphyschem] Add process
- ✨ [cnv.TMADScore] Add TMADScore and TMADScoreSummary
- 🚑 [cnv.TMADScore] Fix wrong `envs.seg_transform` received in script
- 📝 [cnv.TMADScoreSummary] Add report template
- ✨ [cnv.TMADScoreSummary] Support grouping by 2 groups hierarchically
- 💥 [cnv.AneuploidyScore] Change `envs.include_sex` to `envs.excl_chroms` so exclusion of chroms is more flexible
- 🚑 [cnv.AneuploidyScoreSummary] Adjust with of CAA plot based on number of samples
- ✨ [cnv.AneuploidyScoreSummary] Support grouping by 2 groups hierarchically
- ⬆️ Bump pipen-board to 0.6.3

## 0.13.0

- ⬆️ Bump pipen-board to 0.5.8
- ♻️ [scrna_basic] Change detault tag from dev to master for docker image
- 📝 [scrna_basic] Change detault tag from dev to master in docs
- 🔧 [scrna_basic] Change detault tag from dev to master in entry.sh
- 🔧 [scrna_basic] Fix make-examples.sh when running indenpendently
- 🔧 [scrna_basic] Add plugin_opts.report_no_collapse in board.html
- 🚧 [cnvkit_pipeline] Init docker building
- ⚙️ [cnvkit_pipeline] Make examples
- ⚙️ [cnvkit_pipeline] Update example.json for pipen-board
- 🔧 [cnvkit_pipeline] Fix example in docker image
- 📝 [scrna_metabolic_landscape] Update docstrings to adopt pipen-board
- 📝 [utils.misc] Add docstring for run_command
- 🐛 [cnvkit.CNVkitGuessBaits] Use a better way to determine python of `cnvkit.py`

## 0.12.0

- ⬆️ Bump `pipen` to 0.10
- ⬆️ Bump pipen-runinfo to 0.1.1
- ⬆️ Bump pipen-report to 0.12 and pipen-runinfo to 0.2
- ⬆️ Bump pipen-args to 0.10.2
- ⬆️ Bump pipen-board to 0.5.6
- 📝 Use `flag` instead `action=store_true` in docstring
- ✅ [utils.gene] Fix tests
- 🎨 [scrna.SeuratMap2Ref] Add envs.MappingScore
- ✨ [scrna.SeuratMap2Ref] Add report template
- 💄 [scrna.SeuratMap2Ref] Make figures in 2 columns in report
- ✨ [scrna.CellTypeAnnotate] Add ScCATCH for cell type annotation
- 🎨 [scrna.CellTypeAnnotate] Warn when no cell types are given
- 🐛 [cnvkit] Fix when some arguments are `None`
- 📝 [cnvkit_pipeline] Update docstrings to adopt latest pipen-annotate and pipen-board
- 📝 [cnv] Update docstring
- 🚑 [cnv.AneuploidyScoreSummary] Fix when envs.group_col is None but in.metafile is given
- 👷 [scrna_basic] Init docker image building action
- 👷 [scrna_basic] Fix dockhub credentials
- 📝 [scrna_basic] Update docstrings to adopt latest pipen-annotate and pipen-board
- 📝 [scrna_basic] Add documentation
- 🔧 [scrna_basic] Update configuration for docker image building

## 0.11.0

- ⬆️ Bump pipen to 0.9
- ⬆️ Drop support for python3.7
- ➕ Add pipen-board as dependency
- ✨ Add board.toml for pipen-board to run
- 🐛 [cnvkit.CNVkitCoverage] Fix error when generating flat reference
- 🎨 [bed.BedConsensus] Use bedtools genomecov to calculate the consensus regions
- 🐛 [core.filters] Keep list of dict in python as list of list in R
- ✨ [scrna_metabolic_landscape] Allow multiple subsettings for the data
- ✨ [scrna_basic] Initialize the pipeline
- 🐛 [bed.Bed2Vcf] Fix OrderedDiot not found
- 🎨 [cnvkit_pipeline] Import cached_property directly
- 🐛 [scrna.SeuratPerparing] Fix when input contains a single sample
- 🎨 [tests] Use --reuse instead of --former
- 🐛 [vcf.VcfSplitSamples] Fix missing mutations for extract samples
- 🎨 [scrna_metabolic_landscape.MetabolicPathwayHeterogeneity] Add progress indicator
- 🎨 [scrna.SeuratClustering] Allow sample names to be assigned for reference for FindIntegrationAnchors
- 🎨 [scrna_metabolic_landscape.MetabolicPathwayActivity] Add merged heatmaps for subsets
- 🐛 [scrna_metabolic_landscape.MetabolicPathwayIntraSubsets] Fix fetching subsetting_comparison and limit nproc for FGSEA to 1
- 🎨 [scrna_metabolic_landscape.MetabolicPathwayFeatures] Ignore NAs in subsets
- 🎨 [scrna_metabolic_landscape] Adopt pipen-args 0.9.7
- ✨ [scrna.SeuratMap2Ref] Add process
- ➖ [utils] Retire cmdy
- ✨ [bed.BedtoolsMerge] Add process
- 🎨 [core.testing] Use --cache to control of reusing previous run
- 🎨 [csv.BindRows] Allow to add filename
- 📌 [scrna_basic] Adopt pipen-board 0.1.2
- 🐛 [web.Download] Fix when args is Diot
- 🎨 [cnvkit.CNVkitCall] Detach cmdy
- ✨ [bam.BamSplitChroms] Add process
- ✨ [bam.Merge] Add process and test
- 🐛 [core] Fix repr filter in templates for Diot objects
- 🐛 [docs] Add mygene dep for building utils.gene
- ✅ [vcf.TruvariBench] Pin truvari to v3.4.0 for tests

## 0.10.0

- ⬆️ Adopt pipen-report 0.7 for report templates
- ⚡️ Add todot and sortkeys arguments for filter r
- 🐛 Set default lang for processes using bash
- ⚡️ Update docstrings for processes for pipen-cli-config
- ⚡️ [scrna.ExprImpute] Add progress indicators for alra
- 🐛 [scrna.ExprImpute] Set default assay to RNA for rmagic

## 0.9.0

- ⬆️ Bump up pipen to 0.6

## 0.8.0

- 🚀 [vcf.VcfAnno] Add VcfAnno to use vcfanno to annotate VCF files
- ✨ [tcgamaf.Maf2Vcf] Add Variant_Classification and Variant_Type to output vcf
- ✨ [vcf.VcfFix] Allow gziped vcf as input
- 🧹 Remove tests for core pipeline (not needed any more)

## 0.7.1

- ⬆️ Upgrade pipen-filters to 0.2
- 👽️ Adopt pipen-filters 0.2 in reports
- 🔧 Rename `scrna_metabolic` namespace to `scrna_metabolic_landscape` in entry points
- ✨ [scrna.MarkersFinder] Add `each` for cases to run on each value of metadata variable `each`
- ✨ [tcgamaf.Maf2Vcf] Add proc
- ✨ [bcftools.BcftoolsSort] Add proc

## 0.7.0

- 🧑‍💻 [tcr.Immunarch] Allow separating samples for rarefraction analysis
- ✨ [scrna.SeuratClusterStats] Add expression matrix to output
- 🧑‍💻 [tcr.Immunarch] Allow align_x and log scale for rarefraction analysis
- ✨ [cnv.AneuploidyScgitoreSummary] Add heatmaps
- 🧑‍💻 [tcr.Immunarch] Allow separating samples for rarefraction analysis
- ✨ [scrna.SeuratClusterStats] Add expression matrix to output
- 🧑‍💻 [tcr.Immunarch] Allow align_x and log scale for rarefraction analysis
- ✨ [cnv.AneuploidyScoreSummary] Add heatmaps
- 🐛 [cnv.Aneuploidy] Fix when only one arm has signals for a chromosome
- ✨ [cnvkit.CNVkitGuessBaits] Add proc
- ♻️ [cnvkit_pipeline] Refactor and add docs
- 🎨 [cnvkit_pipeline] Use process decorator to define processes
- ✨ [scrna.SeuratClusterStats] Allow groupby other metadata column than Sample in cell stats
- ✨ [scrna.ExprImput] Add ALRA and set as default
- 🎨 [scrna.scrna_metabolic_landscape] Move from scrna_metabolic and use Seurat object directly instead of sce
- 🐛 [scrna.SeuratClustering] Fix when there are fewer cells
- ✨ [scrna.CellTypeAnnotate] Add proc and tests
- ✨ [scrna.SeuratClusterStats] Allow subsetting for cell stats
- ✅ [vcf.Vcf2Bed] Fix test
- ✅ [tests] Add refgenes for testing
- 🐛 [tests] Fix reference preparing
- ✅ [tests] Add sctype db for tests
- ✅ [tests] Try not patch  using lastest poetry
- ✅ [tests] Build test deps and fix tests
- 👷 [tests] Exclude test_scrna_metabolic_landscape from CI
- ⬆️ Upgrade pipen-cli-run to 0.4.1
- ⬆️ Upgrade pipen to 0.3.11

## 0.6.2

- 🎨 [scripts.utils.vcf] Use format keys for samples
- ✨ [vcf.VcfFix] Dedent envs.helpers automatically and allow it to be list of strings
- 🧑‍💻 [tcr.CloneResidency] Add count table and allow grouping samples in the report
- 🧑‍💻 [cnvkit.CNVkitCall] Allow not passing threshold
- 🧑‍💻 [cnvkit.CNVkitCall] Allow setting cutoff to fetch significant genes for enrichment analysis
- 🧑‍💻 [scrna.SeuratPreparing/SeuratClustering] Do QC in SeuratPreparing only and prepare clustering in SeuratClustering
- ✨ [cnvkit_pipeline] Allow customization of colnames in metafile
- 💚 Fix CI (conda-incubator/setup-miniconda#274)

## 0.6.1

- ✨ [cnvkit_pipeline] Allow purity for each sample
- ✨ [tcr.ImmunarchSplitIdents] Add proc
- ✨ [vcf.VcfSplitSamples] Add proc
- 🏗️ [cnvkit.CNVkitCall] Pass purity as input instead of envs
- ✨ [vcf.VcfIntersect] Add proc
- ✨ [vcf.VcfSampleSplits] Add envs.private to keep only private sites for each sample
- 🔧 Fix setup.py file type
- ✅ Fix tests for utils.gene
- 🚨 Ignore template strings in python scripts for pyright

## 0.6.0

- ✨ [cnv] Add AneuploidyScore and AneuploidyScoreSummary
- ✨ [scrna.Write10X] Add Write10X
- ✨ [cnv.AneuploidyScore] Add envs.include_sex
- 🐛 [scrna.SeuratSubset] Fix when envs.groupby is not given
- ✨ [cnvkit.CNVkitHeatmap] Add envs.order for sample order in the heatmap
- ✨ [bam.CNAClinic] Add bam.CNAClinic
- ✨ [bam.CNAClinic] Add report
- ✨ [cnv.AneuploidyScore] Allow a list of thresholds for `envs.cn_transform`
- ✨ [scrna.SeuratSplit] Add scrna.SeuratSplit
- ✏️ [core] Fix typo in core.proc.Pipeline
- 👽️ Refactor pipeline modules with pipen-cli-run 0.3
- 💚 Use mamba in CI

## 0.5.3

- ✨ [scrna.SeuratClusterStats] Allow features to be a file for expression plots
- ✨ [tcr.CloneSizeQQPlot] Add process
- 🩹 [tcr.Immunarch] Fix bad characters in the “Motif Analysis” section in report (#43)

## 0.5.2

- ⬆️ Pump pipen-args to 0.3
- 🩹 [scrna.CellsDistribution] Filter NA `cells.by`

## 0.5.1

- 💚 Fix CI
- 🚨 Add and fix linting
- ⬆️ Pump pipen-report to 0.4.5

## 0.5.0

- ✅ [vcf.VcfFix] Add chrom size fixes
- ✨ [utils.reference] Add bam_index
- 🐛 [bam.CNVpytor] Fix vcf-fix only adds last contig and fix header with snp data
- ✨ [vcf.Vcf2Bed] Add process and test
- 🐛 [bed.BedConsensus] Fix final weighting issue
- 🩹 [All] Use `%>%` instead of `|>` in all R scripts for backward compatibility
- 🐛 [scrna_metabolic] Don't turn "Ident" to "seurat_clusters" for grouping.groupby in config
- 🏗️ [tests] Add prefix "biopipen-" to conda environment names
- ✅ [tests] Enable pipen-report only when necessary

## 0.4.9

- 👷 [test] Reverse immunarch in env_r
- ✨ [bam.CNVpytor] Add filters
- ✨ [cnvkit/cnvkit_pipeline] Add processes and pipeline
- 🐛 [bam.cnvkit] Fix filter direction
- 🚑 [scrna_metabolic] Fix nproc for runFGSEA for MetabolicPathwayHeterogeneity

## 0.4.8

- 🩹 [core] Add default for config.exe.bedtools
- 🩹 [scrna.ScFGSEA] Don't convert sparse matrix to avoid "problem too large" error

## 0.4.7

- 🐛 [scrna.SeuratPreparing] Fix new data preparing when errored

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
