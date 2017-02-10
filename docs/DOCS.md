
# Documentation for bioprocs v0.0.1
A set of procs for bioinformatics using [pyppl](https://github.com/pwwang/pyppl)
## TCGA
###   
pSample2SubmitterID
#### description
  
convert TCGA sample names with submitter id with metadata and sample containing folder
#### input
  
`dir:file`:    the directory containing the samples`mdfile:file`: the metadata file
#### output
  
`outdir`:      the directory containing submitter-id named files
## DEG
###   
pCallByLimmaFromMatrix
#### description
  
Call DEG from expressoin matrix, where column names must in accordant order of <group>
#### input
  
`matfile:file`: the expression matrix`group1`:       columns of group1 (separated by comma)`group2`:       columns of group2 (separated by comma)`group1name`:   the name of group1`group2name`:   the name of group2
#### output
  
degfile:file: the output file containing DEGs
#### args
  
`pval`: the cutoff of DEGs (default: .05)
#### requires
  
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
###   
pCallByLimmaFromFiles
#### description
  
Call DEG from expression files
#### input
  
`expdir:file`:  the directory containing expression files`group1`:       columns of group1 (separated by comma)`group2`:       columns of group2 (separated by comma)`group1name`:   the name of group1`group2name`:   the name of group2   
#### output
  
`degfile:file`: the output file containing DEGs
#### args
  
`pval`: the cutoff of DEGs (default: .05)
#### requires
  
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
