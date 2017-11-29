{% raw %}

## pMutSig

### description
	MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.

	For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
	
	See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)

### input
#### `maffile:file`:
 mutation table  
#### `cvgfile:file`:
 coverage table  
#### `cvrfile:file`:
 covariates table  
#### `mutdict:file`:
 mutation_type_dictionary_file  
#### `chrdir:file`:
  chr_files_hg18 or chr_files_hg19   

### output
#### `outdir:dir`:
 The output directory  

### args
#### `mutsig`:
 The path to `run_MutSigCV.sh`, default: 'mutsig'  
#### `mcr`:
 The Matlab MCR path  

## pVcf2Maf

### description
	Convert a snpEff-annotated somatic mutation vcf file (with normal and tumor samples) to [maf](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) file

### input
#### `infile:file`:
 vcf file  

### output
#### `outfile:file`:
 The maf file  

### args
    `vepdata`: The path of vep data. Default: "" (default data dir of vep)
    `vep`: The path of vep excutable. Default: "vep"
    `vcf2maf`: The path of vcf2maf excutable. Default: "vcf2maf.pl"
    `reffile`: The reference fasta file.
    `nthread`: The number of threads used by vep. Default: 1
    `filtervcf`: The filter vcf
#### `params`:
 Other parameters for `vcf2maf.pl`, default: ""  

## pMergeMafs

### description
	Merge MAF files

### input
#### `indir:file`:
 The directory containing MAF files to be merged  

## pMutsig4Plot

### description
	Prepare somatic mutations for  plotting

### input
#### `msdir:file`:
   The mutsig output directory  

### output
#### `outfile:file`:
  The file for plotting  
	```
	#PANEL: Somatic mutations
	#INFO: MT|PI
	#DESC: Mutation type|Putative impact
	# could also be bordercolor, you can have up to 4 shape features
	#TYPE: shape|bgcolor
	# could also be continuous
	# expressions for set: a,b,c
	#                 norminal: no
	#                 continuous: [0,1]
	#DATA: set|norminal
	#NCOL: 2|2
	#NAME_MT: Frameshift|Missense|Nonsense|Silent|Splice_site|TSS|Nonstop
	#NAME_PI: HIGH|MODERATE|LOW|MODIFIER
	#VALUE_MT: 0|1|20|13|4|17|14
	#EXP_MT: frameshift_variant,inframe_deletion,inframe_insertion|missense_variant,initiator_codon_variant,stop_retained_variant,rare_amino_acid_variant|stop_gained|synonymous_variant|splice_acceptor_variant,splice_donor_variant|start_lost,start_retained|stop_lost
	#
	Sample1	Sample2	Sample3	Sample4	Sample5
	ABC	missense_variant|HIGH	missense_variant|HIGH	...
	...
	```

### args
#### `topn`:
     the cutoff to select genes. If it is >= 1, top N genes will be selected, otherwise, it will be used as pvalue cutoff. Default: .05  

## pMutPlot

### description
	Plot mutations
	```
	|           |             |           |           |---
	|- ftWidth -|  s   s   s  |- pnWidth -|- lgWidth -| snHeight
	|           |             |           |           |---
	    feature1
		feature2
	```

### input
#### `indir:file`:
    The input directory containing plot files  

## pCepip

### input
#### `avinput:file`:
 The avinput file  
#### `cell`:
         The cell  

### output
#### `outfile:file`:
 The cepip result file  

### args
#### `bin-cepip`:
    The jar file path of cepip, default: /data2/junwenwang/shared/tools/cepip/cepip.jar  
{% endraw %}
