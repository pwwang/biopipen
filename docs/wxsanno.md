
## pSnpEff

### description
	This is the default command. It is used for annotating variant filed (e.g. VCF files).

### input
#### `infile:file`:
  The input file   

### output
#### `outdir:file`:
 The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html  

### args
#### `snpEff`:
       The snpEff executable, default: "snpEff"  
#### `params`:
    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"  
#### `genome`:
    The genome used for annotation, default: "hg19"  
#### `informat`:
  The format of input file [vcf or bed], default: "vcf"  
#### `outformat`:
 The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"  
#### `csvStats`:
  Whether to generate csv stats file, default: True.  
#### `htmlStats`:
 Whether to generate the html summary file, default: False.  
#### `javamem`:
   The memory to use. Default: '-Xms1g -Xmx8g'  
