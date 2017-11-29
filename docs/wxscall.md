{% raw %}

## pCNVnator

### description
	Use `CNVnator` to call CNVs from bam file

### input
#### `infile:file`:
  The bam file   

### output
#### `outfile:file`:
 The vcf file  

### args
#### `cnvnator`:
      The CNVnator executable, default: "cnvnator"  
#### `cnv2vcf`:
  The converter executable to convert CNVnator results to vcf, default: "cnvnator2VCF.pl"  
#### `binsize`:
  The bin_size, default: 100  
#### `genome`:
   The genome: default: hg19  
#### `chrom`:
    Chromosome names, default: "" (all chromosomes)  
#### `chrdir`:
   The dir contains reference sequence of chromosomes, default: "" (don't specify)  
	
{% endraw %}
