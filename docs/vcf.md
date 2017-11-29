{% raw %}

## pVcfFilter

### description
	Filter records in vcf file.

### input
#### `infile:file`:
 The input file  

### output
#### `outfile:file`:
 The output file  

### args
#### `filters`:
 The filters, should be a string of lambda function:  
		```
		"lambda record, samples: <expression>"
		* ``record.CHROM`` : 'chr20'
		* ``record.POS``   : 1234567
		* ``record.ID``    : 'microsat1'
		* ``record.REF``   : ''GTC''
		* ``record.ALT``   : [G, GTCT]
		* ``record.QUAL``  : 50
		* ``record.FILTER``: ['PASS']
		* ``record.INFO``  : {'AA': 'G', 'NS': 3, 'DP': 9}
		* samples = record.samples
		* len(samples): 3
		* samples[0].sample: 'NA00001'
		* samples[0]: Call(sample=NA00001, CallData(GT=0/1, GQ=35, DP=4))
		* samples[0].data: calldata(GT='0/1', GQ=35, DP=4)
		* samples[0].data.GT: '0/1'
		```
		- see here for record and samples: https://github.com/jamescasbon/PyVCF
		- Remember if filters() returns True, record remained.
#### `gz`     :
 Whether to gzip the output file. Default: False  
#### `keep`   :
 Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)  

## pVcfAnno

### description
	Annotate the variants in vcf file.
	You have to prepare the databases for each tool.

### input
#### `infile:file`:
 The input vcf file  

### output
#### `outfile:file`:
 The output file (output file of annovar will also be converted to vcf)  
#### `outdir`:
 The output directory, used to fetch some stat/summary files  

### args
#### `tool`:
            The tool used to do annotation. Default: snpeff  
#### `snpeff`:
          The path of snpeff. Default: snpEff  
#### `vep`:
             The path to vep. Default: vep  
#### `gz`:
              Whether to gzip the result file. Default: False  
#### `annovar`:
         The path of annovar. Default: annotate_variation.pl  
#### `annovar_convert`:
 The path of convert2annovar.pl, used to convert vcf to annovar input file. Default: convert2annovar.pl  
#### `genome`:
          The genome for annotation. Default: hg19  
#### `tmpdir`:
          The tmpdir, mainly used by snpeff. Default: <system tmpdir>  
#### `dbpath`:
          The path of database for each tool. Required by 'annovar' and 'vep'  
#### `params`:
          Other params for tool. Default: ''  
#### `snpeffStats`:
     Whether to generate stats file when use snpeff. Default: False  
#### `mem`:
             The memory used by snpeff. Default: '4G'  

## pVcfSplit

### description
	Split multi-sample Vcf to single-sample Vcf files.

### input
#### `infile:file`:
 The input vcf file  
#### `samples`:
     The samples, if not provided, will extract all samples  

### output
#### `outdir:dir`:
  The output directory containing the extracted vcfs  

## pVcfMerge

### description
	Merge single-sample Vcf files to multi-sample Vcf file.

### input
#### `indir:dir`:
 The directory containing multiple vcf files  

### output
#### `outfile:dir`:
  The output multi-sample vcf.  

## pVcf2Maf

### description
	Convert Vcf file to Maf file

### input
#### `infile:file` :
 The input vcf file  

### output
#### `outfile:file`:
 The output maf file  
{% endraw %}
