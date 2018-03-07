# gsea
<!-- toc -->
{% raw %}

## pGMT2Mat

### description
Convert a GMT file to a matrix.
Rownames of GMT file will be the column names of output matrix.

### input
#### `infile:file`:: The input file in GMT format.  

## pExpmat2Gct

### description
Convert expression matrix to GCT file.
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format

### input
#### `expfile:file`:: the input expression matrix file. Samples as columns, genes as rows.  

## pSampleinfo2Cls

### description
Convert sample infomation to cls file.
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for file format
NOTE that the order of samples must be the same as in GMT file in further analysis.

### input
#### `sifile:file`:: the sample information file.  
	- Headers are: [Sample, ]Patient, Group, Batch
	- Rows are samples

## pSSGSEA

### description
Single sample GSEA
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format

### input
#### `gctfile:file`:: the expression file  
#### `gmtfile:file`:: the gmtfile for gene sets  

### output
#### `outdir:file`:: the output directory  
- `report.txt`: the enrichment report for each Gene set.
- `RES_<GeneSet>.png`: the running ES plot for <GeneSet>
- `normP_<GeneSet>.png`: the norminal P value plot for <GeneSet>

## pGSEA

### description
GSEA
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for CLS file format

### input
#### `gctfile:file`:: the expression file  
#### `clsfile:file`:: the class file  
#### `gmtfile:file`:: the gmtfile for gene sets  

### output
#### `outdir:file`:: the output directory  

## pEnrichr

### description
Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list

### input
#### `infile:file`:: The gene list, each per line  

### output
#### `outdir:dir`::  The output directory, containing the tables and figures.  

### args
#### `topn`:: Top N pathways used to plot. Default: 10  
#### `col`:: The columns index containing the genes. Default: 0  
#### `delimit`:: The delimit of input file. Default: '\\t'  
#### `dbs`::  The databases to do enrichment against. Default: KEGG_2016  
  - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
  - Multiple dbs separated by comma (,)
#### `norm`:: Normalize the gene list use [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)  
#### `rmtags`:: Remove pathway tags in the plot. Default: True  
  - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
#### `plot`:: Whether to plot the result. Default: True  
#### `title`:: The title for the plot. Default: "Gene enrichment: {db}"  

## pTargetEnrichr

### description
Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list

### input
#### `infile:file`:: The target genes with regulators  
	- Format: 
	- Header is not required, but may specified in first line starting with `#`
	- If only 3 columns are there, the 3rd column is anyway the relation!
	- If only 4 columns are there, 3rd is target status, 4th is relation!
	  ```
	  #Regulator	Target	Regulator status	Target status	Relation
	  has-mir-22	Gene	+	+	+
	  ```

### output
#### `outdir:dir`::  The output directory, containing the tables and figures.  

### args
#### `dbs`       :: The databases to do enrichment against. Default: KEGG_2016  
  - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
  - Multiple dbs separated by comma (,)
#### `rmtags`    :: Remove pathway tags in the plot. Default: True  
  - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
#### `enrplot`   :: Whether to plot the result. Default: True  
#### `enrn`      :: Top N pathways used to plot. Default: 10  
#### `netplot`   :: Whether to plot the network. Default: True  
#### `netn`      :: Top N pathways used to plot the network. Default: 5  
	- Must <= `enrn`. If `netn` >= `enrn`, `netn` = `enrn`
#### `title`     :: The title for the plot. Default: "Gene enrichment: {db}"  
{% endraw %}
