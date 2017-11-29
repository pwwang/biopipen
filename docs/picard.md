
## pMarkDuplicates

### description
	Identifies duplicate reads.

	This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for additional notes on PCR duplication artifacts. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.

	The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).

	The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following [blog post](https://www.broadinstitute.org/gatk/blog?id=7019) for additional information.

	Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate.

	MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.

	The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.

	If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

### input
#### `infile:file`:
  The bam file   

### output
#### `outfile:file`:
 The marked bam file  

### args
#### `picard`:
     The picard executable, default: "picard"  
#### `params`:
  Other parameters for picard MarkDuplicates, default: ""  
#### `tmpdir`:
  The tmpdir to use. Default: /tmp  

## pAddOrReplaceReadGroups

### description
	Replace read groups in a BAM file.This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.

	For more information about read groups, see the [GATK Dictionary entry](https://www.broadinstitute.org/gatk/guide/article?id=6472). 
	
	This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) (see http://ga4gh.org/#/documentation).

### input
#### `infile:file`:
  The bam file  
#### `rg`:
           The read group information. For example:  
		- "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"

### output
#### `outfile:file`:
 The bam file with read group added  

### args
#### `picard`:
     The picard executable, default: "picard "  
#### `params`:
  Other parameters for picard AddOrReplaceReadGroups, default: ""  

## pCreateSequenceDictionary

### description
	Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records.

	The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

### input
#### `infile:file`:
  The fasta file   

### output
#### `outfile:file`:
 The same fasta file, but with dict file created  

### args
#### `picard`:
     The picard executable, default: "picard"  
#### `params`:
  Other parameters for picard CreateSequenceDictionary, default: ""  

## pCollectWgsMetrics

### description
	Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.

	This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
	
	Note: Metrics labeled as percentages are actually expressed as fractions!

### input
#### `infile:file`:
  The bam file   

### output
#### `outfile:file`:
 The metrics file  

### args
#### `picard`:
     The picard executable, default: "picard"  
#### `params`:
  Other parameters for `picard CollectWgsMetrics`, default: ""  
#### `reffile`:
 The reference file, default: ""  

## pSortSam

### description
	Use `picard SortSam` to sort sam or bam file

### input
#### `infile:file`:
  The sam or bam file to be sorted  

### output
#### `outfile:file`:
 The sorted sam or bam file  

### args
#### `picard`:
     The picard executable, default: "picard"  
#### `order`:
   The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate  
#### `outtype`:
 The type of output file, sam or bam. Default: bam  
#### `params`:
  Other parameters for `picard SortSam`, default: ""  
#### `tmpdir`:
  The tmpdir to use. Default: /tmp  
#### `javamem`:
 The memory for java vm. Default: "-Xms1g -Xmx8g"  

## pIndexBam

### description
	Use `picard BuildBamIndex` to index bam file

### input
#### `infile:file`:
  The bam file   

### output
#### `outfile:file`:
 The same bam file (link) but with .bai file in `proc.outdir`  

### args
#### `picard`:
    The picard executable, default: "picard"  
#### `params`:
  Other parameters for `picard BuildBamIndex`, default: "-Xms1g -Xmx8g"  
