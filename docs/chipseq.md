<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pPeakToRegPotential](#ppeaktoregpotential)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pPeakToRegPotential

### description
	Convert peaks to regulatory potential score for each gene
	The formula is:
```
		             -(0.5 + 4*di/d0)
	PC = sum (pi * e                  )
```
	Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489297/

### input
#### `peakfile:file`:
 The BED/peak file for peaks  
#### `genefile:file`:
 The BED file for gene coordinates  

### output
#### `outfile:file`:
 The regulatory potential file for each gene  
