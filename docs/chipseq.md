
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
