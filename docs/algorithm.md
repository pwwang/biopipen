<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pRWR](#prwr)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pRWR

### description
	Do random walk with restart (RWR)

### input
#### `Wfile:file`:
 The adjecent matrix  
#### `Efile:file`:
 The start vector  

### output
#### `outfile:file`:
 The output of final probabilities  

### args
#### `c`:
       The restart probability. Default: 0.1  
#### `eps`:
     The convergent cutoff || R(i+1) - R(i) ||. Default: 1e-5  
#### `niter`:
   Max iterations to stop. Default: 10000  
#### `normW`:
   Weather to normalize W or not, default True.   
		- Laplacian normalization is used (more to add).
#### `normE`:
   Weather to normalize E or not, default True.   
		- E will be normalized as: E = E/sum(E)
