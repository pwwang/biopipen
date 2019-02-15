## algorithm

!!! hint "pRWR"

    - **description**  
        Do random walk with restart (RWR)

    - **input**  
        - `Wfile:file`: The adjecent matrix  
        - `Efile:file`: The start vector  

    - **output**  
        - `outfile:file`: The output of final probabilities  

    - **args**  
        - `c`: The restart probability. Default: 0.1  
        - `eps`: The convergent cutoff || R(i+1) - R(i) ||. Default: 1e-5  
        - `niter`: Max iterations to stop. Default: 10000  
        - `normW`: Weather to normalize W or not, default True.   
        	- Laplacian normalization is used (more to add).
        - `normE`: Weather to normalize E or not, default True.   
        	- E will be normalized as: E = E/sum(E)

    - **requires**  
        if normW = True, R package `NetPreProc` is required.

!!! hint "pAR"

    - **description**  
        Affinity Regression.
        Ref: https://www.nature.com/articles/nbt.3343
        ```
                b           c        d          d  
            _________    _______    ____       ____
            |       |    |  W  |    |  |       |  |
          a |   D   |  b |_____|  c |Pt|  =  a |Y |   <=>
            |_______|               |__|       |  |
                                               |__|
        
        kronecker(P, YtD)*vec(W) = vec(YtY)             <=>
        X*vec(W) = vec(YtY)
        WPt:
               c        d              d  
            _______    ____          _____
            |  W  |    |  |          |   |
          b |_____|  c |Pt|  --->  b |___|
                          |__|
        
        YtDW:
        WtDtY:
             b           a        d               d    
          _______    _________   ____           _____  
          |  Wt |    |       |   |  |           |   |  
        c |_____|  b |   Dt  | a |Y |    ---> c |___|  
                     |_______|   |  |                 
                                 |__|                  
        ```

    - **input**  
        - `D:file` : The D matrix  
        - `Pt:file`: The Pt matrix  
        - `Y:file`: The Y matrix  
        	- All input files could be gzipped

    - **output**  
        - `W:file`: The interaction matrix  
        - `outdir:dir`: The output directory  

    - **args**  
        - `seed`: The seed for sampling the training set.  
        - `tfrac`: The fraction of samples used for training.  
## bed

!!! hint "pBedSort"

    - **description**  
        Sort bed files

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `tool`: The tool used to sort the file. Default: sort (bedtools, bedops)  
        - `bedtools`: The path to bedtools. Default: bedtools  
        - `bedops_sort`: The path to bedops' sort-bed. Default: sort-bed  
        - `mem`: The memory to use. Default: 8G  
        - `by`: Sort by coordinates("coord", default) or name("name")  
        	- Only available when use tool `sort`
        - `tmpdir`: The tmpdir to use. Default: `$TMPDIR`  
        - `unique`: Remove the dupliated records? Default: True  
        - `params`: Other params for `tool`. Default: {}  

    - **requires**  
        [`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
        [`bedops`](https://github.com/bedops/bedops)

!!! hint "pBedLiftover"

    - **description**  
        Lift over bed files.

    - **input**  
        - `infile:file`: The input bed file  

    - **output**  
        - `outfile:file`: The output file  
        - `umfile:file` : The unmapped file  

    - **args**  
        - `liftover`: The liftover program  
        - `lochain` : the liftover chain file  

    - **require**  
        `liftover` from UCSC

!!! hint "pGff2Bed"

    - **description**  
        Convert GTF/GFF file to BED file

    - **input**  
        - `infile:file`: The input gtf/gff file  

    - **output**  
        - `outfile:file`: The converted bed file  

    - **args**  
        - `attr2name`: The function used to convert attributes from GTF/GFF file to BED field 'name'  
## bedtools

!!! hint "pBedGetfasta"

    - **description**  
        `bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED file.

    - **input**  
        - `infile:file`: The input bed file  

    - **output**  
        - `outfile:file`: The generated fasta file  

    - **args**  
        - `ref`     : The fasta file  
        - `bedtools`: The bedtools executable,                  default: "bedtools"  
        - `params`  : Other parameters for `bedtools getfasta`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedClosest"

    - **description**  
        Similar to intersect, closest searches for overlapping features in A and B. In the event that no feature in B overlaps the current feature in A, closest will report the nearest (that is, least genomic distance from the start or end of A) feature in B. For example, one might want to find which is the closest gene to a significant GWAS polymorphism. Note that closest will report an overlapping feature as the closest that is, it does not restrict to closest non-overlapping feature. The following iconic cheatsheet summarizes the funcitonality available through the various optyions provided by the closest tool.

    - **input**  
        - `afile:file`: The -a file  
        - `bfile:file`: The -b file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools closest`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedClosest2"

    - **description**  
        Multiple b-file version of pBedClosest

    - **input**  
        - `afile:file`: The -a file  
        - `bfiles:files`: The -b files  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools closest`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedFlank"

    - **description**  
        `bedtools flank` will create two new flanking intervals for each interval in a BED file. Note that flank will restrict the created flanking intervals to the size of the chromosome (i.e. no start < 0 and no end > chromosome size).

    - **input**  
        - `infile:file`: The input file  
        - `gfile:file`: The genome size file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools flank`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedIntersect"

    - **description**  
        By far, the most common question asked of two sets of genomic features is whether or not any of the features in the two sets overlap with one another. This is known as feature intersection. bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported. bedtools intersect works with both BED/GFF/VCF and BAM files as input.

    - **input**  
        - `afile:file` : The a file  
        - `bfile:file`: The b file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools intersect`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedIntersect2"

    - **description**  
        Multiple b-file version of pBedIntersect

    - **input**  
        - `afile:file` : The a file  
        - `bfiles:files`: The b files  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools intersect`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedMakewindows"

    - **description**  
        Makes adjacent or sliding windows across a genome or BED file.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `informat`: The format of input file, whether is a "bed" file or "genome" size file. Default: "bed"  
        - `params`: Other parameters for `bedtools makewindows`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedMerge"

    - **description**  
        `bedtools merge` combines overlapping or book-ended features in an interval file into a single feature which spans all of the combined features.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable,               default: "bedtools"  
        - `params`  : Other parameters for `bedtools merge`, default: {}  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedMerge2"

    - **description**  
        A multi-input file model of pBedMerge: Merge multiple input files.

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable,               default: "bedtools"  
        - `params`  : Other parameters for `bedtools merge`, default: {}  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedMultiinter"

    - **description**  
        Identifies common intervals among multiple BED/GFF/VCF files.

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools multiinter`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedRandom"

    - **description**  
        `bedtools random` will generate a random set of intervals in BED6 format. One can specify both the number (-n) and the size (-l) of the intervals that should be generated.

    - **input**  
        - `gfile:file`: The genome size file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable,    default: "bedtools"  
        - `seed`    : The seed for randomization, default: None  
        - `gsize`   : The chromsize file.  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedShift"

    - **description**  
        `bedtools shift` will move each feature in a feature file by a user-defined number of bases. While something like this could be done with an awk '{OFS="\t" print $1,$2+<shift>,$3+<shift>}', bedtools shift will restrict the resizing to the size of the chromosome (i.e. no features before 0 or past the chromosome end).

    - **input**  
        - `infile:file`: The input file  
        - `gfile:file`: The genome size file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools shift`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedShuffle"

    - **description**  
        `bedtools shuffle` will randomly permute the genomic locations of a feature file among a genome defined in a genome file. One can also provide an exclusions BED/GFF/VCF file that lists regions where you do not want the permuted features to be placed. For example, one might want to prevent features from being placed in known genome gaps. shuffle is useful as a null basis against which to test the significance of associations of one feature with another.

    - **input**  
        - `infile:file`: The input file  
        - `gfile:file`: The genome size file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools shuffle`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedSubtract"

    - **description**  
        `bedtools subtract` searches for features in B that overlap A. If an overlapping feature is found in B, the overlapping portion is removed from A and the remaining portion of A is reported. If a feature in B overlaps all of a feature in A, the A feature will not be reported.

    - **input**  
        - `afile:file`: The a file  
        - `bfile:file`: The b file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools subtract`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedWindow"

    - **description**  
        Similar to `bedtools intersect`, `window` searches for overlapping features in A and B. However, window adds a specified number (1000, by default) of base pairs upstream and downstream of each feature in A. In effect, this allows features in B that are near features in A to be detected.

    - **input**  
        - `afile:file`: The a file  
        - `bfile:file`: The b file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bin`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools window`, default: ""  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedGenomecov"

    - **description**  
        `bedtools genomecov` computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
        
        NOTE: only bam file input implemented here.

    - **input**  
        - `infile:file`: The bam file  

    - **output**  
        - `outfile:file`: The result file  

    - **args**  
        - `bedtools`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools genomecov`, default: `Box(bg = True)`  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

!!! hint "pBedCluster"

    - **description**  
        Similar to merge, cluster report each set of overlapping or "book-ended" features in an interval file. In contrast to merge, cluster does not flatten the cluster of intervals into a new meta-interval; instead, it assigns an unique cluster ID to each record in each cluster. This is useful for having fine control over how sets of overlapping intervals in a single interval file are combined.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file with cluster id for each record  

    - **args**  
        - `bedtools`: The bedtools executable, default: "bedtools"  
        - `params`: Other parameters for `bedtools cluster`, default: `Box()`  

    - **requires**  
        [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)
## chipseq

!!! hint "pPeakToRegPotential"

    - **description**  
        Convert peaks to regulatory potential score for each gene
        The formula is:
        ```
                         -(0.5 + 4*di/d0)
        PC = sum (pi * e                  )
        ```
        Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489297/

    - **input**  
        - `peakfile:file`: The BED/peak file for peaks  
        - `genefile:file`: The BED file for gene coordinates  

    - **output**  
        - `outfile:file`: The regulatory potential file for each gene  

    - **args**  
        - `signal`: `pi` in the formula. Boolean value, whether use the peak intensity signale or not, default: `True`,  
        - `genefmt`: The format for `genefile`, default: `ucsc+gz`. It could be:  
        	- ucsc or ucsc+gz: typically, you can download from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
        	- bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 4th column required as gene identity.
        - `peakfmt`: The format for `peakfile`, default: `peak`. It could be:  
        	- peak or peak+gz: (either [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) or [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13), the 7th column will be used as intensity
        	- bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 5th column will be used as intensity.
        - `window`: `2 * d0` in the formula. The window where the peaks fall in will be consided, default: `100000`.   
        	```
        			|--------- window ----------|
        			|---- d0 -----|
        			|--- 50K --- TSS --- 50K ---|
        				^ (peak center)
        				|-- di --|
        	```
## cluster

!!! hint "pDist2Feats"

    - **description**  
        Convert a distance matrix to coordinates, using multidimensional scaling.

    - **input**  
        - `infile:file`: The distance matrix, could be a full distance matrix, a triangle matrix or a pair-wise distance file  
        	- full dist matrix (full):
        	```
        		s1	s2	s3
        	s1	0	1	1
        	s2	1	0	1
        	s3	1	1	0
        	```
        	- triangle matrix (upper/lower), could be also lower triangle
        	```
        		s1	s2	s3
        	s1	0	1	1
        	s2		0	1
        	s3			0
        	```
        	- pair-wise (pair): (assuming auto-pair-wise distance = 0, that is: `s1	s1	0`)
        	```
        	s1	s2	1
        	s1	s3	1
        	s2	s3	1
        	```
        	- Both rownames and header are required.

    - **output**  
        - `outfile:file`: The output coordinate file  

    - **args**  
        - `informat`: The format of the input file: full, triangle or pair. Default: full  
        	- Could also be upper, lower, pair
        - `k`: How many dimensions? Default: 2 (R^2)  

!!! hint "pFeats2Dist"

    - **description**  
        Calculate the distance of pairs of instances

    - **input**  
        - `infile:file`: The input file. Rows are instances, columns are features.  

    - **output**  
        - `outfile:file`: The output distance matrix  

    - **args**  
        - `inopts`: the options for input file  
        	- `cnames`: Whether input file has column names
        	- `rnames`: Whether input file has row names
        - `sim`: Output similarity instead of distance? Default: `False`  
        - `na` : Replace NAs with? Default: `0`  
        - `method`: Method used to calculate the distance. See available below. Default: `euclidean`  
        	- euclidean, manhattan, minkowski, chebyshev, sorensen, gower, soergel, kulczynski_d, 
        	- canberra, lorentzian, intersection, non-intersection, wavehedges, czekanowski, motyka, 
        	- kulczynski_s, tanimoto, ruzicka, inner_product, harmonic_mean, cosine, hassebrook, 
        	- jaccard, dice, fidelity, bhattacharyya, hellinger, matusita, squared_chord, 
        	- quared_euclidean, pearson, neyman, squared_chi, prob_symm, divergence, clark, 
        	- additive_symm, kullback-leibler, jeffreys, k_divergence, topsoe, jensen-shannon, 
        	- jensen_difference, taneja, kumar-johnson, avg
        	- see: https://cran.r-project.org/web/packages/philentropy/vignettes/Distances.html

    - **requires**  
        `r-philentropy`

!!! hint "pCluster"

    - **description**  
        Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering

    - **input**  
        - `infile:file`: The input matrix file. Clustering will be performed against rows. If not, set `args.transpose` = True  

    - **output**  
        - `outfile:file`: The output cluster file  
        - `outdir:dir`: The output directory containing the figures  

    - **args**  
        - `transpose`: Transpose the input matrix. Default: False  
        - `cnames`: Whether the input matrix contains header before transposing. Default: False  
        - `rnames`: Which column is the rownames before transposing. Default: 1  
        - `plot`: Whether plot the cluster. Default: True  
        - `minc`: Min number of clusters to test. Default: 2  
        - `maxc`: Min number of clusters to test. Default: 15  
        	- If number of rows (nrows) <= 15, then max = nrows - 1
        - `methods`: The methods to test. Default: "all"  
        	- Could be any of "agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"
        	- Multiple methods could be separated by comma (,), or put in a list
        	- By default, fanny, model and sota will be excluded because fanny causes error and the latter two are slow. You can manually include them if you want.
        	- Improper methods will be automatically excluded by `args.isCount`
        - `isCount`: Whether the data is count data. Corresponding methods will be tested. Default: False  

    - **requires**  
        [`r-optCluster`](https://rdrr.io/cran/optCluster/man/optCluster.html)
        [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)

!!! hint "pMCluster"

    - **description**  
        Use `r-mclust` to do clustering. Current just do simple clustering with the package

    - **input**  
        - `infile:file`: The input a coordinate file  

    - **output**  
        - `outdir:dir`: The output of final results  

    - **args**  
        - `transpose`: Transpose the input matrix? Default: False  
        - `rnames`: The `row.names` for `read.table` to read the input file, default: True.  
        - `cnames`: The `header` argument for `read.table` to read the input file, default: True.  
        - `caption`: The caption for the `fviz_cluster`, default: "CLARA Clustering".  
        - `minc`: The min # clusters to try, default: 2  
        - `maxc`: The max # clusters to try, default: 15  

    - **requires**  
        [`r-mclust`](https://cran.r-project.org/web/packages/mclust/index.html)
        [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)

!!! hint "pAPCluster"

    - **description**  
        Use `r-apcluster` to do clustering. 

    - **input**  
        - `infile:file`: The input a coordinate file  

    - **output**  
        - `outdir:dir`: The output of final results  

    - **args**  
        - `transpose`: Transpose the input matrix? Default: False  
        - `rownames`: The `row.names` for `read.table` to read the input file, default: 1.  
        - `header`: The `header` argument for `read.table` to read the input file, default: True.  
        - `caption`: The caption for the `fviz_cluster`, default: "APClustering".  

    - **requires**  
        [`r-apcluster`](https://cran.r-project.org/web/packages/apcluster/index.html)
        [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)

!!! hint "pHCluster"

    - **description**  
        Do hierarchical clustering.

    - **input**  
        - `infile:file`: The input files with variants as rows, features as columns.  
        	- NOTE: clustering is performed on rows, rownames are the leaf labels.

    - **output**  
        - `outdir:dir`: The result directory, containing:  
        	- `hclust.merge.txt`: including merge and height information
        	- `hclust.order.txt`: including order and labels information
        	- `hclust.png`:       the dendrogram plot

    - **args**  
        - `fast`: whether to use `fastcluster` package or not, default: False  
        - `gg`: whether to use `ggdendro` or not, default: False  
        - `rownames`: The `row.names` for `read.table` to read the input file, default: 1.  
        - `header`: The `header` argument for `read.table` to read the input file, default: True.  
        - `method`: Which method to use for `hclust`. Default: "complete" (use `?hclust` to check all availables)  
        - `rotate`: Which to rotate the plot or not. Default: False  
        - `transpose`: Whether to transpose the matrix before cluster. Default: False  

    - **requires**  
        [`r-fastcluster`](https://cran.r-project.org/web/packages/fastcluster/index.html) if `args.fast` is True
        [`r-ggdendro`](https://cran.r-project.org/web/packages/ggdendro/index.html) if `args.gg` is True
## cnvkit

!!! hint "pCNVkitPrepare"

    - **description**  
        Generate target files for cnvkit, using probably cnvkit's access, target, autobin commands.

    - **input**  
        - `infiles:files`: The bam files. Indexes are necessary.  
        	- Hint: if indexes are not with the input files, you probably need `pCNVkitPrepare.infile = 'origin'`

    - **output**  
        - `target:file`: The (autobinned) target file  
        - `antitarget:file`: The (autobinned) target file  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `exbaits`: The bait file for exome-sequencing data.  
        	- See https://github.com/AstraZeneca-NGS/reference_data/tree/master/hg19/bed
        - `accfile`: Directly use the access file. Default: generating from the reference file.  
        	- See https://github.com/etal/cnvkit/tree/master/data
        - `nthread`: The number of threads to use. Default: 1  
        - `ref`    : The reference genome.  
        - `params` : The extra parameters for cnvkit's `access`, `target` and `autobin` command. Default:   
        	```python
        	Box(
        		target  = Box({'short-name': True, 'split': True}),
        		access  = Box(s = '5000'),
        		autobin = Box()
        	)
        	```

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitCov"

    - **description**  
        Calculate coverage in the given regions from BAM read depths.

    - **input**  
        - `infile:file`: The bam file  
        - `tgfile:file`: The target file  
        - `atgfile:file`: The antitarget file  

    - **output**  
        - `outfile:file`: The output cnn file  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params`: Other parameters for `cnvkit.py coverage`  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitRef"

    - **description**  
        Compile a copy-number reference from the given files or directory (containing normal samples). If given a reference genome (-f option), also calculate the GC content and repeat-masked proportion of each region.

    - **input**  
        - `infiles:files`: The input reference coverage files  

    - **output**  
        - `outfile:file`: The output reference cnn file  

    - **args**  
        - `ref`   : The reference file.  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params`: Other parameters for `cnvkit.py reference`, default: " --no-edge "  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitFlatRef"

    - **description**  
        Generate reference coverage if there are no normal samples.

    - **input**  
        - `tgfile:file`: The target file  
        - `atgfile:file`: The antitarget file  

    - **output**  
        - `outfile:file`: The output reference cnn file  

    - **args**  
        - `ref`   : The reference file.  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `params`: Other parameters for `cnvkit.py reference`, default: `{}`  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitFix"

    - **description**  
        Combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference. Output a table of copy number ratios (.cnr)

    - **input**  
        - `tgfile:file`: The target coverage file  
        - `atgfile:file`: The antitarget coverage file  
        - `rcfile:file`: The reference cnn file  

    - **output**  
        - `outfile:file`: The cnr file  

    - **args**  
        - `nthread`: The number of threads to use. Default: 1  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `params`: Other parameters for `cnvkit.py fix`, default: " --no-edge "  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitSeg"

    - **description**  
        Infer discrete copy number segments from the given coverage table

    - **input**  
        - `infile:file`: The cnr file   

    - **output**  
        - `outfile:file`: The cns file  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params`: Other parameters for `cnvkit.py segment`, default: ""  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitCall"

    - **description**  
        Given segmented log2 ratio estimates (.cns), derive each segment's absolute integer copy number 

    - **input**  
        - `infile:file`: The cns file   

    - **output**  
        - `outfile:file`: The callcns file  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `params`: Other parameters for `cnvkit.py segment`, default: ""  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkitScatter"

    - **description**  
        Generate scatter plot for CNVkit results.

    - **input**  
        - `cnrfile:file`: The cnr file  
        - `cnsfile:file`: The cns file from call  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params`: Other parameters for `cnvkit.py scatter`  
        - `regions`: The regoins to plot. Default: `['']`  
        	- You can have extra specific regions, format:
        	- `chr5:100-50000000:TERT` or `chr7:BRAF,MET` (genes are used to highlight)

!!! hint "pCNVkitDiagram"

    - **description**  
        Generate diagram plot for CNVkit results.

    - **input**  
        - `cnrfile:file`: The cnr file  
        - `cnsfile:file`: The cns file from call  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params`: Other parameters for `cnvkit.py scatter`  

!!! hint "pCNVkitHeatmap"

    - **description**  
        Generate heatmap plot for CNVkit results.

    - **input**  
        - `cnfiles:files`: The cnr or cns files.  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `params`: Other parameters for `cnvkit.py scatter`  
        - `regions`: The regoins to plot. Default: `['']`  
        	- You can have extra specific regions, format:
        	- `chr5:100-50000000` or `chr7` (genes are used to highlight)

!!! hint "pCNVkitReport"

    - **description**  
        Report CNVkit results

    - **input**  
        - `cnrfile:file`: The file containing copy number ratio  
        - `cnsfile:file`: The file containing copy number segment  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `cnvkit` : The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params` : Extra parameters to the commands.  
        	- `breaks`:       Whether to report breakpoints. Default: True
        	- `gainloss`:     Whether to report gainloss. Default: True
        	- `metrics`:      Whether to report metrics. Default: True
        	- `segmetrics`:   Whether to report segmetrics. Default: True

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkit2Vcf"

    - **description**  
        Output vcf file for cnvkit results

    - **input**  
        - `cnsfile:file`: The cns file  

    - **output**  
        - `outfile:file`: The vcf file  

    - **args**  
        - `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'  
        - `nthread`: The number of threads to use. Default: 1  
        - `params`: Other params for `cnvkit.py export`  

    - **requires**  
        [CNVkit](http://cnvkit.readthedocs.io/)

!!! hint "pCNVkit2Theta"

    - **description**  
        Convert the results to THetA2 interval input.

    - **input**  
        - `cnsfile:file`: The cns file  
        - `cnnfile:file`: The reference cnn file or the cnr file for paired Normal sample. Could be empty.  

    - **output**  
        - `outfile:file`: The interval file for THetA2  

    - **args**  
        - `nthread` : Number threads to use. Default: `1`  
        - `cnvkit`  : The executable of cnvkit. Default: `cnvkit.py`  
        - `params`  : Other params for `cnvkit.py export theta`  
## common

!!! hint "pSort"

    - **description**  
        Sort file using linux command `sort`

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `skip`   : First N lines to skip. Default: `0`
        	- `delimit`: The delimit. Default          : `\t`
        	- `comment`: The comment line mark. Default: `#`
        - `case`: Case-sensitivity. Default: True  
        	- If True, will set $LANG as C
        	- Otherwise, $LANG will be set as en_US.UTF-8
        - `mem`    : The buffer size. Default: 4G  
        - `tmpdir` : The tmpdir.  
        - `unique` : Just keep the unique lines. Default: False  
        - `delimit`: The delimit to separate the fields. Default: '\t'  
        - `params` : The arguments used by `sort`  

!!! hint "pFiles2Dir"

    - **description**  
        A helper process to convert a list of files into a directory, so that some processes can take it as input

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outdir:dir`: The output directory  

!!! hint "pFile2Proc"

    - **description**  
        Convert a file to a proc so it can be used as dependent

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

!!! hint "pStr2File"

    - **description**  
        Save string to a file.

    - **input**  
        - `in:var`: The input string.  

    - **output**  
        - `outfile:file`: The output file.  

!!! hint "pHead"

    - **description**  
        Get the top N lines from a file

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `n`: Top n lines. You may use '-n' to skip last n lines.  

!!! hint "pTail"

    - **description**  
        Get the bottom N lines from a file

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `n`: Bottom n lines. You may use '+n' to skip first n lines.  

!!! hint "pPrepend"

    - **description**  
        Prepend a string to a file

    - **input**  
        - `in:var`: The input string.  
        - `infile:file`: The input file.  

    - **output**  
        - `outfile:file`: The output file.  

!!! hint "pUnique"

    - **description**  
        Make the input file with unique rows (at certain column)

    - **input**  
        - `infile:file`: The input file.  

    - **output**  
        - `outfile:file`: The output file.  

    - **args**  
        - `inopts`: The options for input file  
        	- `delimit`: delimit for columns 
        	- `skip`: skip first lines
        	- `comment`: signs for treating lines as comments
        - `outopts`: The output options  
        	- `head`: Output head or not. Default: `False`
        	- `headPrefix`: The prefix for the head
        	- `headDelimit`: The delimit for the head
        	- `headTransform`: The transform function for the head
        	- `delimit`: The delimit for the data.
        - `col`: The column to compare. Default: `*` (all columns)  
        - `sorted`: Whether the input file is sorted. Default: `False`  

!!! hint "pAddHeader"

    - **description**  
        Add the header of 1st file to 2nd file.

    - **input**  
        - `infile1:file`: The first file containing the header.  
        - `infile2:file`: The second file with the body.  

    - **output**  
        - `outfile:file`: The output file with the header from 1st input file, body from 2nd file.  

    - **args**  
        - `n`: The number of header lines.  

!!! hint "pMergeFiles"

    - **description**  
        Merge files.

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `header`: Whether the input files have header. Default: `False`  
        	- If `True`, input files must have the same header line.

!!! hint "pSplitRows"

    - **description**  
        Split a file by rows, specially usefull to split a job into multithreads/multiprocesses.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outdir:dir`: The output directory including the split files  

    - **args**  
        - `skip`: The skip first n lines. Default: `0`  
        - `cnames`: The column names. If True, the column names will be added to each split file. Default: `True`  
        - `n`: Number of files to split. Default: `8`  
## docx

!!! hint "pDocx"

    - **description**  
        Operating .docx file

    - **input**  
        - `infile`: A input docx file or a string that used as Heading 0 of a new file  
        - `codes:files`: A set of files containing python codes operating the docs file.  
        	- Variable `doc` is used to pass across files.

    - **output**  
        - `outfile:file`: The output docx file  

    - **args**  
        - `bcode`: Some extra BEFORE the content is inserted (after heading).  
        - `acode`: Some extra AFTER the content is inserted.  
        - `error`: What to do when error happens. Default: `exit`  
        	- `ignore` to add nothing to the document
        - `section`: Start a new section of this table. Could be a box of values. Default: `Box()`.  
        	- `Box()` means no new section will start.
        	- `type`: The type of the section. Could be one of:
        		- `CONTINUOUS`
        		- `NEW_COLUMN`
        		- `NEW_PAGE` (default if not provided when `args.section.type` is not `None`)
        		- `EVEN_PAGE` and
        		- `ODD_PAGE`
        	- `orient`: The orientation of the section. Default: `PORTRAIT`
        		- Could also be `LANDSCAPE`
        	- `margin`: The margin of the section, in points. Could be:
        		- Single value for all margins. Default: 36
        		- Paired values for top/bottom and left/right margins
        		- 3 values for top, left/right and bottom margins.
        		- 4 values for top, right, bottom and left margins

    - **requires**  
        [`python-docx`](http://python-docx.readthedocs.io/en/latest)

!!! hint "pTable2DocxCode"

    - **description**  
        Convert a table to docx code

    - **input**  
        - `infile`: The table data file  

    - **output**  
        - `outfile:file`: The code file  

    - **args**  
        - `font`   : The font of table content. Default: `Box(family: 'Arial', size: 10)`  
        - `section`: Start a new section of this table. Could be a box of values. Default: `Box()`.  
        	- `Box()` means no new section will start.
        	- `type`: The type of the section. Could be one of:
        		- `CONTINUOUS`
        		- `NEW_COLUMN`
        		- `NEW_PAGE` (default if not provided when `args.section.type` is not `None`)
        		- `EVEN_PAGE` and
        		- `ODD_PAGE`
        	- `orient`: The orientation of the section. Default: `PORTRAIT`
        		- Could also be `LANDSCAPE`
        	- `margin`: The margin of the section, in points. Could be:
        		- Single value for all margins. Default: 36
        		- Paired values for top/bottom and left/right margins
        		- 3 values for top, left/right and bottom margins.
        		- 4 values for top, right, bottom and left margins
        - `style`  : The style of the table. Default: `Light List Accent 1`  
        	- `Table Normal`
        	- `Colorful Grid`, `Colorful Grid Accent 1~6`
        	- `Colorful List`, `Colorful List Accent 1~6`
        	- `Colorful Shading`, `Colorful Shading Accent 1~6`
        	- `Dark List`, `Dark List Accent 1~6`
        	- `Light Grid`, `Light Grid Accent 1~6`
        	- `Light List`, `Light List Accent 1~6`
        	- `Light Shading`, `Light Shading Accent 1~6`
        	- `Medium Grid 1~3`, `Medium Grid 1~3 Accent 1~6`
        	- `Medium List 1~2`, `Medium List 1~2 Accent 1~6`
        	- `Medium Shading 1~2`, `Medium Shading 1~2 Accent 1~6`
        	- `Table Grid`
        - `title`  : The title of the table, could be heading to table caption, Could be:  
        	- A string, defaulted as `Heading 1` (2nd level heading)
        	- A `tuple` or `list` with first element (string) as the content, 2nd element as the heading level
        - `caption`: The caption of the table. Default: `None`  
        - `bcode`  : Some extra BEFORE the table is inserted.  
        - `acode`  : Some extra AFTER the table is inserted.  
        - `inopts` : The options for reading the input file.  

    - **requires**  
        [`python-docx` v0.8.7](http://python-docx.readthedocs.io/en/latest)

!!! hint "pImage2DocxCode"

    - **description**  
        Convert an image to docx code

    - **input**  
        - `infile`: The image  

    - **output**  
        - `outfile:file`: The code file  

    - **args**  
        - `scale`: The scale of the image. Default: `Box()`  
        	- With keys `width` and/or `height` in points.
        	- If only one is set, the other is automatically scaled.
        - `title`  : The title of the table, could be heading to table caption, Could be:  
        	- A string, defaulted as `Heading 1` (2nd level heading)
        	- A `tuple` or `list` with first element (string) as the content, 2nd element as the heading level
        	- Template available. i.e. `("Title: {{i.infile | fn2}}", 2)`
        - `section`: Start a new section of this table. Could be a box of values. Default: `Box()`.  
        	- `Box()` means no new section will start.
        	- `type`: The type of the section. Could be one of:
        		- `CONTINUOUS`
        		- `NEW_COLUMN`
        		- `NEW_PAGE` (default if not provided when `args.section.type` is not `None`)
        		- `EVEN_PAGE` and
        		- `ODD_PAGE`
        	- `orient`: The orientation of the section. Default: `PORTRAIT`
        		- Could also be `LANDSCAPE`
        	- `margin`: The margin of the section, in points. Could be:
        		- Single value for all margins. Default: 36
        		- Paired values for top/bottom and left/right margins
        		- 3 values for top, left/right and bottom margins.
        		- 4 values for top, right, bottom and left margins
        - `legend`: The legend of the image. Default: `None`  
        	- If starts with `Figure X. `, this part will be bolded.
        	- Template available. i.e. `Figure 1. {{i.infile | fn2}} shows good results.`
        - `align `: How the image is aligned if it is not an inline element. Default: `CENTER`  
        	- To set it inline, set `args.align` as `None`, and add something to the paragraph(`para`) in `args.acode`.
        - `bcode` : Some extra BEFORE the table is inserted.  
        	- `doc` is available for the `document` object
        	- `sec` is available for the `section` object if `args.section` has been set.
        - `acode` : Some extra AFTER the table is inserted.  
        	- `doc` and `sec` are still available
        	- `para`: The paragraph where the run of the image is placed
        	- `run`: The run where the image is placed

    - **requires**  
        [`python-docx`](http://python-docx.readthedocs.io/en/latest)
## eqtl

!!! hint "pMatrixeQTL"

    - **description**  
        Call eQTLs using Matrix eQTL

    - **input**  
        - `snpfile:file`: The genotype file, rows are snps and columns are samples  
        - `expfile:file`: The expression file, rows are genes  
        - `covfile:file`: The covariant file, rows are covariants  

    - **output**  
        - `outfile:file`: The matrix eqtl output file  

    - **args**  
        - `model`: The model to use, either modelLINEAR(default) or modelANOVA  
        - `pval` : The pvalue cutoff (if `cisopts.dist` > 0, will be used as pval for trans-eQTL)  
        	- Set this to 0 and `cisopts.cispv` to do ciseqtl analysis only.
        - `fdr`  : Calculate FDR or not (default: True)  
        - `cisopts`: Options for calling cis-, trans-eQTL  
        	- `snppos` : The snp position file (columns are: snp, chr, pos)
        	- `genepos`: The gene position file (columns are: gene, chr, start, end)
        	- `dist`   : The distance to define cis-eQTL. (default: 0 (don't do cis-, trans- calling)
        	- `cispv`  : The pvalue cutoff for cis-eQTL (`pval` will not work). Default: `1e-3` 

    - **requires**  
        [`Matrix-eQTL (R)`](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)		
## fastx

!!! hint "pFastq2Expr"

    - **description**  
        Use Kallisto to get gene expression from pair-end fastq files.

    - **input**  
        - `fqfile1:file`: The fastq file1.  
        - `fqfile2:file`: The fastq file2.  

    - **output**  
        - `outfile:file`: The expression file  
        - `outdir:dir`  : Output direcotry with expression and other output files  

    - **args**  
        - `params`  : Other parameters for `kallisto quant`. Default: `Box()`  
        - `idxfile` : The kallisto index file. Default: `params.kallistoIdx`  
        - `kallisto`: The path to `kallisto`. Default: `params.kallisto`  
        - `nthread` : # threads to use. Default: `1`  

!!! hint "pFastqSim"

    - **description**  
        Simulate reads

    - **input**  
        - `seed`: The seed to generate simulation file  
        	- None: use current timestamp.

    - **output**  
        - `fq1:file`: The first pair read file  
        - `fq2:file`: The second pair read file  

    - **args**  
        - `tool`: The tool used for simulation. Default: wgsim (dwgsim)  
        - `len1`: The length of first pair read. Default: 100  
        - `len2`: The length of second pair read. Default: 100  
        - `num`: The number of read PAIRs. Default: 1000000  
        - `gz`: Whether generate gzipped read file. Default: True  
        - `wgsim`: The path of wgsim. Default: wgsim  
        - `dwgsim`: The path of wgsim. Default: dwgsim  
        - `ref`: The reference genome. Required  
        - `params`: Other params for `tool`. Default: ""  

    - **requires**  
        [`wgsim`](https://github.com/lh3/wgsim)

!!! hint "pFastQC"

    - **description**  
        QC report for fastq file

    - **input**  
        - `fq:file`: The fastq file (also fine with gzipped)  

    - **output**  
        - `outdir:dir`: The output direcotry  

    - **args**  
        - `tool`: The tool used for simulation. Default: fastqc  
        - `fastqc`: The path of fastqc. Default: fastqc  
        - `nthread`: Number of threads to use. Default: 1  
        - `params`: Other params for `tool`. Default: ""  

    - **requires**  
        [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

!!! hint "pFastMC"

    - **description**  
        Multi-QC based on pFastQC

    - **input**  
        - `qcdir:file`: The direcotry containing QC files  

    - **output**  
        - `outdir:dir`: The output direcotry  

    - **args**  
        - `tool`: The tool used for simulation. Default: multiqc  
        - `multiqc`: The path of fastqc. Default: multiqc  
        - `params`: Other params for `tool`. Default: ""  

    - **requires**  
        [`multiqc`](http://multiqc.info/)

!!! hint "pFastqTrim"

    - **description**  
        Trim pair-end FASTQ reads

    - **input**  
        - `fq1:file`: The input fastq file  
        - `fq2:file`: The input fastq file  

    - **output**  
        - `outfq1:file`: The trimmed fastq file  
        - `outfq2:file`: The trimmed fastq file  

    - **args**  
        - `tool`        : The tools used for trimming. Default: trimmomatic (cutadapt|skewer)  
        - `cutadapt`    : The path of seqtk. Default: cutadapt  
        - `skewer`      : The path of fastx toolkit trimmer. Default: skewer  
        - `trimmomatic` : The path of trimmomatic. Default: trimmomatic  
        - `params`      : Other params for `tool`. Default: ""  
        - `nthread`     : Number of threads to be used. Default: 1  
        - Not for cutadapt
        - `gz`          : Whether gzip output files. Default: True  
        - `mem`         : The memory to be used. Default: 4G  
        - Only for trimmomatic
        - `minlen`      : Discard trimmed reads that are shorter than `minlen`. Default: 18  
        - For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
        - `minq`        : Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3  
        - `cut5`        : Remove the 5'end reads if they are below qulity. Default: 3  
        - `cut3`        : Remove the 3'end reads if they are below qulity. Default: 3  
        - Not for skewer
        - `adapter1`    : The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  
        - `adapter2`    : The adapter for pair-end sequence. Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA  

    - **requires**  
        [`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
        [`skewer`](https://github.com/relipmoc/skewer)
        [`trimmomatic`](https://github.com/timflutre/trimmomatic)

!!! hint "pFastqSETrim"

    - **description**  
        Trim single-end FASTQ reads

    - **input**  
        - `fq:file`: The input fastq file  

    - **output**  
        - `outfq:file`: The trimmed fastq file  

    - **args**  
        - `tool`        : The tools used for trimming. Default: trimmomatic (cutadapt|skewer)  
        - `cutadapt`    : The path of seqtk. Default: cutadapt  
        - `skewer`      : The path of fastx toolkit trimmer. Default: skewer  
        - `trimmomatic` : The path of trimmomatic. Default: trimmomatic  
        - `params`      : Other params for `tool`. Default: ""  
        - `nthread`     : Number of threads to be used. Default: 1  
        - Not for cutadapt
        - `gz`          : Whether gzip output files. Default: True  
        - `mem`         : The memory to be used. Default: 4G  
        - Only for trimmomatic
        - `minlen`      : Discard trimmed reads that are shorter than `minlen`. Default: 18  
        - For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
        - `minq`        : Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3  
        - `cut5`        : Remove the 5'end reads if they are below qulity. Default: 3  
        - `cut3`        : Remove the 3'end reads if they are below qulity. Default: 3  
        - Not for skewer
        - `adapter`     : The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  

    - **requires**  
        [`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
        [`skewer`](https://github.com/relipmoc/skewer)
        [`trimmomatic`](https://github.com/timflutre/trimmomatic)

!!! hint "pFastqSE2Sam"

    - **description**  
        Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

    - **args**  
        - `tool`: The tool used for alignment. Default: bwa (bowtie2|ngm)  
        - `bwa`: Path of bwa, default: bwa  
        - `ngm`: Path of ngm, default: ngm  
        - `bowtie2`: Path of bowtie2, default: bowtie2  
        - `rg`: The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}  
        - `id` will be parsed from filename with "_LX_" in it if not given
        - `sm` will be parsed from filename
        - `ref`: Path of reference file  
        - `params`: Other params for tool, default: ''  

!!! hint "pFastq2Sam"

    - **description**  
        Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

    - **args**  
        - `tool`   : The tool used for alignment. Default: bwa (bowtie2, ngm, star)  
        - `bwa`    : Path of bwa, default: bwa  
        - `ngm`    : Path of ngm, default: ngm  
        - `star`   : Path of ngm, default: STAR  
        - `bowtie2`: Path of bowtie2, default: bowtie2  
        - `rg`: The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}  
        - `id` will be parsed from filename with "_LX_" in it if not given
        - `sm` will be parsed from filename
        - `ref`    : Path of reference file  
        - `refgene`: The GTF file for STAR to build index. It's not neccessary if index is already been built. Default: ''  
        - `params` : Other params for tool, default: ''  
## gatk

!!! hint "pRealignerTargetCreator"

    - **description**  
        The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such that mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
        Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect. There are 2 steps to the realignment process:
        - Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
        - Running the realigner over those intervals (see the IndelRealigner tool)
        For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).

    - **input**  
        - `bamfile:file`: The aligned bam file  
        - `reffile`: The reference file  

    - **brings**  
        - `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: A list of target intervals to pass to the IndelRealigner.  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `picard`: The picard executable, default: "picard"  
        - `params`: Other parameters for RealignerTargetCreator, default: ""  
        - `samtools`: The samtools executable, default: "samtools"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pIndelRealigner"

    - **description**  
        The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such at mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
        Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect.
        There are 2 steps to the realignment process:
        - Determining (small) suspicious intervals which are likely in need of realignment (see the RealignerTargetCreator tool)
        - Running the realigner over those intervals (IndelRealigner)
        For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).

    - **input**  
        - `bamfile:file`: The aligned bam file  
        - `intfile:file`: Intervals file output from RealignerTargetCreator  
        - `reffile:file`: The reference file  

    - **brings**  
        - `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: A realigned version of input BAM file.  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `picard`: The picard executable, default: "picard"  
        - `params`: Other parameters for IndelRealigner, default: ""  
        - `samtools`: The samtools executable, default: samtools  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pBaseRecalibrator"

    - **description**  
        Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. These scores are per-base estimates of error emitted by the sequencing machines. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes. This tool performs the first step described above: it builds the model of covariation and produces the recalibration table. It operates only at sites that are not in dbSNP; we assume that all reference mismatches we see are therefore errors and indicative of poor base quality. This tool generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and context). Assuming we are working with a large amount of data, we can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, number of observations, number of mismatches, empirical quality score).

    - **input**  
        - `bamfile:file`: A BAM file containing data that needs to be recalibrated.  
        - `reffile:file`: The reference file  

    - **brings**  
        - `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: A GATKReport file with many tables:  
        	- The list of arguments
        	- The quantized qualities table
        	- The recalibration table by read group
        	- The recalibration table by quality score
        	- The recalibration table for all the optional covariates

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `params`: Other parameters for BaseRecalibrator, default: ""  
        - `knownSites`: The known polymorphic sites to mask out, required  
        - `samtools`: The samtools executable, default: samtools  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pPrintReads"

    - **description**  
        PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
        Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.

    - **input**  
        - `bamfile:file`: A BAM file.  
        - `recaltable:file`: The GATKReport file  
        - `reffile:file`: The reference file  

    - **brings**  
        - `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: A single processed bam file.  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `params`: Other parameters for PrintReads, default: ""  
        - `samtools`: The samtools executable, default: samtools  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pHaplotypeCaller"

    - **description**  
        PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
        Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.

    - **input**  
        - `bamfile:file`: A BAM file.  
        - `reffile:file`: The reference file  

    - **brings**  
        - `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: Either a VCF or gVCF file with raw, unfiltered SNP and indel calls.  

    - **args**  
        - `gatk`    : The gatk executable, default: "gatk"  
        - `params`  : Other parameters for HaplotypeCaller, default: ""  
        - `samtools`: The samtools executable, default: samtools  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  
        - `nthread`: Corresponding to -nct option  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pSelectVariants"

    - **description**  
        Often, a VCF containing many samples and/or variants will need to be subset in order to facilitate certain analyses (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain requirements, displaying just a few samples in a browser like IGV, etc.). SelectVariants can be used for this purpose.
        There are many different options for selecting subsets of variants from a larger callset:
        - Extract one or more samples from a callset based on either a complete sample name or a pattern match.
        - Specify criteria for inclusion that place thresholds on annotation values, e.g. "DP > 1000" (depth of coverage greater than 1000x), "AF < 0.25" (sites with allele frequency less than 0.25). These - criteria are written as "JEXL expressions", which are documented in the article about using JEXL expressions.
        - Provide concordance or discordance tracks in order to include or exclude variants that are also present in other given callsets.
        - Select variants based on criteria like their type (e.g. INDELs only), evidence of mendelian violation, filtering status, allelicity, and so on.
        There are also several options for recording the original values of certain annotations that are recalculated when a subsetting the new callset, trimming alleles, and so on.

    - **input**  
        - `vcffile:file`: A variant call set from which to select a subset.  
        - `reffile:file`: The reference file  

    - **brings**  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: A new VCF file containing the selected subset of variants.  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `params`: Other parameters for SelectVariants, default: ""  
        - `samtools`: The samtools executable, default: samtools  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pVariantFiltration"

    - **description**  
        This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.
        The most common way of specifying filtering criteria is by using JEXL queries. See the article on JEXL expressions in the documentation Guide for detailed information and examples.

    - **input**  
        - `vcffile:file`: A variant call set from which to select a subset.  
        - `reffile:file`: The reference file  

    - **brings**  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: A filtered VCF.  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk -T VariantFiltration"  
        - `params`: Other parameters for VariantFiltration, default: ""  
        - `samtools`: The samtools executable, default: samtools  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pMuTect2"

    - **description**  
        MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect ([Cibulskis et al., 2013](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html)) with the assembly-based machinery of HaplotypeCaller. The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller.
        NOTE: only Tumor/Normal variant calling implemented in bioprocs

    - **input**  
        - `tumor:file`: the tumor bam file  
        - `normal:file`: the normal bam file  
        - `reffile:file`: the reference file  

    - **brings**  
        - `tumor`: `{{tumor | bn}}.bai` the index file of tumor  
        - `normal`: `{{normal | bn}}.bai` the index file of normal  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: The vcf file containing somatic mutations  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `samtools`: The samtools executable, default: samtools  
        - `params`: Other parameters for MuTect2, default: ""  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if index files of input files are not found
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.

!!! hint "pMuTect2Interval"

    - **description**  
        Use interval file model of MuTect2

    - **input**  
        - `tumor:file`: the tumor bam file  
        - `normal:file`: the normal bam file  
        - `reffile:file`: the reference file  

    - **brings**  
        - `tumor`: `{{tumor | bn}}.bai` the index file of tumor  
        - `normal`: `{{normal | bn}}.bai` the index file of normal  
        - `reffile#fai`: `{{reffile | bn}}.fai`  
        - `reffile#dict`: `{{reffile | bn}}.dict`  

    - **output**  
        - `outfile:file`: The vcf file containing somatic mutations  

    - **args**  
        - `gatk`: The gatk executable, default: "gatk"  
        - `samtools`: The samtools executable, default: samtools  
        - `params`: Other parameters for MuTect2, default: ""  
        - `picard`: The picard executable, default: "picard"  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [GATK](https://software.broadinstitute.org/gatk)
        [samtools](http://www.htslib.org/) if index files of input files are not found
        [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.
## gene

!!! hint "pGenePromoters"

    - **description**  
        Alias of `seq.pPromoters`.

!!! hint "pGeneNameNorm"

    - **description**  
        Normalize gene names using MyGeneinfo.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `inopts`: options for reading input file.  
        - `outopts`: options for writing output file.  
        	- `query`: Output the original query column? Default: `False`
        - `notfound`: What if a symbol is not found. Default: ignore  
        	- skip  : skip the record(don't write it to output file)
        	- ignore: use the original name;
        	- error : report error
        - `genecol` : the column index containing the gene names  
        - `frm`     : the original format. Default: 'symbol, alias'  
        - `to`      : the output gene name format. Default: 'symbol'  
        - `genome`  : the genome. Default: 'hg19'  
        - `cachedir`: The cache directory  

!!! hint "pIPI"

    - **description**  
        Convert gene symbol to IPI protein accession and vice versa.
        One gene symbol could map to multiple IPIs, which will be separated by pipe (|)

    - **input**  
        - `infile:file` : The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `notfound`: What if a record is not found: Default: `ignore`  
        	- `skip`  : skip the record(don't write it to output file)
        	- `ignore`: use the original name;
        	- `error` : report error
        - `genecol`: The column index containing the gene/protein record  
        - `ipidb`: The IPI xref database (see http://ftp.ebi.ac.uk/pub/databases/IPI/last_release/current/).  
        - `fromipi`: Whether the input is IPI or genes  
        - `inopts`: The options for input file  
        - `outopts`: The options for output file  

!!! hint "pGeneTss"

    - **description**  
        Get gene TSS in BED format.

    - **input**  
        - `infile:file`: The input file containing genes  

    - **output**  
        - `outfile:file`: The output BED file  

    - **args**  
        - `notfound`: What if the gene is not found. Default: skip.  
        	- error: report error
        - `header`: Whether the input file contains header. Default: False  
        - `skip`: Skip N lines of input file. Default: 0  
        	- This has highest priority of header and comment
        - `comment`: The comment line start sign. Default: #  
        - `delimit`: The delimit of input file if it has multiple column. Default: `\\t`  
        - `col`: The column index contains the genes. Default: 0  
        - `frm`: The format of the genes. Default: `symbol, alias`  
        - `tmpdir`: The tmpdir used to store mygene cache files.  
        - `genome`: In which genome to fetch the coordinates. Default: hg19  

!!! hint "pGeneBody"

    - **description**  
        Get gene body region in BED format

    - **input**  
        - `infile:file`: The input file containing genes  

    - **output**  
        - `outfile:file`: The gene body region  

    - **args**  
        - `notfound`: What if a gene is not found when transfer the gene names to gene symbols  
        	- error: report error
        	- skip (default): skip it
        - `inmeta`: The metadata for input file, mainly to indicate where the GENE column is.  
        - `inopts`: Input options for reading input file.  
        	- skip: number of lines to skip. Default: 0
        	- comment: the starting string for comment lines. Default: #
        	- delimit: The delimit for the input file. Default: '\\t'
        frm: The gene name format in the input file. Default: 'symbol, alias'
        tmpdir: The tmpdir to cache the gene name conversion.
        genome: The genome used to do the conversion.
## genomeplot

!!! hint "pInteractionTrack"

    - **description**  
        Gererate genomic interaction track for Gviz

    - **input**  
        - `name`: The name of the track  
        - `infile:file`: The input file.   
        	- See the `type` argument for `makeGenomicInteractionsFromFile` from `GenomicInteractions` r-package
        - `region`: the region, just chromosome!  

    - **output**  
        - `outfile:file`: The dumped track data  

    - **args**  
        - `intype`: Input file type. Default: auto  
        	- Identified by extension
        	- One of "chiapet.tool", "bed12", "bedpe", "hiclib", "homer", "bam", "two.bams".
        - `params`: The display params  

!!! hint "pGeneTrack"

    - **description**  
        Generate gene track using ucsc data source

    - **input**  
        - `name`: The name of the track  

    - **output**  
        - `outfile:file`: The file to save the track  

    - **args**  
        - `genome`: The genome  
        - `params`: use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args  

    - **requires**  
        [r-Gviz](https://rdrr.io/bioc/Gviz)

!!! hint "pAnnoTrack"

    - **description**  
        The annotation track of Gviz

    - **input**  
        - `name`: the name of the track  
        - `infile:file`: the file for the track. (wig, bigWig or bedGraph, bam, need to be indexed!)  
        - `chrom`: the chrom  

    - **output**  
        - `outfile:file`: the rds file for the track  

    - **args**  
        - `genome`: The genome  
        - `params`: See `displayPars(DataTrack())` for all available display params  

    - **requires**  
        [r-Gviz](https://rdrr.io/bioc/Gviz/man/DataTrack-class.html)

!!! hint "pDataTrack"

    - **description**  
        The data track of Gviz

    - **input**  
        - `name`: the name of the track  
        - `infile:file`: the file for the track. (wig, bigWig or bedGraph, bam, need to be indexed!)  
        - `chrom`: the chrom  

    - **output**  
        - `outfile:file`: the rds file for the track  

    - **args**  
        - `genome`: The genome  
        - `params`: See `displayPars(DataTrack())` for all available display params  

    - **requires**  
        [r-Gviz](https://rdrr.io/bioc/Gviz/man/DataTrack-class.html)

!!! hint "pUcscTrack"

    - **description**  
        Generate track from ucsc

    - **input**  
        - `name`     : the name of the track  
        - `track`    : the UCSC track  
        - `trackType`: the Gviz track  
        - `region`   : the region  

    - **output**  
        - `outfile:file`: the dumped track  

    - **args**  
        - `params`: use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.  

    - **requires**  
        [r-Gviz](https://rdrr.io/bioc/Gviz)

!!! hint "pGenomePlot"

    - **description**  
        plot the genomic features

    - **input**  
        - `trkfiles:files`: the list of track dumped files  
        - `region`: the region, in format of `chr1:1-1000`  
        - `highlight`: the highlight regions, informat of start1-end1; start2-end2; ...  

    - **output**  
        - `outfile:file`: the figure  

    - **args**  
        - `genome`  : The genome  
        - `showIdeo`: Show ideogram track? Default: True  
        - `showAxis`: Show axis? Default: True  
        - `showGenes`: Show geneTrack? Default: True  
        - `params`: The params  
        	- `genneral`:  General params for plotTracks
        	- `geneTrack`: The params for geneTrack

    - **requires**  
        [r-Gviz](https://rdrr.io/bioc/Gviz)
## gsea

!!! hint "pGMT2Mat"

    - **description**  
        Convert a GMT file to a matrix.
        Rownames of GMT file will be the column names of output matrix.

    - **input**  
        - `infile:file`: The input file in GMT format.  

    - **output**  
        - `outfile:file`: output matrix file  

!!! hint "pExprMat2GCT"

    - **description**  
        Convert expression matrix to GCT file.
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format

    - **input**  
        - `expfile:file`: the input expression matrix file. Samples as columns, genes as rows.  

    - **output**  
        - `outfile:file`: the gct file  

!!! hint "pSampleinfo2CLS"

    - **description**  
        Convert sample infomation to cls file.
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for file format
        NOTE that the order of samples must be the same as in GMT file in further analysis.

    - **input**  
        - `sifile:file`: the sample information file.  
        	- Headers are: [Sample, ]Patient, Group, Batch
        	- Rows are samples

    - **output**  
        - `outfile:file`: the cls file  

!!! hint "pSSGSEA"

    - **description**  
        Single sample GSEA
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format

    - **input**  
        - `gctfile:file`: the expression file  
        - `gmtfile:file`: the gmtfile for gene sets  

    - **output**  
        - `outdir:file`: the output directory  
        - `report.txt`: the enrichment report for each Gene set.
        - `RES_<GeneSet>.png`: the running ES plot for <GeneSet>
        - `normP_<GeneSet>.png`: the norminal P value plot for <GeneSet>

    - **args**  
        - `weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75  
        - `nperm`: Number of permutations. Default: 10000  

!!! hint "pGSEA"

    - **description**  
        GSEA
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
        Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for CLS file format

    - **input**  
        - `gctfile:file`: the expression file  
        - `clsfile:file`: the class file  
        - `gmtfile:file`: the gmtfile for gene sets  

    - **output**  
        - `outdir:file`: the output directory  

    - **args**  
        - `weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75  
        - `nperm`: Number of permutations. Default: 1000  

!!! hint "pEnrichr"

    - **description**  
        Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list

    - **input**  
        - `infile:file`: The gene list, each per line  

    - **output**  
        - `outdir:dir`: The output directory, containing the tables and figures.  

    - **args**  
        - `top`: Top N pathways used to plot. Default: 10  
        - `genecol`: The columns index containing the genes. Default: 0  
        - `inopts`: The input options.  
        	- `delimit`: The delimit of input file. Default: `\t`
        	- `skip`:    Skip first N lines. Default: `0`
        	- `comment`: Line comment mark. Default: `#`
        	- Other parameters fit `bioprocs.utils.tsvio.TsvReader`
        - `libs`: The databases to do enrichment against. Default: KEGG_2016  
          - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
          - Multiple dbs separated by comma (,)
        - `plot`: Whether to plot the result. Default: True  
        - `devpars`: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`  

!!! hint "pTargetEnrichr"

    - **description**  
        Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list

    - **input**  
        - `infile:file`: The target genes with regulators  
        	- Format (RegulatorStatus and TargetStatus are optional):
        	  ```
        	  Regulator	Target	RegulatorStatus	TargetStatus	Relation
        	  has-mir-22	Gene	+	+	+
        	  ```

    - **output**  
        - `outdir:dir`: The output directory, containing the tables and figures.  

    - **args**  
        - `dbs`       : The databases to do enrichment against. Default: KEGG_2016  
          - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
          - Multiple dbs separated by comma (,)
        - `rmtags`    : Remove pathway tags in the plot. Default: True  
          - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
        - `enrplot`   : Whether to plot the result. Default: True  
        - `enrn`      : Top N pathways used to plot. Default: 10  
        - `netplot`   : Whether to plot the network. Default: True  
        - `netn`      : Top N pathways used to plot the network. Default: 5  
        	- Must <= `enrn`. If `netn` >= `enrn`, `netn` = `enrn`
        - `title`     : The title for the plot. Default: "Gene enrichment: {db}"  

    - **requires**  
        [`python-mygene`](https://pypi.python.org/pypi/mygene/3.0.0)
        [`graphviz`](https://pypi.python.org/pypi/graphviz)
## hic

!!! hint "pPartners"

    - **description**  
        Find the interaction partners of the regions in input file.

    - **input**  
        - `regfile:file`: The region file for the regions to find partners for.  
        - `intfile:file`: The interaction file  

    - **output**  
        - `outfile:file`: The regions with partners.  

    - **args**  
        - `regtype`: The type of region file. Default: `auto` (tell from file extension)  
        	- Could also be `bed` or `bedx`
        - `inttype`: The type of interaction file. Default: `auto`  
        	- Could also be `bedpe`, `chiapet.tool`, `hiclib` and `bed12`
## marray

!!! hint "pCELDir2Matrix"

    - **description**  
        Convert CEL files to expression matrix
        File names will be used as sample names (colnames)

    - **input**  
        - `indir:file`: the directory containing the CEL files, could be gzipped  
        	- If you have files, then use `pFiles2Dir` first
        - `sifile:File`: the sample infor file, extensions for samples are not necessary.  
        	- So that it can be also used by `pMArrayDEG`

    - **output**  
        - `outfile:file`: the expression matrix file  
        - `outdir:dir`: the directory containing expr file and plots  

    - **args**  
        - `pattern`  : The pattern to filter files. Default `'*'`  
        - `norm`     : The normalization method. Default: rma (mas5)  
        - `transfm`  : The extra tranformer for the expression values after the nomralization.  
        	- `Note the expression values have been done with log`
        - `cdffile`  : The cdffile. Default: ''  
        - `fn2sample`: The transformer to transform file name  

!!! hint "pMArrayDEG"

    - **description**  
        Detect DEGs from microarray data.

    - **input**  
        - `efile:file`: The expression matrix file  
        - `gfile:file`: The sample information file  

    - **output**  
        - `outfile:file`: The file with DEGs  
        - `outdir:dir`  : The directory containing files and plots  

    - **args**  
        - `tool`    : The tool to use. Default: `limma`  
        - `annofile`: The annotation file for the probes. Default: ``  
        	- If not provided, raw probe name will be used.
        - `filter`: The filter of the data. Default: `[0, 0]`  
        	- 1st element: the `rowSums` of the expression matrix
        	- 2nd element: how many samples(columns) have to reach the `rowSums` given by the 1st element
        - `pval`  : The pval cutoff for DEGs. Default: `0.05`  
        	- You may specify which kind of measurement to use
        	- `p:0.05`: using cutoff 0.05 for pvalue
        	- `q:0.05`: using cutoff 0.05 for qvalue(fdr)
        	- If no prefix is used, then it defaults to pvalue
        - `plot`  : What kind of plots to generate.   
        	- `mdsplot` : The MDS plot
        	- `volplot` : The volcano plot
        	- `maplot`  : The MA plot
        	- `heatmap` : The heatmap. (use `args.hmrows` to determine how many genes to plot)
        - `ggs`   : The extra ggplot element for each plot (should be able to concatenate by `+`).  
        	- `maplot`  : for MA plot
        	- `heatmap` : for heatmap. Default: `Box(theme = {'axis.text.y': 'r:element_blank()'})`
        	- `volplot` : for volcano plot
        - `devpars`: The parameters for plotting device. Default: `Box(res = 300, width = 2000, height = 2000)`  

    - **requires**  
        `r-limma`
## misc

!!! hint "pGEP70"

    - **description**  
        Calculate GEP70 scores for multiple mylenoma 70-gene-signatures and do survival plot.
        Or add one gene to it to see the survival plot.
        see: https://www.ncbi.nlm.nih.gov/pubmed/17105813

    - **input**  
        - `exprfile:file`: The input file with expression matrix  
        	- Columns are samples, rows are genes
        	- make sure the expression values are log2-scale normalized
        - `survfile:file`: The survival data file (header is required).  
        	- col1: rownames
        	- col2: the survival time
        	- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
        - `gene`: An extra gene to be added to the plot.  
        	- If not provided, only do the GEP70 plot.

    - **output**  
        - `outdir:file`: The output directory, containing:  
        	- The survival data file.
        	- The GEP70 plot and results
        	- The GEP70 + gene plot and results if `i.gene` is provided.

    - **args**  
        - `gep70`: The GEP70 genes.   
        	- Column 1: genes 
        	- Column 2: how is the regulated (up or down)
        - `inunit`  : The time unit in input file. Default: days  
        - `outunit` : The output unit for plots. Default: days  
        - `params`  : The params for `ggsurvplot`. Default: `Box({'risk.table': True, 'conf.int': True, 'font.legend': 13, 'pval': 'Log-rank p = {pval}'})`  
        	- You may do `ylim.min` to set the min ylim. Or you can set it as 'auto'. Default: 0. 
        - `ggs`     : Extra ggplot2 elements for main plot. `ggs.table` is for the risk table.  
        - `devpars` : The device parameters for png. Default: `{res:300, height:2000, width:2000}`  

!!! hint "pNCBI"

    - **description**  
        The NCBI E-Utils

    - **input**  
        - `term`: The term or the id argument for esearch or efetch  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `prog`   : The program to use, esearch (Default) or efetch  
        - `apikey` : The api key for E-utils  
        	- Without API key, we can only query 3 time in a second
        	- With it, we can do 10.
        - `sleep`  : Sleep sometime after job done. Default: `0.15`  
        	- Because of the limit of # queries/sec, we need to sleep sometime after the job is done
        	- At the same time, we also have to limit # jobs to run at the same time. typically: `pNCBI.forks = 10`
        - `db`     : The database to query. Default: `pubmed`. Available databases:  
        	- annotinfo, assembly, biocollections, bioproject, biosample, biosystems, blastdbinfo, books, 
        	- cdd, clinvar, clone, dbvar, gap, gapplus, gds, gencoll, gene, genome, geoprofiles, grasp, gtr, 
        	- homologene, ipg, medgen, mesh, ncbisearch, nlmcatalog, nuccore, nucest, nucgss, nucleotide, 
        	- omim, orgtrack, pcassay, pccompound, pcsubstance, pmc, popset, probe, protein, proteinclusters, 
        	- pubmed, pubmedhealth, seqannot, snp, sparcle, sra, structure, taxonomy, unigene
        - `joiner` : The delimit to use if the field is a list  
        - `record` : A function to transform the record.  

    - **requires**  
        [python-eutils](https://github.com/biocommons/eutils)
## mlearn

!!! hint "pRegressTrain"

    - **description**  
        Train a regression model

    - **input**  
        - `infile:file`: The input file (regression done on columns)  
        - `fmfile:file`: The file defining the regression formula.  
        	```
        	Case1   lm   Y ~ X1 + X2
        	Case2   glm  Y ~ X1 + X2 + X1:X2
        	```
        	- `Case` column can be omitted, if so, `Case1`, `Case2` will be added automatically.
        	- No head.
        - `casefile:file`: Defining cases:  
        	```
        	R1   Case1
        	R2   Case1
        	... ...
        	Rm   Case1
        	Rm   Case2
        	Rm+1 Case2
        	... ...
        	Rn   Case2
        	```
        	- One row can bed used for multiple cases
        	- No head.

    - **output**  
        - `outfile:file`: The output file with the summary of the models  
        - `outdir:dir`  : The output directory containing the summaries, models and plots  

    - **args**  
        - `plot`   : Plot the model? Default: `False`  
        - `inopts` : Input options for input file. Default: `Box(cnames = True, rnames = True)`  
        - `devpars`: The device parameters for the plot. Default: `Box(res = 300, height = 4000, width = 4000)`  

!!! hint "pRegressPred"

    - **description**  
        Use a trained linear regression model to predict

    - **input**  
        - `infile:file`: The input file   
        - `model:file` : The trained model by `pRegressTrain`.  
        	- It can also be a directory containing models generated by `pRegressTrain`.
        	- We will scan the model with file name: `*.<model>.rds`
        	- Models have to be on the same column.

    - **output**  
        - `outfile:file`: The output file, a copy of infile with Case, predicted value and items in `args.out` added.  
        - `outdir:dir`  : The output directory  

    - **args**  
        - `inopts` : The input options.  
        - `out`    : What to include in the output file  
        	- `prob`: Output the probabilities. Default: `True`
        	- `auc` : Output AUCs. Default: `True` (only when Y column exists in infile)
        	- If `auc` is True, then ROC is anyway plotted.
        - `plot`   : Parameters to plot the ROC. Use `False` to disable. Default: `False`  
        	- Enabled only when Y column exists in infile
        	- Could be a dict (Box) with keys:
        	- `labels`:  show labels or not
        	- `showAUC`: show AUC or not
        	- `combine`: combine all roc in one figure?
        - `ggs`    : The ggs items to be added for the plot. Default:   
        	```Box({
        		'style_roc': {},
        		# show legend at bottom right corner
        		'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]} 
        	})```
        - `devpars`: The device parameters for the plot. Default: `Box(res = 300, height = 2000, width = 2000)`  
        - `inopts` : Options for reading the input file. Default: `Box(cnames = True, rnames = True, delimit = "\t")`  

!!! hint "pLogitRegTrain"

    - **description**  
        Train a linear regression model

    - **input**  
        - `infile:file`: The input file (Last column as Y)  

    - **output**  
        - `outmodel:file`: The output model (RData file)  
        - `outdir:dir`   : The output directory containing model, plots and other files  

    - **args**  
        - `plot`   : Whether plot the glm probability. Default: `True`  
        - `formula`: The formula to perform the regression. Default: `None`.  
        	- If `None`, will use all first N-1 columns as features.
        - `inopts` : The input options.  
        - `yval`   : The type of y values. Default: `categ`  
        	- `categ`  : categorical values
        	- `prob`   : probabilities
        	- `numeric`: numeric values

!!! hint "pLogitRegPredict"

    - **description**  
        Use a trained linear regression model to predict

    - **input**  
        - `infile:file`: The input file   
        - `model:file` : The trained model by `pLogitRegTrain`  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `inopts` : The input options.  
        - `outprob`: Also output probabilities? Default: True  

!!! hint "pRandomForestTrain"

    - **description**  
        Train a random forest model

    - **input**  
        - `infile:file`: The input file (Last column as Y)  

    - **output**  
        - `outmodel:file`: The output model (RData file)  
        - `outdir:dir`   : The output directory containing model, plots and other files  

    - **args**  
        - `plot`   : Whether plot the feature importance. Default: `True`  
        - `formula`: The formula to perform the regression. Default: `None`.  
        	- If `None`, will use all first N-1 columns as features.
        - `inopts` : The input options.  
        - `na`     : Replace NAs with? Default: `0`  
        - `devpars`: The device parameters for the plot. Default: `Box(res = 300, height = 2000, width = 2000)`  

    - **requires**  
        `r-randomForst`

!!! hint "pDecisionTreeTrain"

    - **description**  
        Train a decision tree model

    - **input**  
        - `infile:file`: The input file (Last column as Y)  

    - **output**  
        - `outmodel:file`: The output model (RData file)  
        - `outdir:dir`   : The output directory containing model, plots and other files  

    - **args**  
        - `plot`   : Whether plot the feature importance. Default: `True`  
        - `formula`: The formula to perform the regression. Default: `None`.  
        	- If `None`, will use all first N-1 columns as features.
        - `inopts` : The input options.  
        - `na`     : Replace NAs with? Default: `0`  
        - `devpars`: The device parameters for the plot. Default: `Box(res = 300, height = 2000, width = 2000)`  

    - **requires**  
        `r-rpart`

!!! hint "pCrossValid"

    - **description**  
        Do cross validation on a model using R carent package.

    - **input**  
        - `infile:file`: The input data file.  

    - **output**  
        - `outmodel:file`: The trained model in RDS format  
        - `outdir:dir`   : The output directory containing output model and plots.  

    - **args**  
        - `inopts` : The options to read the input file.  
        - `ctrl`   : Arguments for `trainControl`. See `?trainControl`. Default: `Box(method = '', savePredictions = True, classProbs = True)`  
        - `train`  : Arguments for `train` other than `data` and `trControl`. Default: `Box(form = None, method = '', metric = 'ROC')`  
        	- see `?train`
        - `seed`   : The seed. Default: `None`  
        - `nthread`: # threads to use. Default: `1`  
        - `plots`  : Do types of plots. Default: `['model', 'roc']`  
        	- `varimp` also available
        	- You can also concatenate them using comma (`,`)

    - **requires**  
        `r-caret`
## network
## picard

!!! hint "pMarkDuplicates"

    - **description**  
        Identifies duplicate reads.
        
        This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for additional notes on PCR duplication artifacts. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.
        
        The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).
        
        The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following [blog post](https://www.broadinstitute.org/gatk/blog?id=7019) for additional information.
        
        Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate.
        
        MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.
        
        The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.
        
        If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

    - **input**  
        - `infile:file`: The bam file   

    - **output**  
        - `outfile:file`: The marked bam file  

    - **args**  
        - `picard`: The picard executable, default: "picard"  
        - `params`: Other parameters for picard MarkDuplicates, default: ""  
        - `tmpdir`: The tmpdir to use. Default: /tmp  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/)

!!! hint "pAddOrReplaceReadGroups"

    - **description**  
        Replace read groups in a BAM file.This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
        
        For more information about read groups, see the [GATK Dictionary entry](https://www.broadinstitute.org/gatk/guide/article?id=6472). 
        
        This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) (see http://ga4gh.org/#/documentation).

    - **input**  
        - `infile:file`: The bam file  
        - `rg`: The read group information. For example:  
        	- "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"

    - **output**  
        - `outfile:file`: The bam file with read group added  

    - **args**  
        - `picard`: The picard executable, default: "picard "  
        - `params`: Other parameters for picard AddOrReplaceReadGroups, default: ""  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/)

!!! hint "pCreateSequenceDictionary"

    - **description**  
        Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records.
        
        The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

    - **input**  
        - `infile:file`: The fasta file   

    - **output**  
        - `outfile:file`: The same fasta file, but with dict file created  

    - **args**  
        - `picard`: The picard executable, default: "picard"  
        - `params`: Other parameters for picard CreateSequenceDictionary, default: ""  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/)

!!! hint "pCollectWgsMetrics"

    - **description**  
        Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.
        
        This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
        
        Note: Metrics labeled as percentages are actually expressed as fractions!

    - **input**  
        - `infile:file`: The bam file   

    - **output**  
        - `outfile:file`: The metrics file  

    - **args**  
        - `picard`: The picard executable, default: "picard"  
        - `params`: Other parameters for `picard CollectWgsMetrics`, default: ""  
        - `reffile`: The reference file, default: ""  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/)

!!! hint "pSortSam"

    - **description**  
        Use `picard SortSam` to sort sam or bam file

    - **input**  
        - `infile:file`: The sam or bam file to be sorted  

    - **output**  
        - `outfile:file`: The sorted sam or bam file  

    - **args**  
        - `picard`: The picard executable, default: "picard"  
        - `order`: The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate  
        - `outtype`: The type of output file, sam or bam. Default: bam  
        - `params`: Other parameters for `picard SortSam`, default: ""  
        - `tmpdir`: The tmpdir to use. Default: /tmp  
        - `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"  

    - **requires**  
        [picard](http://broadinstitute.github.io/picard/command-line-overview.html)

!!! hint "pIndexBam"

    - **description**  
        Use `picard BuildBamIndex` to index bam file

    - **input**  
        - `infile:file`: The bam file   

    - **output**  
        - `outfile:file`: The same bam file (link) but with .bai file in `proc.outdir`  

    - **args**  
        - `picard`: The picard executable, default: "picard"  
        - `params`: Other parameters for `picard BuildBamIndex`, default: "-Xms1g -Xmx8g"  

    - **requires**  
        [picard](http://broadinstitute.github.io/picard/command-line-overview.html)
## plink

!!! hint "pPlinkMiss"

    - **description**  
        Find samples and snps with missing calls, calculate the call rates and plot them.

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outdir:dir`: The output directory. Default: `{{i.indir | fn}}.miss`  
        	- `.imiss`: The miss calls for samples
        	- `.lmiss`: The miss calls for snps
        	- `.samplecr.fail`: The samples fail sample call rate cutoff (`args.samplecr`)
        	- `.snpcr.fail`: The SNPs fail snp call rate cutoff (`args.snpcr`)

    - **args**  
        - `plink`: The path to plink.  
        - `samplecr`: The sample call rate cutoff. Default: `.95`  
        - `snpcr`: The SNP call rate cutoff. Default: `.95`  
        - `plot`: Whether plot the distribution of the call rates? Default: `True`  
        - `devpars`: The device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`  

!!! hint "pPlinkSexcheck"

    - **description**  
        Check inconsistency between sex denoted and from genotypes.

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outdir:dir`: The output directory. Default: `{{i.indir | fn}}.sexcheck`  
        	- `.sexcheck`: The orginal sex check report from `plink`
        	- `.sex.fail`: The samples that fail sex check.

    - **args**  
        - `plink`: The path to plink.  

!!! hint "pPlinkHet"

    - **description**  
        Calculate the heterozygosity of each sample.

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outdir:dir`: The output directory. Default: `{{i.indir | fn}}.het`  
        	- `.het`: The heterozygosity file generated by `plink`.
        	- `.het.fail`: The samples fail sample heterozygosity cutoff (`args.cutoff`)

    - **args**  
        - `plink`: The path to plink.  
        - `cutoff`: The sample heterozygosity cutoff. Default: `3` (mean-3SD ~ mean+3SD)  
        - `plot`: Whether plot the distribution of the heterozygosity? Default: `True`  
        - `devpars`: The device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`  

!!! hint "pPlinkHWE"

    - **description**  
        Hardy-Weinberg Equilibrium report and filtering.

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outdir:dir`: The output directory. Default: `{{i.indir | fn}}.hwe`  
        	- `.hwe`: The HWE report by `plink`
        	- `.hardy.fail`: The SNPs fail HWE test

    - **args**  
        - `plink`: The path to plink.  
        - `cutoff`: The HWE p-value cutoff. Default: `1e-5`  
        - `plot`: Whether plot the distribution of the HWE p-values? Default: `True`  
        - `devpars`: The device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`  

!!! hint "pPlinkIBD"

    - **description**  
        Estimate the degree of recent shared ancestry individual pairs,
        the identity by descent (IBD)

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outdir:dir`: The output directory. Default: `{{i.indir | fn}}.ibd`  
        	- `.genome`: The original genome report from `plink`
        	- `.ibd.png`: The heatmap of PI_HAT

    - **args**  
        - `plink`: The path to plink.  
        - `indep`: To give a concrete example: the command above that specifies 50 5 0.2 would a) consider a window of 50 SNPs, b) calculate LD between each pair of SNPs in the window, b) remove one of a pair of SNPs if the LD is greater than 0.5, c) shift the window 5 SNPs forward and repeat the procedure.  
        - `pihat`: The PI_HAT cutoff. Default: 0.1875 (see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5007749/)  
        - `plot` : Whether plot the PI_HAT heatmap? Default: `True`  
        - `devpars`: The device parameters for the plot. Default: `Box(res=300, width=2200, height=1600)`  

!!! hint "pPlinkRemove"

    - **description**  
        Remove failed samples and/or SNPs
        The samples/SNPs to be removed should be generated by one of:
        `pPlinkHet`, `pPlinkHWE`, `pPlinkIBD` or `pPlinkMiss`

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  
        - `pdir:dir` : The output directory from one of the processes listed in description  
        	- It could also be the `.fail` file generated by those processes

    - **output**  
        - `outdir:dir`: The output directory containing the `.bed/.bim/.fam` after filtering.  

    - **args**  
        - `plink`: The path to plink.  

!!! hint "pPlink2Vcf"

    - **description**  
        Convert plink binaries into VCF file.

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outfile:file`: The output vcf file.  

    - **args**  
        - `plink`: The path to plink.  
        - `gz`   : Whether bgzip the output vcf file. Default: `False`  
        - `samid`: What to use as sample ID. Default: `both`  
        	- `both`: use `<FID>_<IID>` as sample id
        	- `fid` : use `<FID>` as sample id
        	- `iid` : use `<IID>` as sample id

!!! hint "pPlink2GTMat"

    - **description**  
        Convert plink binaries into genotype matrix.

    - **input**  
        - `indir:dir`: The input directory containing .bed/.bim/.fam files  

    - **output**  
        - `outfile:file`: The output genotype matrix file.  

    - **args**  
        - `plink`: The path to plink.  
        - `samid`: What to use as sample ID. Default: `both`  
        	- `both`: use `<FID>_<IID>` as sample id
        	- `fid` : use `<FID>` as sample id
        	- `iid` : use `<IID>` as sample id
## plot

!!! hint "pPlot"

    - **description**  
        Use ggplot2 to generate plots

    - **input**  
        - `infile:file`: The input data file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `cnames` : Whether the input file has colnames. Default: True  
        - `rnames` : Whether the input file has rownames. Default: False  
        - `aes`    : The default aes. Default: {'x':1, 'y':2} (corresponding to colnames)  
        - `helper` : Some helper codes to generate `params` and `ggs`  
        - `devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
        - `ggs`    : The extra ggplot elements.  

!!! hint "pScatter"

    - **description**  
        Use ggplot2 geom_point to generate plots

    - **infile**  
        - `infile:file`: The input data file  

    - **outfile**  
        - `outfile:file`: The output file  

    - **args**  
        - `cnames` : Whether the input file has colnames. Default: True  
        - `rnames` : Whether the input file has rownames. Default: False  
        - `x`      : The x aes. Default: 1 (corresponding to colnames)  
        - `y`      : The y aes. Default: 2 (corresponding to colnames)  
        - `helper` : Some helper codes to generate `params` and `ggs`  
        - `devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
        - `params` : The extra params for `geom_point`  
        - `ggs`    : The extra ggplot elements.  

!!! hint "pPoints"

    - **description**  
        Alias for pScatter

!!! hint "pHisto"

    - **description**  
        Use ggplot2 geom_histogram to generate histograms

    - **infile**  
        - `infile:file`: The input data file  

    - **outfile**  
        - `outfile:file`: The output file  

    - **args**  
        - `cnames` : Whether the input file has colnames. Default: True  
        - `rnames` : Whether the input file has rownames. Default: False  
        - `x`      : The x aes. Default: 1 (corresponding to colnames)  
        - `helper` : Some helper codes to generate `params` and `ggs`  
        - `devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
        - `params` : The extra params for `geom_point`  
        - `ggs`    : The extra ggplot elements.  

!!! hint "pFreqpoly"

    - **description**  
        Use ggplot2 geom_freqpoly to generate frequency polygon plot.

    - **infile**  
        - `infile:file`: The input data file  

    - **outfile**  
        - `outfile:file`: The output file  

    - **args**  
        - `cnames` : Whether the input file has colnames. Default: True  
        - `rnames` : Whether the input file has rownames. Default: False  
        - `x`      : The x aes. Default: 1 (corresponding to colnames)  
        - `helper` : Some helper codes to generate `params` and `ggs`  
        - `devpars`: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
        - `params` : The extra params for `geom_point`  
        - `ggs`    : The extra ggplot elements.  

!!! hint "pBoxplot"

    - **description**  
        Generate box plot

    - **input**  
        - `infile:file`: The data file  

    - **output**  
        - `outpng:file`: The output figure  

    - **args**  
        - `inopts` : Input options to read the input file  
        	- `cnames` :   Whether the input file has header. Default: `True`
        	- `rnames` :   Whether the input file has row names. Default: `False`
        	- `delimit`:   The seperator. Defualt: `\\t`
        - `x`      : The `ind` (index) column. Only for `args.stacked = True`. Default: `2`  
        - `y`      : The `values` column. Only for `args.stacked = True`. Default: `1`  
        - `helper` : Some raw codes to help to construct the matrix and arguments.  
        - `stacked`: Whether the input file is stacked  
        	- Stacked file looks like:
        	  ```
        	  values	ind
        	  1.1	col1
        	  1.2	col1
        	  ...
        	  .8	col2
        	  .9	col2
        	  ...
        	  3.2	col3
        	  ...
        	  ```
        	- Unstacked file looks like:
        	  ```
        	  col1	col2	col3
        	  1.1	.8	3.2
        	  1.2	.9	2.2
        	  ```
        - `params`: Other parameters for `geom_boxplot`, default: `Box()`  
        - `ggs`   : Extra ggplot2 statements  

!!! hint "pBar"

    - **description**  
        Generate bar/col plot

    - **input**  
        - `infile:file`: The data file  

    - **output**  
        - `outpng:file`: The output figure  

    - **args**  
        - `inopts` : Input options to read the input file  
        	- `cnames` :   Whether the input file has header. Default: `True`
        	- `rnames` :   Whether the input file has row names. Default: `False`
        	- `delimit`:   The seperator. Defualt: `\\t`
        - `x`      : The `ind` (index) column. Only for `args.stacked = True`. Default: `2`  
        - `y`      : The `values` column. Only for `args.stacked = True`. Default: `1`  
        - `helper` : Some raw codes to help to construct the matrix and arguments.  
        - `stacked`: Whether the input file is stacked  
        	- see `pBoxplot.args.stacked`
        - `params`: Other parameters for `geom_bar`, default: `Box()`  
        - `ggs`   : Extra ggplot2 statements  

!!! hint "pHeatmap"

    - **description**  
        Plot heatmaps.

    - **input**  
        - `infile:file`: The input matrix file  

    - **output**  
        - `outfile:file`: The heatmap  

    - **args**  
        - `ggs`: The ggplot items for heatmap  
        - `devpars`: The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
        - `dendro`: The parameters for control of the dendrogram. Default: `{'dendro': True}`  
        	- `dendro`: `True`: plot dendros for both rows and cols; `col`: only plot dendro for cols; `row`: only plot dendro for rows
        	- `rows`: The rownames to subset the rows and control the order of rows. Must a list. Only works when not plotting dendro for rows.
        	- `cols`: The colnames to subset the cols and control the order of cols. Must a list. Only works when not plotting dendro for cols.
        - `header`: The input file has header? Default: True  
        - `rownames`: The input file has rownames? Default: 1  
        - `rows`: Row selector  
        	- `all`: All rows
        	- `top:N`: Top N rows (original data ordered in descending order). N defaults to 100
        	- `bottom:N`: Bottom N rows. N defaults to 100
        	- `both:N`: Top N rows and bottom N rows. N defaults to 50
        	- `random:N`: Random N rows. N defaults to 50
        	- `random-both:N`: Random N rows from top part and N rows from bottom part. N defaults to 50
        - `cols`: Col selector (see `rows`).  

!!! hint "pHeatmap2"

    - **description**  
        Plot heatmaps using R package ComplexHeatmap. Example:
        ```
        bioprocs plot.pHeatmap2 
        	-i.infile MMPT.txt 
        	-i.annofiles:l:o PatientAnno.txt 
        	-args.params.row_names_gp 'r:fontsize5' 
        	-args.params.column_names_gp 'r:fontsize5' 
        	-args.params.clustering_distance_rows pearson 
        	-args.params.clustering_distance_columns pearson 
        	-args.devpars.width 5000 
        	-args.devpars.height 5000 
        	-args.draw.merge_legends 
        	-args.params.heatmap_legend_param.title AUC 
        	-args.params.row_dend_reorder 
        	-args.params.column_dend_reorder 
        	-args.params.top_annotation 'r:HeatmapAnnotation(Mutation = as.matrix(annos[,(length(groups)+1):ncol(annos)]), Group = as.matrix(annos[,groups]), col = list(Mutation = c(`0`="grey", `1`="lightgreen", `2`="green", `3`="darkgreen")), annotation_name_gp = fontsize8, show_legend = c(Group=F))' 
        	-args.params.right_annotation 'r:rowAnnotation(AUC = anno_boxplot(as.matrix(data), outline = F))' 
        	-args.helper 'fontsize8 = gpar(fontsize = 12); fontsize5 = gpar(fontsize = 8); groups = c("Group1", "Group2", "Group3")' 
        	-args.seed 8525
        ```

    - **input**  
        - `infile:file`: The input data file for the main heatmap.  
        - `annofiles:files`: The annotation files.  
        	- For now, they should share the same `args.anopts`

    - **output**  
        - `outfile:flie`: The plot.  

    - **args**  
        - `devpars`: The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
        - `draw`   : The parameters for `ComplexHeatmap::draw`  
        - `params` : Other parameters for `ComplexHeatmap::Heatmap`  
        - `anopts` : The options to read the annotation files.  
        - `inopts` : The options to read the input files.  
        - `seed`   : The seed. Default: `None`  
        - `helper` : The raw R codes to help defining some R variables or functions.  

    - **requires**  
        `R-ComplexHeatmap` (tested on c269eb425bf1b2d1713b9d5e68bf9f08bd8c7acb)

!!! hint "pScatterCompare"

    - **description**  
        Plot scatter plot to compare values of first 2 columns of input data

    - **input**  
        - `infile:file`: The input file containing a matrix with at least 2 columns  
        	- Other columns are groups used to group the scatter points
        	- Data must be normalized to [0, 1]

    - **output**  
        - `outfile:file`: The output plot  

    - **args**  
        - `ggs`: Extra expressions for ggplot. Note if geom_point is included, original geom_point will be ignored.  
        - `devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
        - `rownames`: Whether the input file has row names. Default: True  
        - `regr`: Whether draw the regression line. Default: False  
        - `corr`: The method to calculate the correlation. Default: `pearson`  
        	- Could be: `pearson`, `spearman` or `kendall`
        	- If it's neither of the three, no correlations will show.

!!! hint "pROC"

    - **description**  
        Generate ROC curves and output AUC.

    - **input**  
        - `infile:file`: The input matrix file.  
        	- Col1: rownames if args.rnames is True else label (0, 1 class)
        	- Col2: prediction values from model1
        	- ...

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `inopts`: The options for input file. Default: `Box(rnames = True, cnames = True)`  
        - `params`: The parameters for `plot.roc` from `utils/plot.r`  
        - `ggs`   : Additaional ggplot terms. Default:   
        	```python
        	Box({
        	    'style_roc': {},
        	    # show legend at bottom right corner
        	    'theme#auc': {'legend.position': [1, 0], 'legend.justification': [1, 0]} 
        	})
        	```
        - `devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  

!!! hint "pVenn"

    - **description**  
        Venn/UpsetR plots.

    - **input**  
        - `infile:file`: The input matrix  
        	- format:
        	```
        		category1	category2	category3
        	[e1]	0	1	1
        	[e2]	0	0	1
        	...
        	[eN]	1	0	0
        	```
        	rownames are not necessary but colnames are.
        - `metafile:file`: The metadata file for each category for upset plot.  
        	- format:
        	```
        		col1	col2	...	colN
        	category1	x1	y1	...	z1
        	category2	x2	y2	...	z2
        	...	...
        	categoryN	xN	yN	...	zN
        	```

    - **output**  
        - `outfile:file`: The plot  

    - **args**  
        - `tool`    : Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))  
        - `rnames`  : Whether input file has rownames. Default: False  
        - `params`  : Other params for `venn.diagram` or `upset`. Default: {}  
        - `devpars` : The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  

    - **requires**  
        [`r-VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram)
        [`r-UpSetR`](https://www.rdocumentation.org/packages/UpSetR)

!!! hint "pPie"

    - **description**  
        Plot piechart

    - **input**  
        - `infile:file`: The input file. Could be either:  
        	- Direct numbers of each category.
        	```
        	Group1	Group2
        	50	50
        	```
        	- Presence of each items in the category.
        	```
        		Group1	Group2
        	Item1	1	0
        	Item2	0	1
        	...
        	```

    - **output**  
        - `outfile:file`: the output plot  

    - **args**  
        - `rnames` : Whether the input file has row names. Default: `False`  
        - `ggs`    : Extra expressions for ggplot.  
        - `devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  

!!! hint "pManhattan"

    - **description**  
        Manhattan plot.

    - **input**  
        - `infile:file`: The input file. First 6 columns should be BED6, and then column:  
        	- 7: The raw pvalue.
        	- 8: The x-axis labels for the records.
        	- For example:
        		```
        		chr19	45604163	45604163	rs7255060	0	+	3.238E-03	+200K
        		chr19	45595277	45595277	rs10417101	0	+	3.870E-03	+200K
        		chr19	45394336	45394336	rs71352238	0	+	6.440E-03	-50K
        		chr19	45615857	45615857	rs6509194	0	+	1.298E-02	+250K
        		chr19	45594170	45594170	rs3178166	0	+	3.617E-02	+200K
        		chr19	45361574	45361574	rs3852856	0	+	2.070E-02	-100K
        		chr19	45220205	45220205	rs57090948	0	+	4.384E-02	-200K
        		chr19	45396219	45396219	rs157582	0	+	9.472E-03	-50K
        		chr19	45210634	45210634	rs10421830	0	+	1.375E-02	-250K
        		chr19	45228502	45228502	rs10422350	0	+	4.121E-02	-200K
        		```
        - `hifile:file`: The file with the record names (one per line) to highlight in the plot.  

    - **output**  
        - `outfile:file`: The plot. Default: `{{i.infile | fn}}.manht.png`  

    - **args**  
        - `inopts` : Options to read the input file. Default: `Box(cnames = False, rnames = False)`  
        - `ggs`    : Extra expressions for ggplot.  
        - `devpars`: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
## power

!!! hint "pSurvivalPower"

    - **description**  
        Do power analysis for survival analysis.
        See http://www.sample-size.net/sample-size-survival-analysis/

    - **input**  
        - `infile:file`: The input file, could be either:  
         	- detailed suvival data with [`patient`, ]`time`, `status`, `variable1`, `variable2`, ...; or
        	- ratios with `variable`, `survrate1`, `survrate2`, `ssratio`, where `survrate1` and
        		`survrate2` are survival rates in group1 and group2, respectively,
        		and `ssratio` is sample size ratio in group1/group2
        - `ngroup1`: The size of 1st group, for detailed input file  
        - `ngroup2`: The size of 2nd group, for detailed input file  
        - `ngroup3`: The size of 3rd group, for detailed input file  
        - `ngroup4`: The size of 4th group, for detailed input file  

    - **output**  
        - `outfile:file`: The output file with columns:  
        	- Variable: the variable (including paired groups)
        	- Alpha: the alpha value
        	- Beta: the beta value (1-power)
        	- SSize1: the sample size for group1
        	- SSize2: the sample size for group2
        	- Total: the total sample size

    - **args**  
        - `rnames`: Whether the detailed input file has row names. Default: `True`  
        - `plot`  : Plot the results? Default: `False`  
        - `ngroup`: Number of groups to divide into for detailed input file. Default: `2`  
        - `intype`: The input file type. Either `detailed` (default) or `ratio`  
        - `alphas`: The alpha values (two-sided). You need to *2 for one-sided. Default:  
        	- `[.005, .01, .05, .1]`
        - `betas` : 1-power. Default: `[.05, .1, .2]`  
## rank

!!! hint "pRank"

    - **description**  
        Convert values to ranks.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file with ranks.  

    - **args**  
        - `na`: Where to put the `NA` values.  
        	- `"first"` : Put `NA` first
        	- `"last"`  : Put `NA` last (default)
        	- `"remove"`: Remove `NA` values
        	- `"keep"`  : keep `NA` values
        - `tie`: How to deal with ties  
        	- `"average"` : Use average ranks (default)
        	- `"first"`   : Use the ranks come first
        	- `"last"`    : Use the ranks come last
        	- `"random"`  : Use the random ranks
        	- `"max"`     : Use the max ranks
        	- `"min"`     : Use the min ranks
        - `byrow`: Calculate ranks by row (instead of by column)? Default: `True`  
        - `reverse`: Take the reverse rank? Default: `True`  
        	- Large number gets higher rank (smaller rank index)
        	- `args.na` remains the same.
        - `inopts`: The input options:  
        	- `cnames`: Whether the input file has header. Default: `True`
        	- `rnames`: Whether the input file has row names. Default: `True`
        	- `delimit`: The separator of columns. Default: `\t`

!!! hint "pRankProduct"

    - **description**  
        Calculate the rank product of a set of ranks. Refer to [here](https://en.wikipedia.org/wiki/Rank_product)

    - **input**  
        - `infile:file`: The input file  
        	- Format:
        	```
        				Case1	Case2	...
        	Feature1	8.2  	10.1 	...
        	Feature2	2.3  	8.0  	...
        	...
        	```
        	- Or instead of values, you can also have ranks in the input file:
        	```
        				Rank1	Rank2	...
        	Feature1	2    	1    	...
        	Feature2	3    	2    	...
        	...
        	```

    - **output**  
        - `outfile:file`: The output file with original ranks, rank products and p-value if required  

    - **args**  
        - `informat`: The input format of the values. Whether they are real values (value) or ranks (rank). Default: value  
        	- Records will be ordered descendingly by value (Larger value has higher rank (lower rank index)).
        - `pval`: Whether to calculate the p-value or not. Default: True  
        - `plot`: Number of rows to plot. Default: 0 (Don't plot)  
        - `cex`: Font size for plotting. Default: 0.9  
        - `cnheight`: Colname height. Default: 80  
        - `rnwidth`: Rowname width. Default: 50  
        - `devpars`: device parameters for the plot. Default: `Box(res=300, width=2000, height=2000)`  
        - `inopts`: Options for reading the input file. Default: `Box(cnames=True, rnames=True, delimit="\t")`  
## resource

!!! hint "pTxt"

    - **description**  
        Download CSV format files.

    - **input**  
        - `in`: The name of the resource  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `cols`: Select the columns to keep. Default: '' (all cols)  
        - `rowfilter`: Filter rows. For example, to filter out rows not start with 'Chr':  
        	- `"lambda x: not x[0].startswith('Chr')"`
        	- Note that rowfilter applied before cols filter.
        - `urls`: Available resources and their urls.  
        - `gz`: Whether to gzip the output file.  

    - **requires**  
        [`curl`](https://en.wikipedia.org/wiki/CURL)
## rnaseq

!!! hint "pExprDir2Matrix"

    - **description**  
        Convert expression files to expression matrix
        File names will be used as sample names (colnames)
        Each gene and its expression per line.
        Suppose each expression file has the same rownames and in the same order.

    - **input**  
        - `indir:file`: the directory containing the expression files, could be gzipped  

    - **output**  
        - `outfile:file`: the expression matrix file  
        - `outdir:dir`: the directory containing expr file and plots  

    - **args**  
        - `pattern` : The pattern to filter files. Default `'*'`  
        - `namefunc`: Transform filename (no extension) as column name. Default: "function(fn) fn"  
        - `header`  : Whether each expression file contains header. Default: `False`  
        - `exrows`  : Rows to be excluded, regular expression applied. Default: `["^Sample", "^Composite", "^__"]`  
        - `boxplot` : Whether to plot a boxplot. Default: False  
        - `heatmap` : Whether to plot a heatmap. Default: False  
        - `histplot`: Whether to plot a histgram. Default: False  
        - `devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`  
        - `boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`  
        	- See ggplot2 documentation.
        - `heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`  
        - `histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`  

!!! hint "pExprFiles2Mat"

    - **description**  
        Merge expression to a matrix from single samples.

    - **input**  
        - `infiles:files`: The expression file from single samples, typically with 2 columns: gene and expression value  

    - **output**  
        - `outfile:file`: the expression matrix file  

    - **args**  
        - `fn2sample`: Transform filename (no extension) as column name. Default: "function(fn) fn"  
        - `inopts`   : Options to read input files. Default: `Box(rname = True, cnames = True)`  

!!! hint "pExprStats"

    - **description**  
        Plot the expression out.

    - **input**  
        - `infile:file`: The expression matrix (rows are genes and columns are samples).  
        - `gfile:file` : The sample information file. Determines whether to do subgroup stats or not.  
        	- If not provided, not do for all samples

    - **output**  
        - `outdir:dir`: The directory containing the plots  
        	- If `args.filter` is given, a filtered expression matrix will be generated in `outdir`.

    - **args**  
        - `inopts`: Options to read `infile`. Default: `Box(cnames = True, rnames = True)`  
        - `tsform`: An R function in string to transform the expression matrix (i.e take log).  
        - `filter`: An R function in string to filter the expression data.  
        - `plot`  : Which plot to do? Default: `Box(boxplot = True, histogram = True, qqplot = True)`  
        - `ggs`   : The ggs for each plot. Default:  
        	- `boxplot   = Box(ylab = {0: "Expression"})`,
        	- `histogram = Box(labs = {'x': 'Expression', 'y': '# Genes'})`,
        	- `qqplot    = Box()`
        - `params` : The params for each ggplot function.  
        - `devpars`: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`  

!!! hint "pBatchEffect"

    - **description**  
        Remove batch effect with sva-combat.

    - **input**  
        - `expr:file`: The expression file, generated by pExprDir2Matrix  
        - `batch:file`: The batch file defines samples and batches.  

    - **output**  
        - `outfile:file`: the expression matrix file  
        - `outdir:dir`: the directory containing expr file and plots  

    - **args**  
        - `tool`    : The tool used to remove batch effect. Default `'combat'`  
        - `hmrows`  : How many rows to be used to plot heatmap  
        - `plot`: Whether to plot  
        	- `boxplot`   : Whether to plot a boxplot. Default: False
        	- `heatmap`   : Whether to plot a heatmap. Default: False
        	- `histogram` : Whether to plot a histgram. Default: False
        - `devpars`    : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`  
        - `ggs`: The ggplot parameters  
        	- `boxplot`  : The ggplot parameters for boxplot. Default: `Box(ylab = {0: "Log2 Intensity"})`
        	- `heatmap`  : The ggplot parameters for heatmap. Default: `Box(theme = {'axis.text.y': 'r:element_blank()'})`
        	- `histogram`: The ggplot parameters for histgram. Default: `Box(labs = {'x': "Log2 Intensity", "y": "Density"})`

!!! hint "pUnitConversion"

    - **description**  
        Convert expression between units  
        See here: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ and  
        https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#fpkm  
        Available converstions:
        - `count -> cpm, fpkm/rpkm, fpkm-uq/rpkm-rq, tpm, tmm`
        - `fpkm/rpkm -> count, tpm, cpm`
        - `tpm -> count, fpkm/rpkm, cpm`
        - `cpm -> count, fpkm/rpkm, tpm`
        NOTE that during some conversions, `sum(counts/effLen)` is approximated to `sum(counts)/sum(effLen) * length(effLen))`

    - **input**  
        - `infile:file`: the expression matrix  
        	- rows are genes, columns are samples

    - **output**  
        - `outfile:file`: the converted expression matrix  

    - **args**  
        - `inunit`: The unit of input expression values. Default: `count`  
        - `outunit`: The unit of output expression values. Default: `tpm`  
        - `meanfl`: A file containing the mean fragment length for each sample by rows, without header.  
        	- Or a fixed universal estimated number. Default: `1` (used by TCGA)
        - `nreads`: The estimatied total number of reads for each sample. Default: `30000000`  
        	- Or you can pass a file with the number for each sample by rows, without header.
        	- In converting `fpkm/rpkm -> count`: it should be total reads of that sample
        	- In converting `cpm -> count`: it should be total reads of that sample
        	- In converting `tpm -> count`: it should be total reads of that sample
        	- In converting `tpm -> cpm`: it should be total reads of that sample
        	- In converting `tpm -> fpkm/rpkm`: it should be `sum(fpkm)` of that sample
        - `inform`: Transform the input to the unit specified. (sometimes the values are log transformed)  
        	- For example, if the `inunit` is `tpm`, but it's actually `log2(expr+1)` transformed,
        	- Then this should be `function(expr) 2^(expr - 1)`.
        	- Default: `None` (No transformation has been done)
        - `outform`: Transform to be done on the output expression values. Default: `None`  
        - `refexon`: the exome gff file, for RPKM/FPKM  
        	- `gene_id` is required for gene names

    - **requires**  
        [edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
        [coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen

!!! hint "pRNASeqDEG"

    - **description**  
        Detect DEGs for RNA-seq data

    - **input**  
        - `efile:file`: The expression matrix  
        - `gfile:file`: The group information  
        	- Like:
        	```
        	Sample	[Patient	]Group
        	sample1	[patient1	]group1
        	sample2	[patient1	]group1
        	sample3	[patient2	]group1
        	sample4	[patient2	]group2
        	sample5	[patient3	]group2
        	sample6	[patient3	]group2
        	```

    - **output**  
        - `outfile:file`: The DEG list  
        - `outdir:file`: The output directory containing deg list and plots  

    - **args**  
        - `tool`  : The tool used to detect DEGs. Default: 'deseq2' (edger is also available).  
        - `inopts`: Options to read `infile`. Default: `Box(cnames = True, rnames = True)`  
        - `cutoff`: The cutoff used to filter the results. Default: `0.05`  
        	- `0.05` implies `{"by": "p", "value": "0.05", "sign": "<"}`
        - `plot`  : The plots to do. Default:   
        	- `mdsplot`: True, MDS plot
        	- `volplot`: True, volcano plot
        	- `maplot `: True, MA plot
        	- `heatmap`: True, heatmap for DEGs
        	- `qqplot `: True, The QQplot for pvalues
        - `ggs`   : The ggs for each plot. Default:  
        	- `heatmap`: `Box(theme = {'axis.text.y': 'r:element_blank()'})`
        	- Not available for `mdsplot`.
        	- Others are empty `Box()`s
        - `params`: Parameters for each plot. Default:   
        	- `volplot`: `Box(pcut = 0.05, fccut = 2)`
        	- `maplot` : `Box(pcut = 0.05)`
        - `devpars`: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`  
        - `mapfile`: Probe to gene mapping file. If not provided, assume genes are used as rownames.  

!!! hint "pCoexp"

    - **description**  
        Get co-expression of gene pairs in the expression matrix.
## sambam

!!! hint "pSam2Bam"

    - **description**  
        Deal with mapped sam/bam files, including sort, markdup, and/or index

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output bam file  

    - **args**  
        - `tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools)  
        - `sambamba`         : The path of the sambamba. Default: sambamba  
        - `picard`           : The path of the picard. Default: picard  
        - `biobambam_bamsort`: The path of the biobambam's bamsort. Default: bamsort  
        - `samtools`: The path of the samtools. Default: samtools  
        - `sort`    : Do sorting? Default: True  
        	- If input is sam, tool is biobambam, this should be True
        - `index`  : Do indexing? Default: True  
        - `markdup`: Do duplicates marking? Default: False  
        	- `rmdup` for samtools will be called
        - `rmdup`  : Do duplicates removing? Default: False  
        - `tmpdir` : The tmp dir used to store tmp files. Default: <system default tmpdir>  
        - `sortby` : Sort by coordinate or queryname. Default: coordinate  
        - `nthread`: Default: 1  
        - `infmt`  : The format of input file. Default: <detect from extension> (sam|bam)  
        - `params` : Other parameters for `tool`. Defaut: ""  
        - `mem`    : The max memory to use. Default: "16G"  
        	- Unit could be G/g/M/m
        	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it

    - **requires**  
        [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)

!!! hint "pBamMarkdup"

    - **description**  
        Mark/remove duplicates for bam files

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output bam file  

    - **args**  
        - `tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools|bamutil)  
        - `sambamba`         : The path of sambamba. Default: sambamba  
        - `picard`           : The path of picard. Default: picard  
        - `biobambam_bamsort`: The path of biobambam's bamsort. Default: bamsort  
        - `samtools`         : The path of samtools. Default: samtools  
        - `bamutil`          : The path of bamutil. Default: bam  
        - `rmdup`            : Do duplicates removing? Default: False  
        - Samtools will anyway remove the duplicates
        - `tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>  
        - `nthread`          : Default: 1  
        - Not available for samtools and picard
        - `params`           : Other parameters for `tool`. Defaut: ""  
        - `mem`              : The max memory to use. Default: "16G"  
        - Unit could be G/g/M/m
        - Will be converted to -Xmx4G, and -Xms will be 1/8 of it

    - **requires**  
        [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)
        [bamutil](http://genome.sph.umich.edu/wiki/BamUtil#Programs)

!!! hint "pBamRecal"

    - **description**  
        Recalibrate a bam file

    - **input**  
        - `infile:file`: The bam file  

    - **output**  
        - `outfile:file`: The output bam file  

    - **args**  
        - `tool`    : The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)  
        - `gatk`    : The path of gatk, including java path. Default: `gatk`  
        - `samtools`: The path of samtools. Default: `samtools`  
        - `bamutil` : The path of bamutil. Default: `bam`  
        - `picard`  : The path of picard. Default: `picard`  
        - `params`  : Other parameters for `bam recab`. Default         : ""  
        	`RealignerTargetCreator` : Other parameters for `gatk RealignerTargetCreator`. Defaut: ""
        	`IndelRealigner`         : Other parameters for `gatk IndelRealigner`. Defaut: ""
        	`BaseRecalibrator`       : Other parameters for `gatk BaseRecalibrator`. Defaut: ""
        	`PrintReads`             : Other parameters for `gatk PrintReads`. Defaut: ""
        - `mem`: The max memory to use. Default: "32G"  
        - `knownSites`: The known polymorphic sites to mask out. Default: "" (Required for GATK)  
        - `ref`: The reference file. Required.  
        	- Will be converted to -Xmx4G, and -Xms will be 1/8 of it

    - **requires**  
        [gatk](https://software.broadinstitute.org/gatk)
        [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`

!!! hint "pBamReadGroup"

    - **description**  
        Add or replace read groups of a bam file

    - **input**  
        - `infile:file`: The bam file  

    - **output**  
        - `outfile:file`: The output bam file  

    - **args**  
        - `tool`                         : The tool used. Default: `picard` (picard|bamutil)  
        - `picard`                       : The path of picard. Default: `picard`  
        - `bamutil`                      : The path of bamutil. Default: `bam`  
        - `rg`                           : The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}  
        - `id` will be parsed from filename with "_LX_" in it if not given
        - `sm` will be parsed from filename
        - `params`                       : Other parameters for `tool`. Defaut: ""  
        - `mem`                          : The max memory to use. Default: "4G"  
        - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
        - `tmpdir`                       : The temporary directory. Default: <system tmpdir>  

    - **requires**  
        [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
        [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`

!!! hint "pBamReorder"

    - **description**  
        Reorder a sam/bam file by a given reference file using `picard ReorderSam`

    - **input**  
        - `infile:file`: The sam/bam file  

    - **output**  
        - `outfile:file`: The output bam file  

    - **args**  
        - `picard`                       : The path of picard. Default: `picard`  
        - `ref`                          : The reference file. Required  
        - `params`                       : Other parameters for `picard ReorderSam`. Defaut: ""  
        - `mem`                          : The max memory to use. Default: "4G"  
        - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
        - `tmpdir`                       : The temporary directory. Default: <system tmpdir>  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)

!!! hint "pBamMerge"

    - **description**  
        Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.

    - **input**  
        - `infiles:files`: Input sam/bam files to be merged  

    - **output**  
        - `outfile:file`: The merged bam file  

    - **args**  
        - `tool`     : The tool used to merge. Default: bamutil (picard|samtools|sambamba)  
        - `picard`   : The path of picard. Default: `picard`  
        - `bamutil`  : The path of bamutil. Default: `bam`  
        - `samtools` : The path of samtools. Default: `samtools`  
        - `sambamba` : The path of sambamba. Default: `sambamba`  
        - `params`   : Other parameters for `tool`. Defaut: ""  
        - `mem`      : The max memory to use. Default: "4G"  
        - Will be converted to -Xmx4G, and -Xms will be 1/8 of it, just for picard
        - `tmpdir`   : The temporary directory. Default: <system tmpdir>  
        - `nthread`  : # threads to use. Default: 1  
        - For picard, if nthread>1, USE_THREADING=true, otherwise USE_THREADING=false

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)

!!! hint "pBam2Gmut"

    - **description**  
        Call germline (snps and indels) from a call-ready bam file.

    - **input**  
        - `infile:file`: The input bam file  

    - **output**  
        - `outfile:file`: The vcf file containing the mutations  

    - **args**  
        - `tool`      : The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)  
        - `gatk`      : The path of gatk. Default: gatk  
        - `vardict`   : The path of vardict. Default: vardict  
        - `snvsniffer`: The path of snvsniffer. Default: SNVSniffer  
        - `samtools`  : The path of samtools. Default: samtools (used to generate reference index)  
        - `platypus`  : The path of platypus. Default: platypus  
        - `strelka`   : The path of strelka. Default: configureStrelkaGermlineWorkflow.py  
        - `cfgParams` : The params for `strelka` configuration. Default: ""  
        - `picard`    : The path of picard. Default: picard  
        - `mem`       : The memory to be used. Default: 32G  
        	- will be converted to -Xms4G -Xmx32G for java programs
        - `ref`: The reference file. Required.  
        - `gz`: Gzip output file? Default: False  
        - `tmpdir`: The temporary directory. Default: <system tmpdir>  
        - `params`: Other params for `tool`. Default: ""  

    - **requires**  
        [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
        [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
        [vardict](https://github.com/AstraZeneca-NGS/VarDict)
        [snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
        [platypus](http://www.well.ox.ac.uk/platypus)
        [strelka@2.7.1+](https://github.com/Illumina/strelka)

!!! hint "pBamPair2Smut"

    - **description**  
        Call somatic mutations from tumor-normal bam pair.

    - **input**  
        - `tumor:file`: The tumor bam file  
        - `normal:file`: The normal bam file  

    - **output**  
        - `outfile:file`: The vcf file  

    - **args**  
        - `tool`: The tool used to call mutations. Default: gatk (somaticsniper, strelka, snvsniffer, virmid, varidct)  
        - `gatk`: The path to gatk. Default: gatk  
        - `somaticsniper`: The path to gatk. Default: bam-somaticsniper  
        - `strelka`: The path to gatk. Default: configureStrelkaSomaticWorkflow.py  
        - `snvsniffer`: The path to gatk. Default: SNVSniffer  
        - `virmid`: The path to gatk. Default: virmid  
        - `vardict`: The path to gatk. Default: vardict  
        - `samtools`: The path to gatk. Default: samtools  
        - `picard`: The path to gatk. Default: picard  
        - `configParams`: The configuration parameters for `configureStrelkaSomaticWorkflow.py`. Default: `{}`  
        - `params`: The parameters for main programs. Default: `{}`  
        - `meme`: The memory. Default: 24G  
        - `ref`: The reference genom. Default: `params.ref.value`  
        - `gz`: Whether gzip the output vcf file. Default: False  
        - `nthread`: The number of threads to use. Default: 1  
        - `tmpdir`: The temporary directory. Default: `params.tmpdir.value`  

    - **requires**  
        [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
        [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
        [vardict](https://github.com/AstraZeneca-NGS/VarDict)
        [snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
        [platypus](http://www.well.ox.ac.uk/platypus)
        [strelka@2.7.1+](https://github.com/Illumina/strelka)

!!! hint "pBam2Cnv"

    - **description**  
        see `bioaggr.wxs.aBam2SCnv` and `bioaggr.wxs.aBam2GCnv`

!!! hint "pBamStats"

    - **description**  
        Get read depth from bam files.

    - **input**  
        - `infile:file`: The input bam file  

    - **output**  
        - `outfile:file`: The output statistic file  
        - `outdir:dir`: The directory containing result files and figures.  

    - **args**  
        - `tool`: The tool used to do the job. Default: bamstats  
        - `bamstats`: The path to bamstats. Default: bamstats  
        - `params`: Other params to main program. Default: `{}`  
        - `mem`: The memory to be used. Default: 16G  
        - `plot`: Whether plot the result. Default: True  

!!! hint "pBam2Fastq"

    - **description**  
        Convert sam/bam files to pair-end fastq files.

    - **input**  
        - `infile:file`: The sam/bam file.  
        	- Sam files only available for biobambam, picard

    - **output**  
        - `fqfile1:file`: The 1st match of paired reads  
        - `fqfile2:file`: The 2nd match of paired reads  

    - **args**  
        - `tool`     : The tool to use. Default: biobambam (bedtools, samtools, picard)  
        - `biobambam`: The path of bamtofastq of biobambam. Default: bamtofastq  
        - `bedtools` : The path of bedtools. Default: bedtools  
        - `samtools` : The path of samtools. Default: samtools  
        - `picard`   : The path of picard. Default: picard  
        - `mem`      : The memory to be used by picard. Default: 8G  
        - `gz`       : Whether gzip the output files. Default: True  
        - `params`: : Other params for `tool`. Default: ''  
        - `tmpdir`   : The tmpdir. Default: `__import__('tempfile').gettempdir()`  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)
        [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

!!! hint "pBam2FastqSE"

    - **description**  
        Convert sam/bam files to single-end fastq files.

    - **input**  
        - `infile:file`: The sam/bam file.  
        	- Sam files only available for biobambam, picard

    - **output**  
        - `fqfile:file`: The fastq file  

    - **args**  
        - `tool`     : The tool to use. Default: biobambam (bedtools, samtools, picard)  
        - `biobambam`: The path of bamtofastq of biobambam. Default: bamtofastq  
        - `bedtools` : The path of bedtools. Default: bedtools  
        - `samtools` : The path of samtools. Default: samtools  
        - `picard`   : The path of picard. Default: picard  
        - `mem`      : The memory to be used by picard. Default: 8G  
        - `gz`       : Whether gzip the output files. Default: True  
        - `params`: : Other params for `tool`. Default: ''  
        - `tmpdir`   : The tmpdir. Default: `__import__('tempfile').gettempdir()`  

    - **requires**  
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)
        [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

!!! hint "pBam2Counts"

    - **description**  
        Extract read counts from RNA-seq bam files.

    - **input**  
        - `infile:file`: The input bam files  

    - **outfile**  
        - `outfile:file`: The count file  

    - **args**  
        - `tool`: The tool used to extract counts. Default: ht-seq  
        - `htseq`: The path of htseq-count.  
        - `params`: Other params for main program.  
        - `refgene`: The reference gene in GTF format.  

    - **requires**  
        [`htseq`](https://htseq.readthedocs.io/)

!!! hint "pBamIndex"

    - **description**  
        Index bam files.

    - **input**  
        - `infile:file`: The input bam file  

    - **output**  
        - `outfile:file`: The symbolic link to the input file  
        - `outidx:file` : The index file  

    - **args**  
        - `samtools`: Path to samtools. Default: `params.samtools`  
        - `params`  : Other parameters for samtools. Default: `Box(b = True)`  
        - `nthread` : # threads to use. Default: `1`  
## seq

!!! hint "pConsvPerm"

    - **description**  
        Generate a null distribution of conservation scores.

    - **input**  
        - `seed`: The seed to generate the random regions. Default: None  

    - **output**  
        - `outfile:file`: A file with mean conservation scores sorted descendingly.  

    - **args**  
        - `len`: The length of a random region. Default: 50  
        - `nperm`: Number of permutations. Default: 1000  
        - `gsize`: The chrom size file.  
        - `bedtools`: The path of bedtools.  
        - `bwtool`: The path of bwtool.  
        - `consvdir`: The directory containing bigwig files of conservation scores  
        	- The bigwig file should start with chr name: chrN.*

    - **requires**  
        [bwtool](https://github.com/CRG-Barcelona/bwtool)
        [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

!!! hint "pConsv"

    - **description**  
        Get the conservation scores of regions.
        It uses wigFix to find the conservation scores.
        But first you have to convert those wigFix.gz files to bigWig files using ucsc-wigToBigWig

    - **input**  
        - `bedfile:file`: The bedfile with regions in the same chromosome  
        - `permfile:file`: The permutaiton file generated by `pConsvPerm`, used to calculate p-values  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `consvdir`: The bigwig directory, the bigwig files must be named as "chrN.*.bw"  
        	- For example: `chr1.phyloP30way.bw`
        - `bwtool`: The path of bwtool executable. Default: `bwtool`  
        - `bedtools`: The path of bedtools executable. Default: `bedtools`  
        - `pval`: Whether calculate pvalue of each region. Default: False  
        	- In this case, the `i.permfile` can be ignored.

    - **requires**  
        [bwtool](https://github.com/CRG-Barcelona/bwtool)
        [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

!!! hint "pPromoters"

    - **description**  
        Get the promoter regions in bed format of a gene list give in infile.
        Gene names are supposed to be normalized by `gene.pGeneNameNorm`

    - **input**  
        - `infile:file`: the gene list file  

    - **output**  
        - `outfile:file`: the bed file containing the promoter region  

    - **args**  
        - `region`: The region to output. Default: `Box(up = 2000, down = None, withbody = False)`  
        	- `up`: The upstream distance to TSS.
        	- `down`: The downstream distance to TSS. Defaults to `args.region.up` if `None`
        	- `withbody`: Include gene body in the output region? Default: `False`
        - `notfound`: How to deal with records can't be found. Default: `skip`  
        	- `skip` : Skip the record
        	- `error`: Report error and exit
        - `refgene`: The ref gene file. Default: `@params.refgene`  
        - `inopts` : The options for input file. Default: `Box(cnames = False, genecol = 0, delimit = "\t")`  
        	- `cnames`: Whether the input file has header
        	- `genecol`: The 0-based index or the colname of gene column.
        	- `delimit`: The delimit of the input file.
        - `outopts`: The options for output file. Default: `Box(cnames = False)`  
## snp

!!! hint "pRs2Bed"

    - **description**  
        Find coordinates for SNPs in BED format.

    - **input**  
        - `snpfile:file`: the snp file, each snp per line  

    - **output**  
        - `outfile:file`: the result file, could be a 3 or 6-col bed file.  

    - **args**  
        - `dbsnp`: The dbsnp vcf file  
        - `notfound`: What to do if the snp is not found. Default: skip  
        - `inopts`: The input options for input file  
        ``

!!! hint "pSnp2Avinput"

    - **description**  
        Convert SNP list to avinput to ANNOVAR.

    - **input**  
        - `snpfile:file`: the snp file, each snp per line  

    - **output**  
        - `outfile:file`: the result avinput file  

    - **args**  
        - `genome`: default: hg19  
        - `snpver`: default: snp147  

    - **requires**  
        [`python-cruzdb`](https://github.com/brentp/cruzdb)
## snparray

!!! hint "pGistic"

    - **description**  
        Runing GISTIC to get CNV results.
        see: ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm

    - **input**  
        - `segfile:file`: Segmentation File  
        - `mkfile:file` : Markers File  
        - `alfile:file` : Array List File  
        - `cnvfile:file`: CNV File  

    - **output**  
        - `outdir:dir`: The output directory  
        	- All Lesions File (all_lesions.conf_XX.txt, where XX is the confidence level)
        	- Amplification Genes File (amp_genes.conf_XX.txt, where XX is the confidence level)
        	- Deletion Genes File (del_genes.conf_XX.txt, where XX is the confidence level)
        	- Gistic Scores File (scores.gistic)
        	- Segmented Copy Number (raw_copy_number.pdf)

    - **args**  
        - `gistic`: The path to gistic.  
        - `genome`: The genome used to select refgene file from refgenefiles.  
        - `mcr`: The mcr path  
        - `params`: Other params for gistic  

!!! hint "pSNP6Genotype"

    - **description**  
        Call genotypes from GenomeWideSNP_6 CEL file

    - **input**  
        - `celfile:file`: the CEL file  

    - **output**  
        - `outfile:file`: the outfile containing probe name and genotypes  
        - format: `<Probe name>\t<genotype>`
        - `<genotype>` = 0: AA, 1: AB, 2: BB

    - **requires**  
        [bioconductor-crlmm](http://bioconductor.org/packages/release/bioc/html/crlmm.html)

!!! hint "pGenoToAvInput"

    - **description**  
        Convert the genotype called by pSNP6Genotype to [ANNOVAR input file](http://annovar.openbioinformatics.org/en/latest/user-guide/input/#annovar-input-file) using dbSNP identifiers.	

    - **input**  
        - `genofile:file`: the genofile generated by pSNP6Genotype, must be sorted by probe names  
        - `annofile:flie`: the annotation file downloaded from http://www.affymetrix.com/support/technical/annotationfilesmai.affx  
        	- Could be in .gz format

    - **output**  
        - `outfile:file`: the avinput file  

    - **requires**  
        [python-read2](https://github.com/pwwang/read2)
## sql

!!! hint "pCreateTable"

    - **description**  
        Create tables in the database

    - **input**  
        - `dsn`: The dsn to connect to the database  
        	- currently support `sqlite:file=...`
        - `schema:file`: The schema file  
        	- could be a pure schema file:
        	```
        	Field	Type	Statement
        	ID	INT	PRIMARY KEY
        	...
        	```
        	- or a data file with header

    - **output**  
        - `dsn`: The dsn  

    - **args**  
        - `intype`: The input file schema file or a data file. Default: `schema`  
        - `drop`: Force creating the table (drop the pre-existing table)  
        - `delimit`: The delimit of input file. Default: `\\t`  

!!! hint "pImportData"

    - **description**  
        Create tables and import the data

    - **input**  
        - `dsn`: The dsn to connect to the database  
        	- currently support `sqlite:file=...`
        - `datafile:file`: The schema file  
        	- must have header

    - **output**  
        - `dsn`: The dsn  

    - **args**  
        - `delimit`: The delimit of input file. Default: `\\t`  

!!! hint "pUpdateTable"

    - **description**  
        Update table using sql.

    - **input**  
        - `dsn`: The dsn to connect to the database  
        	- currently support `sqlite:file=...`

    - **output**  
        - `dsn`: The dsn  

    - **args**  
        - `sql`: The sql to update the table (list)  

!!! hint "pSelectTable"

    - **description**  
        Select data from table and dump it.

    - **input**  
        - `dsn`: The dsn to connect to the database  
        	- currently support `sqlite:file=...`

    - **output**  
        - `outfile:file`: The dumped file  

    - **args**  
        - `sql`: The sql to select data from the table (list)  
## stats

!!! hint "pMetaPval"

    - **description**  
        Calculate a meta-pvalue using different methods

    - **input**  
        - `infile:file`: The infile containing multiple pvalues for each entries.   
        	- Could be two types (see `args.intype`)
        	- `matrix`: A matrix with rows as entries (rownames are optional), columns as cases
        		```
        		        Case1   Case2   ...  CaseN
        		Entry1  1.33e-2 NA      ...  1.77e-10
        		Entry2  2.66e-2 4.22e-5 ...  1.71e-3
        		... ...
        		EntryM  NA      0.00013 ...  4.11e-3
        		```
        	- `melt`: Rows are entries from each case 
        	  `args.rnames` should be `False`, but entry names are required in the first column,
        	  column names are optional (have to be 3 columns)
        		```
        		Entry   Pvalue   Case
        		Entry1  1.33e-2  Case1
        		Entry1  1.77e-10 CaseN
        		Entry2  2.66e-2  Case1
        		Entry2  4.22e-5  Case2
        		Entry2  1.71e-3  CaseN
        		... ...
        		EntryM  0.00013  Case2
        		EntryM  4.11e-3  CaseN
        		```

    - **output**  
        - `outfile:file`: The output file containing the meta-pvalues. Default: `{{i.infile | fn}}.meta{{i.infile | ext}}`  

    - **args**  
        - `intype`: The type of the input file. Default: `matrix` (see `i.infile`)  
        - `inopts`: The input options to read the input file. Default: `Box(rnames = True, cnames = True)`  
        - `method`: The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method, aka `fisher`)  
        	- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
        	- See: https://www.rdocumentation.org/packages/metap/versions/0.8
        - `na`    : How to deal with `NA` p-values. Default: `skip` (just don't count it)  
        	- Or a numeric value to replace it with (e.g.: `1`).

    - **requires**  
        [`r-matep`](https://www.rdocumentation.org/packages/metap/)

!!! hint "pSurvival"

    - **description**  
        Survival analysis

    - **input**  
        - `infile:file`: The input file (header is required).  
        	- col1: rownames if args.inopts.rnames = True
        	- col2: the survival time
        	- col3: the status. 0/1 for alive/dead or 1/2 for alive dead
        	- col4: var1.
        	- ... other variables

    - **output**  
        - `outfile:file`: The outfile containing the pvalues  
        - `outdir:dir`  : The output directory containing the pval files and plots  

    - **args**  
        - `inunit`    : The time unit in input file. Default: days  
        - `outunit`   : The output unit for plots. Default: days  
        - `nthread`   : Number of threads used to perform analysis for groups. Default: 1  
        - `inopts`    : The options for input file  
        	- `rnames`: Whether input file has row names. Default: True
        - `combine`   : Whether combine groups in the same plot. Default: `Box()`  
        	- `nrow`: The number of rows. Default: 1
        	- `ncol`: The number of cols. Default: 1
        - `devpars`   : The device parameters for png. Default: `{res:300, height:2000, width:2000}`  
        	- The height and width are for each survival plot. If args.combine is True, the width and height will be multiplied by `max(arrange.ncol, arrange.nrow)`
        - `covfile`   : The covariant file. Require rownames in both this file and input file.  
        - `ngroups`   : Number of curves to plot (the continuous number will divided into `ngroups` groups.  
        - `params`    : The params for `ggsurvplot`. Default: `Box({'risk.table': True, 'conf.int': True, 'font.legend': 13, 'pval': '{method}\np = {pval}'})`  
        	- You may do `ylim.min` to set the min ylim. Or you can set it as 'auto'. Default: 0. 
        - `ggs`       : Extra ggplot2 elements for main plot. `ggs.table` is for the risk table.  
        - `pval`      : The method to calculate the pvalue shown on the plot. Default: True (logrank)  
        	- Could also be `waldtest`, `likeratio` (Likelihoold ratio test)
        - `method`    : The method to do survival analysis.   

    - **requires**  
        [`r-survival`](https://rdrr.io/cran/survival/)
        [`r-survminer`](https://rdrr.io/cran/survminer/)

!!! hint "pPostSurvival"

    - **description**  
        Statistic comparison between groups after survival analysis.

    - **input**  
        - `infile:file`: The result file from `pSurvival`  
        - `survfile:file`: The survival data. See format of infile of `pSurvival`  

    - **output**  
        - `outfile:file`: The output excel file.  

    - **args**  
        - `covfile`: The covariant file. Require rownames in both this file and input file.  
        - `methods`: A list of testing methods  
        	- `wilcox`: Wilcox rank sum test
        	- `t`: t-test
        	- `chisq`: chisquare-test
        - `inopts`: The input options for `i.survfile`.  
        	- `rnames`: whether the file has row names. This has to be True if `args.covfile` provided.

!!! hint "pBin"

    - **description**  
        Bin the data in columns.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file. Default: `{{i.infile | stem}}.binned{{i.infile | ext}}`  

    - **args**  
        - `inopts`: The input options.  
        	- `delimit`: The delimiter. Default: `\t`
        	- `rnames`: Whether input file has row names. Default: `False`
        	- `cnames`: Whether input file has column names. Default: `True`
        	- Other arguments available for `read.table`
        - `binopts`: The default bin options:  
        	- `nbin`: Number of bins.
        	- `step`: The step of binning.
        	- `nan`:  What to do if the value is not a number. Default: `skip`
        		- `skip/keep`: Keep it
        		- `as0`: Treat it as 0
        	- `out`: The out value. Default: `step`
        		- `step`: Use the step breaks
        		- `lower/min`: Use the min value of the records in the bin
        		- `upper/max`: Use the max value of the records in the bin
        		- `mean`: Use the mean value of the records in the bin
        		- `median`: Use the median value of the records in the bin
        		- `binno`: Use the bin number (empty bins will be skipped).
        - `cols`: The detailed bin options for each column.   
        	- If not provided (`None`), all columns will use `binopts`. 
        	- If column specified, only the specified column will be binned.
        	- Column indices can be used. It's 1-based.

!!! hint "pQuantileNorm"

    - **description**  
        Do quantile normalization

    - **input**  
        - `infile:file`: The input matrix  

    - **output**  
        - `outfile:file`: The output matrix. Default: `{{i.infile | bn}}`  

!!! hint "pChiSquare"

    - **description**  
        Do chi-square test.

    - **input**  
        - `infile:file`: The input file.  

    - **output**  
        - `outfile:file` : The output file containing Xsquare, df, pval and method  
        - `obsvfile:file`: The observation matrix  
        - `exptfile:file`: The expectation matrix  

    - **args**  
        - `intype`: The type of the input file:  
        	- `count` (default): The contingency table
        	```
        	#         | Disease | Healthy |
        	# --------+---------+---------+
        	#   mut   |   40    |   12    |
        	# non-mut |   23    |   98    |
        	# --------+---------+---------+
        	```
        	- `raw`: The raw values:
        	```
        	# Contingency table rows: Mut, Non
        	# Contingency table cols: Disease, Healthy
        	#
        	#         | S1 | S2 | ... | Sn |
        	# --------+----+----+-----+----+
        	# Disease | 1  | 0  | ... | 1  |
        	# Healthy | 0  | 1  | ... | 0  |
        	# --------+----+----+-----+----+
        	# Mut     | 1  | 0  | ... | 1  |
        	# Non     | 0  | 1  | ... | 0  |
        	```
        - `ctcols`: The colnames of contingency table if input file is raw values  
        	- You may also specify them in the head of the input file

!!! hint "pFisherExact"

    - **description**  
        Do fisher exact test.

    - **input**  
        - `infile:file`: The input file.  

    - **output**  
        - `outfile:file` : The output file containing confInt1, confInt2, oddsRatio, pval, alternative and method.  

    - **args**  
        - `intype`: The type of the input file:  
        	- `count` (default): The contingency table
        	```
        	#         | Disease | Healthy |
        	# --------+---------+---------+
        	#   mut   |   40    |   12    |
        	# non-mut |   23    |   98    |
        	# --------+---------+---------+
        	```
        	- `raw`: The raw values:
        	```
        	# Contingency table rows: Disease, Healthy
        	# Contingency table cols: Mut, Non
        	#
        	#    | Disease Healthy | Mut  Non  |
        	# ---+--------+--------+-----+-----+
        	# S1 |    1   |    0   |  0  |  1  |
        	# S2 |    0   |    1   |  1  |  0  |
        	# .. |   ...  |   ...  | ... | ... |
        	# Sn |    0   |    1   |  0  |  1  |
        	#
        	```
        - `ctcols`: The colnames of contingency table if input file is raw values  
        	- You may also specify them in the head of the input file

!!! hint "pPWFisherExact"

    - **description**  
        Do pair-wise fisher exact test.
        Commonly used for co-occurrence/mutual-exclusivity analysis.
        P-value indicates if the pairs are significantly co-occurred or mutually exclusive.
        Co-occurrence: Odds ratio > 1
        Mutual-exclusivity: Odds ratio < 1

    - **input**  
        - `infile:file`: The input file.  

    - **output**  
        - `outfile:file` : The output file containing confInt1, confInt2, oddsRatio, pval, qval, alternative and method.  

    - **args**  
        - `intype`: The type of the input file:  
        	- `pairs`: The contingency table
        	```
        	#
        	# A+	B+	4
        	# A-	B-	175
        	# A+	B-	12
        	# A-	B+	1
        	#
        	```
        	- `raw` (default): The raw values:
        	```
        	#
        	#    | A | B | ... | X |
        	# ---+---+---+-----+---+
        	# S1 | 1 | 0 | ... | 1 |
        	# S2 | 0 | 1 | ... | 0 |
        	# .. | 0 | 0 | ... | 1 |
        	# Sn | 0 | 1 | ... | 1 |
        	#
        	```
        - `padj`: The p-value adjustment method, see `p.adjust.methods` in R. Default: `BH`  
        - `rnames`: If the input file has rownames for `raw` input type.  

!!! hint "pMediation"

    - **description**  
        Do mediation analysis

    - **input**  
        - `infile:file`: The input file (a matrix or data.frame). Example:  
        	```
        	     V1   V2   V3
        	S1   1    2    3
        	S2   4    1    8
        	... ...
        	Sn   3    3    1
        	```
        - `casefile:file`: The mediation options. Example:  
        	```
        	Case1   lm      V3~V2|V1
        	Case2   lm,glm  V3~V1|V2
        	```
        	- No column names, but implies `Case`, 'Model' and `Formua`.
        	- `\t` as delimiter.
        	- This file is optional. If it is not provided, `args.case` is required.
        	- If this file is provided, `args.case` is ignored
        	- For different models, model for Mediator comes last. For Case2, the models will be:
        		- lm(V3 ~ V2 + V1) and glm(V2 ~ V1)

    - **output**  
        - `outfile:file`: The result file.  
        - `outdir:dir`  : The output directory containing output file and plots.  

    - **args**  
        - `inopts`: The options for input file. Default: `Box(cnames = True, rnames = True)`  
        	- `cnames`: Whether the input file has column names
        	- `rnames`: Whether the input file has row names
        - `medopts`: The options for mediation analysis.  
        	- `boot`: Use bootstrap?
        	- `sims`: How many time simulations?
        - `cov`: The covariate file. Default: ``  
        - `pval`: The pvalue cutoff. Default: `0.05`  
        - `fdr` : Method to calculate fdr. Use `False` to disable. Default: `True` (`BH`)  
        - `plot`: Parameters for `plot.mediate`? Use `False` to disable plotting. Default: `Box()`  
        	- Only case with pvalue < `args.pval` will be plotted.
        	- To plot all cases, use `args.pval = 1`
        - `nthread`: Number of threads to use for different cases. Default: `1`  
        - `devpars`: device parameters for the plot. Default: `Box(res = 300, width = 2000, height = 2000)`  
        - `case`   : Define cases, each case should have `model` and `fmula`.  
        	- If you only have one case, then it could be: `Box(model = 'lm', fmula = 'Y~X|M')`
        	  In this case, `{{i.infile | fn2}}` will be used as case name
        	- For multiple cases, this should be a dict of cases: 
        	  `Box(Case1 = Box(model='lm', fmula='Y~X|M'), Case2 = ...)`

!!! hint "pLiquidAssoc"

    - **description**  
        Do liquid association analysis

    - **input**  
        - `infile:file`: The input file with input data, where LA will be done on rows.  
        	```
        	     S1   S2 ... ... Sn
        	G1   1    2  ... ... 9
        	G2   3    1  ... ... 1
        	... ...
        	Gm   9    2  ... ... 3
        	```
        - `casefile:file`: Defining the groups (X, Y, Z) and the cases. If case (3rd col) is not provided, all will be treated as one case.  
        	- Group "Z" is required. You can also specify group "X", then the rest will be group "Y"
        	```
        	G1   X   Case1
        	G2   X   Case1
        	Gx   Z   Case1
        	Gy   Z   Case1
        	```

    - **output**  
        - `outfile:file`: The results of the analysis  
        - `outdir:dir`  : The output directory containing the result file and plots.  

    - **args**  
        - `inopts` : The options for reading input file  
        - `zcat`   : Whether the group "Z" is categorical. Default: `False`  
        	- If it is, then `stein lemma` is not suitable, we will calculate LA manually (`E(g'(z)`)
        - `pval`   : The pval cutoff. Default: `0.05`  
        - `fdr`    : The method to calculate FDR. Use `False` to disable. Default: `True` (BH)  
        - `nthread`: The number of threads to use. Default: `1` (WCGNA requires)  
        - `plot`   : Whether do plotting or not. Default: `False`  
        - `devpars`: device parameters for the plot. Default: `Box(res = 300, width = 2000, height = 2000)`  

    - **requires**  
        [r-fastLiquidAssociation](https://github.com/pwwang/fastLiquidAssociation)

!!! hint "pHypergeom"

    - **description**  
        Do hypergeometric test.

    - **input**  
        - `infile:file`: The input file, could be raw data (presence (1) and absence (0) of elements) or number of overlapped elements and elements in each category.  
        	- Set `args.intype` as `raw` if it is raw data. The population size `args.N` is required
        	- Set `args.intype` as `numbers` (or any string except `raw`) if it is numbers. You can specified explicit header: `k` = overlapped elements, `m` = size of set 1, `n` = size of set 2 and `N` = the population size. If `N` not included, then `args.N` is required

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `intype`: the type of input file. Default: `raw`. See `infile:file`  
        - `inopts`: The options for input file.  
        	- `cnames`: Whether the input file has column names
        	- `rnames`: Whether the input file has row names
        - `N`: The population size. Default: `None`  

!!! hint "pChow"

    - **description**  
        Do Chow-Test

    - **input**  
        - `infile:file`: The input file for data to do the regressions. Example:  
        	```
        		    X1  X2  X3  X4 ... Y
        		G1  1   2   1   4  ... 9
        		G2  2   3   1   1  ... 3
        		... ...
        		Gm  3   9   1   7  ... 8
        	```
        - `groupfile:file`: Specify the groups to compare. You may also specify the cases. The Chow-Test will be done between the group for each case. Example:  
        	```
        		G1	Group1	Case1
        		G2	Group1	Case1
        		... ...
        		Gs	Group2	Case1
        		Gt	Group2	Case1
        		Gt	Group1	Case2
        		... ...
        		Gu	Group1	Case2
        		... ...
        		Gz	Group2	Case2
        	```
        	- In such case, the test will be done between Group1 and Group2 for Case1 and Case2, respectively.
        	- Instances can be resued (Gt in the example)
        	- If cases not provided, all will be treated as one case.
        - `casefile:file`: Define the formula (which columns to use for each case). Example:  
        	```
        	Case1	Y ~ X1
        	Case2	Y ~ X2
        	```

    - **output**  
        - `outfile:file`: The result file of chow test. Default: `{{i.infile | fn}}.chow/{{i.infile | fn}}.chow.txt`  
        - `outdir:dir`: The output directory, containing the output file, results of regressions and plots.  

    - **args**  
        - `inopts`: The options for input file.  
        	- `cnames`: Whether the input file has column names. Default: `True`
        	- `rnames`: Whether the input file has row names. Default: `True`
        - `cov`: The covariate file. `inopts.rnames` required and this file should have row names too. Default: `''`  
        - `fdr`   : Calculate FDR or not. Use `False` to disable. If `True` will use `BH` method, otherwise, specify the method (see `R`'s `p.adjust`).  
        - `pval`: The pvalue cutoff. Default: `0.05`  
        - `plot`: Whether plot the regressions. Default: `False`  
        - `ggs` : The extra ggs for the plot.  
        - `devpars`: device parameters for the plot. Default: `Box(res = 300, width = 2000, height = 2000)`  

!!! hint "pCorr"

    - **description**  
        Calculate the correlation coefficient for the input matrix

    - **input**  
        - `infile:file`: The input file of data to calculate correlations.  

    - **output**  
        - `outfile:file`: The output file containing the correlation coefficients  
        - `outdir:dir`  : The output directory containing the outfile and the plot  

    - **args**  
        - `outfmt`: The output format. Could be `matrix` or `pairs` (default)  
        - `metohd`: The method used to calculate the correlation coefficient. Default: `pearson`. Could also be `spearman` or `kendall`  
        - `byrow`: Calculate the correlation coefficient by row or by col. Default: `True`  
        - `inopts`: The input options:  
        	- `cnames`: Whether the input file has header. Default: `True`
        	- `rnames`: Whether the input file has row names. Default: `True`
        	- `delimit`: The separator of columns. Default: `\t`
        - `plot`: Whether output a correlation plot. Default: `False`  
        - `params`: The params for `plot.heatmap` in `utils/plot.r`  
        - `ggs`: The extra ggplot2 statements.  
        - `devpars`: The parameters for the plot device. Default: `Box(height = 2000, width = 2000, res = 300)`  

    - **requires**  
        R packages: `ggplot2` and `reshape`

!!! hint "pCorr2"

    - **description**  
        Calculate correlation coefficient between instances of two files
        Don't do it between instances within the same file.

    - **input**  
        - `infile1:file`: The first file. See input of `pCorr`  
        - `infile2:file`: The second file.  
        	- must have same number of columns with `infile1`

    - **output**  
        - `outfile:file`: The output file.  
        - `outdir:dir`  : The output directory containing output file and other files:  
        	- pvalues/fdr file and plots

    - **args**  
        - `pval`  : Whether output pvalue. Default: `False`  
        - `fdr`   : Whether output qvalue. Default: `False`  
        - `outfmt`: The output format. `pairs` (default) or `matrix`  
        - `plot`  : Whether plot a heatmap or not. Default: `False`  
        - `params`: The params for `plot.heatmap` in `utils/plot.r`  
        - `ggs`: The extra ggplot2 statements.  
        - `devpars`: The parameters for the plot device. Default: `Box(height = 2000, width = 2000, res = 300)`  

!!! hint "pDiffCorr"

    - **description**  
        Test correlation differences using Fisher Z method.

    - **input**  
        - `infile:file`: The entire dataset used to calculate correlations. Rownames and colnames are required. Example:  
        	```
        		    S1  S2  S3  S4 ... Sn
        		G1  1   2   1   4  ... 9
        		G2  2   3   1   1  ... 3
        		... ...
        		Gm  3   9   1   7  ... 8
        	```
        - `samfile:file`: The sample groups, between which you want to compare the correlations. You can also specify one sample to multiple groups, and assign with different cases. Example:  
        	```
        		S1	Healthy
        		S2	Healthy
        		S3	Disease
        		... ...
        		Sn	Disease
        	```
        - `casefile:file`: Assign the cases to compare. If not provided, it will do for every possible combination. Example:  
        	```
        		Healthy	Disease
        	```
        - `groupfile:file`: Specify groups for rows, then the correlation will be only done within the pairs, each of which is from different groups (only 2 allowed). If not provided, it will investigate for every possible row pairs. Example:  
        	```
        		G1	Kinase
        		G2	Kinase
        		... ...
        		Gm	TF
        	```

    - **output**  
        - `outfile:file`: The pairs under different cases that their correlations have been changed significantly  
        - `outdir:dir`: The output directory containing plots and other output files.  

    - **args**  
        - `inopts`: The input options for `infile`. `cnames` and `rnames` have to be `True`  
        - `method`: The method to calculate the correlation. Default: `pearson`  
        - `pval`  : The pvalue cutoff to define correlation change significance. Default: `0.05`  
        - `fdr`   : Calculate FDR or not. Use `False` to disable. If `True` will use `BH` method, otherwise, specify the method (see `R`'s `p.adjust`).  
        - `fdrfor`: Do FDR calculation for each case (`case`) or all instances (`all`). Default: `case`  
        - `plot`  : Plot the correlation for cases? Default: `False`  
        - `ggs`   : `ggs` items for the plot.  
        - `devpars`: The device parameters for the plot.  

!!! hint "pBootstrap"

    - **description**  
        Do bootstrapping resampling

    - **input**  
        - `infile:file`: The input data file  

    - **output**  
        - `outfile:file`: The output file with the bootstrapped statistics values  
        	- depends on the `args.stats` function
        - `outdir:dir`: The directory to save the outfile and figures.  

    - **args**  
        - `inopts`: The options to read the input file. Default: `Box(cnames = True, rnames = True)`  
        - `params`: Other parameters for `boot` function from R `boot` package  
        - `nthread`: # of threads(cores) to use. Default: `1`  
        - `n`: Sampling how many times? Default: `1000`  
        - `stats`: The function to generate statistics for output. Default: `function(x) x`  
        	- Default to use all data
        	- This function can return a multiple statistics in a vector
        	- The argument `x` is the data generate for each sampling. 
        	- Unlink the `statistic` argument from `boot`, to make it convenient, we don't put the `index` here.
        - `plot`: Plot the statistics? Default: `all` (plot all statistics)  
        	- You may also specify indices. For example: `[1, 2]` to plot the 1st and 2nd statistics
        	- Use `False` or `None` to disable plotting
        - `devpars`: The device parameters for the plot.  

!!! hint "pPCA"

    - **description**  
        Perform PCA analysis. Example:
        ```
        bioprocs stats.pPCA 
        	-i.infile Cellline_t.txt 
        	-i.annofile CLAnno.txt 
        	-args.plots.clplot.repel 
        	-args.plots.clplot.shape 3 
        	-args.plots.clplot.ggs.geom_point 'r:list(aes(shape = Cellline), color = "#2b6edb", data = anno)' 
        	-args.seed 8525 
        	-args.plots.cluster.centers 2 
        	-args.plots.clplot.show-clust-cent 0 
        	-args.plots.cluster.npcs 2
        ```

    - **input**  
        - `infile:file`: The matrix to do the analysis  
        	- Columns are the features

    - **output**  
        - `outfile:file`: The file with the components  
        - `oudir:dir`   : The directory c  

    - **args**  
        - `devpars`: The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
        - `anopts` : The options to read the annotation files.  
        - `inopts` : The options to read the input files.  
        - `na`     : How to deal with `NA` values. Default: `0`  
        	- A logistic/boolean value will remove them (use `complete.cases`)
        	- Otherwise, it will be replaced by the given value.
        - `seed`   : The seed. Default: `None`  
        - `plots`  : Use R package `factoextra` to do Plots. You can use `False` for each to disable each plot.  
        	- `scree`  : Scree plot, see `?fviz_screeplot`. Default: `Box(ncp = 20`)
        	- `var`    : Var plot, see `?fviz_pca_var`. Default: `Box(repel = False)`
        	- `bi`     : Biplot,   see `?fviz_pca_biplot`. Default: `Box(repel = False)`
        	- `clplot` : Cluster plot, see `?fviz_cluster`. Default: `Box(repel = False, main = "", ggs = Box())`
        		- The extra `ggs` is used to extend the plot. See example in description.
        	- `cluster`: Cluster options for the cluster plot. Default: `Box(npcs  = .8, method = 'kmeans')`
        		- `npcs`: # of PCs to use for clustering. `npcs` < 1 will be treated as variance contribution of PCs. For example, `0.8` will take first N PCs will contribute 80% of variance. Default: `.8`
        		- `method`: Clustering method. Available methods would be `kmeans` and methods supported by `cluster` package.
        		- Other arguments for the clustering function.

    - **requires**  
        [`R-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html) for plots
## tabix

!!! hint "pTabix"

    - **description**  
        Use tabix to extract information.

    - **input**  
        - `infile`: a local or remote file  
        - `region`: a region or a file containing regions  

    - **output**  
        - `outfile:file`: The information extracted from the input file  

    - **args**  
        - `tabix`: The path to `tabix`  
        - `params`: Other params for `tabix`  

!!! hint "pTabixIndex"

    - **description**  
        Generate tabix index file.

    - **input**  
        - `infile:file`: the input file  
        	- Could be bgzipped.

    - **output**  
        - `outfile:file`: The bgzipped file  
        - `outidx:file`: The tabix index file  

    - **args**  
        - `tabix`: The path to `tabix`  
        - `params`: Other params for `tabix`  
## tcga

!!! hint "pDownload"

    - **description**  
        Download TCGA use `gdc-client` and a manifest file

    - **input**  
        - `manifile:file`: the manifest file  

    - **output**  
        - `outdir:file`: the directory containing downloaded file  

    - **args**  
        - `params`    : other params for `gdc-client download`, default: `{'no-file-md5sum': True}`  
        - `gdc_client`: the executable file of `gdc-client`,    default: "gdc-client"  
        - `nthread`   : Number of threads to use. Default     : `1`  
        - `token`     : The token file if needed.  

!!! hint "pSample2SubmitterID"

    - **description**  
        convert TCGA sample names with submitter id with metadata and sample containing folder

    - **input**  
        - `indir:file`: the directory containing the samples  
        - `mdfile:file`: the metadata file  

    - **output**  
        - `outdir:file`: the directory containing submitter-id named files  

    - **args**  
        - `method`: How the deal with the files. Default: `symlink`  
        	- We can also do `copy`
        - `nthread`: Number threads to use. Default: `1`  

!!! hint "pGtFiles2Mat"

    - **description**  
        Convert TCGA genotype files to a matrix.

    - **input**  
        - `infiles:files`: The input genotypes files  

    - **output**  
        - `outfile:file`: The output matrix file  

    - **args**  
        - `rsmap`  : The rsid probe mapping file. If not provided, will use the probe id for matrix rownames.  
        - `fn2sam` : How to convert filename(without extension) to sample name. Default: `None`  

!!! hint "pClinic2Survival"

    - **description**  
        Convert TCGA clinic data to survival data
        The clinic data should be downloaded as "BCR Biotab" format

    - **input**  
        - `infile:file`: The clinic data file downloaded from TCGA  

    - **output**  
        - `outfile:file`: The output file  
        - `covfile:file`: The covariate file  

    - **args**  
        - `cols`: The column names:  
        	- `time_lastfollow`: The column names of last follow up. Default: `['days_to_last_followup']`
        	- `time_death`: The column names of time to death. Default: `['days_to_death']`
        	- `status`: The columns of vital status. Default: `['vital_status']`
        	- `age`: The columns of days to birth. Default: `['days_to_birth']`
        - `covs`: The covariates to output. Default:  
        	- `gender`, `race`, `ethnicity`, `age`
        - `mat`: An expression or genotype matrix with samples as column names, used to get sample names for patient instead of short ones. Default: `None`  

!!! hint "pConvertExpFiles2Matrix"

    - **description**  
        convert TCGA expression files to expression matrix, and convert sample name to submitter id

    - **input**  
        - `dir:file`: the directory containing the samples  
        - `mdfile:file`: the metadata file  

    - **output**  
        - `outfile:file`: the output matrix  

    - **requires**  
        [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)

!!! hint "pConvertMutFiles2Matrix"

    - **description**  
        convert TCGA mutation files (vcf.gz) to mut matrix, and convert sample name to submitter id

    - **input**  
        - `dir:file`: the directory containing the samples  
        - `mdfile:file`: the metadata file  

    - **output**  
        - `outfile:file`: the output matrix  
## tfbs

!!! hint "pMotifScan"

    - **description**  
        Scan motif along the given sequences.

    - **input**  
        - `tffile:file`: The infile containing TF name and motif name.  
        	- If only one column is give, will be used as both TF and motif name
        	- If there are 2+ columns, 1st column will be motif name, 2nd column will be TF name
        - `sfile:file`: The sequence file  

    - **output**  
        - `outdir:file`: The output dir  

    - **args**  
        - `tools`   : The tool used to scan the motif. Default: 'meme'  
        - `meme`    : The path of MEME's fimo. Default: 'fimo'  
        - `motifs`  : The motif database in MEME format.  
        - `pval`    : The pvalue cutoff. Default: 1e-4  
        - `cleanmname`: Whether to clean motif name. Default: True  
        - `ucsclink`: The ucsc link template. Default: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}`  
        - `nthread` : Number of threads used to scan, only available when you have multiple mids. Default: 1  
        - `params`  : Other parameters for `fimo`  

    - **requires**  
        [`fimo` from MEME Suite](http://meme-suite.org/tools/fimo)

!!! hint "pAtSnp"

    - **description**  
        Scan motifs on Snps to detect binding affinity changes.

    - **input**  
        - `tffile:file`: The tf-motif file with 1st column the motif and 2nd the tf  
        - `snpfile:file`: The snp file.   
        	- Could be a bed file with first 6 columns the bed6 format and 7th the alleles.
        	- Alleles including reference allele should be seperated by `,`
        	- It also could be a snp file required by `atSNP` package.
        	- See: https://rdrr.io/github/kinsigne/atSNP_modified/man/LoadSNPData.html

    - **output**  
        - `outfile:file`: The output file  
        - `outdir:dir`  : The output directory containing the output file and plots.  

    - **args**  
        - `tfmotifs`: The motif database. Defaut: `params.tfmotifs`  
        - `genome`  : The reference genome to get the sequences. Default: `params.genome`  
        - `fdr`     : Do fdr or not. Or could be the p-value adjustment method. Default: `True` (using `BH` method)  
        - `pval`    : The pvalue cutoff for output and plot.  
        - `plot`    : Do plot or not. Default: `True`  
        - `nthread` : How many threads to use. Default: `1`  
        - `depvars` : The device parameters for plotting. Default: `Box(res = 300, width = 2000, height = 2000)`  

    - **requires**  
        `r-atSNP`
## tsv

!!! hint "pMatrixR"

    - **description**  
        Operate a matrix and save the new matrix to file.

    - **input**  
        - `infile:file`: The input file containing the matrix  

    - **output**  
        - `outfile:file`: The output matrix  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `cnames`: Whether the input file has cnames. Default: True
        	- `rnames  `: Whether the input file has rnames. Default: True
        	- `delimit`: The delimit. Default: `\t`
        	- `skip`: First N lines to skip. Default: `0`
        - `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`  
        - `code`: The R code to operating the matrix. (the matrix is read in variable `mat`)  

!!! hint "pTranspose"

    - **description**  
        Transpose a matrix

    - **input**  
        - `infile:file`: The input matrix file  

    - **output**  
        - `outfile:file`: The Transposed matrix file. Default: `{{i.infile | bn}}`  

    - **args**  
        - `inopts`: Input options for input file.  

!!! hint "pPaired"

    - **description**  
        Subset each input file and make sure they have paired columns.

    - **input**  
        - `infile1:file`: The first file  
        - `infile2:file`: The second file  

    - **outfile**  
        - `outfile1:file`: The paired file for infile1. Default: `{{i.infile1 | fn}}.paired{{i.infile1 | ext}}`  
        - `outfile2:file`: The paired file for infile2. Default: `{{i.infile2 | fn}}.paired{{i.infile2 | ext}}`  

    - **args**  
        - `inopts1`: reading options for input file1  
        - `inopts2`: reading options for input file2  

!!! hint "pCbind"

    - **description**  
        Cbind the rest of files to the first file.

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outfile:file`: The output matrix  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `cnames`: Whether the input file has cnames. Default: True
        	   - or [True, True, False] corresponding to the file order
        	- `rnames  `: Whether the input file has rnames. Default: True
        	- `delimit`: The delimit. Default: `\t`
        	- `skip`: First N lines to skip. Default: `0`
        - `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`  
        - `fn2cname`: The function (r) used to convert file name to column name.  
        - `fill`: Do `cbind.fill` instead of `cbind`. Default: `True`  
        	- Set it to `False` if the row names are in the same order
        - `na`: Replacement for missing values. Default: `NA`  

!!! hint "pRbind"

    - **description**  
        Rbind the rest of files to the first file.

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outfile:file`: The output matrix  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `cnames`: Whether the input file has cnames. Default: True
        	   - or [True, True, False] corresponding to the file order
        	- `rnames  `: Whether the input file has rnames. Default: True
        	- `delimit`: The delimit. Default: `\t`
        	- `skip`: First N lines to skip. Default: `0`
        - `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`  
        - `na`: Replacement for missing values. Default: `NA`  
        - `fn2rname`: The function (r) used to convert file name to row name.  
        - `fill`: Do `rbind.fill` instead of `rbind`. Default: `True`  
        	- Set it to `False` if the row names are in the same order
        - `na`: Replacement for missing values. Default: `NA`  

!!! hint "pCsplit"

    - **description**  
        Split a matrix by columns and save them into files.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outdir:dir`: The directory containing the output column files  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `cnames`: Whether the input file has cnames. Default: True
        	- `rnames  `: Whether the input file has rnames. Default: True
        	- `delimit`: The delimit. Default: `\t`
        	- `skip`: First N lines to skip. Default: `0`
        - `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`  
        - `size`: The chunk size (how many columns to split into one file). Default: `1`  

!!! hint "pRsplit"

    - **description**  
        Split a matrix by rows and save them into files.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outdir:dir`: The directory containing the output row files  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `cnames`: Whether the input file has cnames. Default: True
        	- `rnames  `: Whether the input file has rnames. Default: True
        	- `delimit`: The delimit. Default: `\t`
        	- `skip`: First N lines to skip. Default: `0`
        - `params`: Other params for `read.table`. Default: `{"check.names": "FALSE", "quote": ""}`  
        - `size`: The chunk size (how many rows to split into one file). Default: `1`  

!!! hint "pTsv"

    - **description**  
        Read, Transform, filter a TSV file.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `delimit`: The delimit. Default: `\t`
        	- `comment`: The comment sign. Default: `#`
        	- `skip`: First N lines to skip. Default: `0`
        	- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict). If not specified, metadata will be generated automatically.
        - `outopts`: The output options for outfile:  
        	- `delimit`: The delimit for records. Default: `\t`
        	- `head`: Output header or not. Default: `False`
        	- `headDelimit`: The delimit for header. Default: `\t`
        	- `headPrefix`: The prefix for header. Default: ``
        	- `headTransform`: The transformer for header. Default: `None`
        	- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict, '+' as an element or key is allowed to indicate extra meta from the reader). If not specified, metadata will be borrowed from the reader. 
        - `ops`: A ops function to transform the row. Argument is an instance of `readRecord`  
        - `opshelper`: A helper function for `args.ops`  

!!! hint "pTsvJoin"

    - **description**  
        Read files simultaneously.
        NOTE: only one file allows multiple lines with same value to compare, and that file should be the first one. For example: 
        ```
        File1:
        1	1
        1	2
        1	3
        File2:
        1	1
        2	2
        3	3
        ```
        If you compare the first column, File1 has to put at the begining for input.

    - **input**  
        - `infiles:files`: The input files  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `inopts`: The input options for infile:  
        	- `skip`   : First N lines to skip. Default: `0`
        	- `delimit`: The delimit. Default          : `\t`
        	- `comment`: The comment line mark. Default: `#`
        	- `cnames`   : Whether input file has head. Default: `True`
        - `outopts`:   
        	- `delimit`      : The delimit. Default: `\t`
        	- `cnames`         : Whether to output the head? Default: `False`
        - `match`: The match function.   
        - `do`: The do function. Global vaiable `fout` is available to write results to output file.  
        - `helper`: Some helper codes.  

    - **requires**  
        [`python-simread`](https://github.com/pwwang/simread)

!!! hint "pTsvSql"

    - **description**  
        Query tsv file using SQL. (see: http://harelba.github.io/q/examples.html)

    - **input**  
        - `infile:file` : The input tsv file  
        - `sqlfile:file`: The file containing the SQLs. If provided, `args.sql` will be ignored.   

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `sql`: If SQL to execute. Use `-` for table name  
        - `inopts`: Options for input file.  
        	- `cnames`: Input file has header? Default: `True`
        	- `delimit`: The delimit of input file. Default: `\t`
        	- `encoding`: Encoding of input file. Default: `UTF-8`
        	- `gz`: Whether input file is gzipped. Default: `auto` (detected from file extension)
        - `outopts`: Output options.  
        	- `cnames`: Inherited from `args.inopts`
        	- `delimit`: Inherited from `args.inopts`
        	- `encoding`: Inherited from `args.inopts`

    - **requires**  
        [`q`](http://harelba.github.io/q/index.html)
        This process is built on `q 1.7.1`

!!! hint "pMergeFiles"

    - **description**  
        Merge files in the input directory

    - **input**  
        - `indir:file`: The input directory  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `inopts`: The options for input file.  
        	- defaults: skip: 0, comment: #, delimit '\\t'
        - `outopts`: The options for output file. Defaults:  
        	- head: False (not output head line)
        	- headPrefix: `#` (The prefix for head line)
        	- headDelimit: `\\t` (The delimit for head line)
        	- headTransform: `None` (The callback for head line)
        	- delimit: `\\t` (The delimit for output line)

!!! hint "pMergeRows"

    - **description**  
        Merge repeated rows

    - **input**  
        - `infile:file`: The input file (has to be sorted by the repeated columns)  

    - **output**  
        - `outfile:file`: The output file. Default: `{{i.infile | bn}}`  

    - **args**  
        - `inopts`: The options for input file.  
        	- defaults: skip: 0, comment: #, delimit '\\t'
        - `outopts`: The options for output file. Defaults:  
        	- head: False (not output head line)
        	- headPrefix: `#` (The prefix for head line)
        	- headDelimit: `\\t` (The delimit for head line)
        	- headTransform: `None` (The callback for head line)
        	- delimit: `\\t` (The delimit for output line)
        - `match`: The function to return a value to decide whether the row is repeated, argument is a `TsvRecord`.  
        - `do`   : The merge function in python, argument is a list of `TsvRecord`s or `list`s if `args.inopts.ftype` is `nometa`  
## tumhet

!!! hint "pSciClone"

    - **description**  
        Run sciClone for subclonal analysis.

    - **input**  
        - `vfvcfs:files`: The VCF files of mutations of each sample  
        - `cnvcfs:files`: The VCF files of copy number variations of each sample  

    - **output**  
        - `outdir:dir`: The output directory.  

    - **args**  
        - `params`  : Other parameters for original `sciClone` function. Default: `Box()`  
        - `exfile`  : The regions to be excluded. In BED3 format  
        - `vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`  
        - `cnsamcol`: The index of the target sample in copy number VCF file, 1-based. Default: `1`  
        - `varcount`: An R function string to define how to get the variant allele count. Default: `function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])`  
        	- If this function returns `NULL`, record will be skipped.
        	- It can use the sample calls (`fmt`) and also the record info (`info`)
        	- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
        	- Don't include `info` if not necessary. This saves time.
        	- This function can return the variant count directly, or 
        	- an R `list` like: `list(count = <var count>, depth = <depth>)`.
        	- By default, the `depth` will be read from `fmt$DP`
        - `cncount` : An R function string to define how to get the copy number. Default: `function(fmt) fmt$CN`  
        	- Similar as `varcount`
        	- Returns copy number directly, or
        	- an R `list` like: `list(cn = <copy number>, end = <end>, probes = <probes>)`
        	- `end` defines where the copy number variation stops
        	- `probes` defines how many probes cover this copy number variantion.

!!! hint "pPyClone"

    - **description**  
        Run PyClone for subclonal analysis

    - **input**  
        - `vfvcfs:files`: The VCF files of mutations of each sample  
        - `cnvcfs:files`: The VCF files of copy number variations of each sample  

    - **output**  
        - `outdir:dir`: The output directory.  

    - **args**  
        - `params`  : Other parameters for original `PyClone run_analysis_pipeline` function. Default: `Box()`  
        - `vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`  
        - `cnsamcol`: The index of the target sample in copy number VCF file, 1-based. Default: `1`  
        - `varcount`: A python lambda string to define how to get the variant allele count. Default: `lambda fmt: fmt.get("AD") and fmt.get("AD")[1]`  
        	- If this function returns `None`, record will be skipped.
        	- It can use the sample calls (`fmt`) and also the record info (`info`)
        	- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
        	- This function can return the variant count directly, or 
        	- a `dict` like: `dict(count = <var count>, depth = <depth>)`.
        	- By default, the `depth` will be read from `fmt.DP`
        - `cncount` : An python lambda string to define how to get the copy number. Default: `lambda fmt: fmt.get("CN")`  
        	- Similar as `varcount`
        	- Returns copy number directly, or
        	- a `dict` like: `dict(cn = <copy number>, end = <end>)`
        	- `end` defines where the copy number variation stops

!!! hint "pQuantumClone"

    - **description**  
        Run QuantumClone: https://academic.oup.com/bioinformatics/article/34/11/1808/4802225

    - **input**  
        - `vfvcfs:files`: The input vcf files with mutations  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `params`  : other parameters for `QuantumClone`'s `One_step_clustering`  
        - `vfsamcol`: The index of the target sample in mutation VCF file, 1-based. Default: `1`  
        - `varcount`: An R function string to define how to get the variant allele count. Default: `function(fmt) as.integer(unlist(strsplit(fmt$AD, ","))[2])`  
        	- If this function returns `NULL`, record will be skipped.
        	- It can use the sample calls (`fmt`) and also the record info (`info`)
        	- Both `function(fmt) ...` and `function(fmt, info) ...` can be used.
        	- Don't include `info` if not necessary. This saves time.
        	- This function can return the variant count directly, or 
        	- an R `list` like: `list(count = <var count>, depth = <depth>)`.
        	- By default, the `depth` will be read from `fmt$DP`
        - `nthread` : # threads to use. Default: `1`  

!!! hint "pTheta"

    - **description**  
        Run THetA2 for tumor purity calculation
        Set lower MIN_FRAC if interval is not enough and NO_CLUSTERING if it raises 
        "No valid Copy Number Profiles exist", but have to pay attention to the results. 
        (see: https://groups.google.com/forum/#!topic/theta-users/igrEUol3sZo)

    - **args**  
        - `affysnps`: The affymetrix Array snps, or other candidate snp list, in BED6-like format  
        	- The first 6 columns should be in BED6 format
        	- The 7th column is reference allele, and 8th column is mutation allele.

    - **install**  
        `conda install -c bioconda theta2`
        `conda install -c bioconda bam-readcount`
## vcf

!!! hint "pVcfFilter"

    - **description**  
        Filter records in vcf file.

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `filters`: A dict of filters like `{<filtername>: <filter>}`. `<filter>` should be a string of lambda function:  
        	```
        	"lambda record, samples: <expression>"
        	* ``record.CHROM`` : 'chr20'
        	* ``record.POS``   : 1234567
        	* ``record.ID``    : 'microsat1'
        	* ``record.REF``   : ''GTC''
        	* ``record.ALT``   : [G, GTCT]
        	* ``record.QUAL``  : 50
        	* ``record.FILTER``: ['PASS'] # NO!, PASS should be []
        	* ``record.INFO``  : {'AA': 'G', 'NS': 3, 'DP': 9}
        	* samples = record.samples
        	* len(samples): 3
        	* samples[0].sample: 'NA00001'
        	* samples[0]: Call(sample=NA00001, CallData(GT=0/1, GQ=35, DP=4))
        	* samples[0].data: calldata(GT='0/1', GQ=35, DP=4)
        	* samples[0].data.GT: '0/1'
        	```
        	- see here for record and samples: https://github.com/jamescasbon/PyVCF
        	- Remember if filters() returns True, record filtered.
        	- For builtin filters, you may specify them as `{<filter>: <param>}`
        	- You can also use `!` to specify a negative builtin filter: `{!<filter>: <param>}`
        	- Bulitin filters: 
        		- SNPONLY: keeps only SNPs (`{"!SNPONLY": None}` means filtering SNPs out)
        		- BIALTONLY: keeps only mutations with bi-allele
        		- QUAL: keeps only site quality >= param (`{'QUAL': 30}`)
        - `gz`     : Whether to gzip the output file. Default: False  
        - `keep`   : Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)  

    - **requires**  
        [`pyvcf`](https://github.com/jamescasbon/PyVCF)

!!! hint "pVcfUnique"

    - **description**  
        Remove duplicate mutations from a VCF file.
        Because in most case where we want to remove the duplicate mutations, it might be 
        because other program not accepting them. In this case, we don't put a filter on 
        the records, but just remove them instead.

    - **input**  
        - `infile:file`: The input vcf file.  

    - **output**  
        - `outfile:file`: The output vcf file. Default: `{{i.infile | fn2}}.vcf{% if args.gz %}.gz{% endif %}`  

    - **args**  
        - `upart`: The unique parts. Could be part of: `['CHROM', 'POS', 'ID', 'REF', 'ALT']`  
        - `keep` : Which record to keep.  
        	- `bisnp` : Snp with bi-allele
        	- `snp`   : Snps
        	- `bialt` : Bi-allele mutations
        	- `first` : The first record (default)
        	- `last`  : The last record
        	- `random`: A random record
        	- Multiple ways can be used: `first, snp` is to select first snp (in case multiple snps duplicated)
        - `gz`: Bgzip the output vcf file or not. Default: `False`  

    - **requires**  
        `pyvcf`

!!! hint "pVcfRemoveFilter"

    - **description**  
        Remove one or more filters in vcf files

    - **input**  
        - `infile:file`: The input vcf file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `rmfilter`: The filters to remove. If None, ALL filters will be removed!  
        	- A `list` of filter names.

!!! hint "pVcf"

    - **description**  
        Use pyvcf to manipulate vcf file

    - **input**  
        - `infile:file`: The input vcf file  

    - **output**  
        - `outfile:file`: The output vcf file  

    - **args**  
        - `helper`: The helper code injected to script  
        	- Since lambda function can't do assignment and manipulation so you can write some help function here
        - `readerops`: A lambda function (must be quoted) to manipulate the reader (vcf.Reader instance)  
        - `recordops`: A lambda function (must be quoted) to manipulate the record (vcf.Record instance)  
        - `gz`: Gzip the ouput file  

!!! hint "pVcfAnno"

    - **description**  
        Annotate the variants in vcf file.
        You have to prepare the databases for each tool.

    - **input**  
        - `infile:file`: The input vcf file  

    - **output**  
        - `outfile:file`: The output file (output file of annovar will also be converted to vcf)  
        - `outdir`: The output directory, used to fetch some stat/summary files  

    - **args**  
        - `tool`: The tool used to do annotation. Default: snpeff  
        - `snpeff`: The path of snpeff. Default: snpEff  
        - `vep`: The path to vep. Default: vep  
        - `gz`: Whether to gzip the result file. Default: False  
        - `annovar`: The path of annovar. Default: annotate_variation.pl  
        - `annovar_convert`: The path of convert2annovar.pl, used to convert vcf to annovar input file. Default: convert2annovar.pl  
        - `genome`: The genome for annotation. Default: hg19  
        - `tmpdir`: The tmpdir, mainly used by snpeff. Default: <system tmpdir>  
        - `dbpath`: The path of database for each tool. Required by 'annovar' and 'vep'  
        - `params`: Other params for tool. Default: ''  
        - `snpeffStats`: Whether to generate stats file when use snpeff. Default: False  
        - `mem`: The memory used by snpeff. Default: '4G'  

    - **requires**  
        [`annovar`](http://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
        [`snpeff`](http://snpeff.sourceforge.net/SnpEff_manual.html#intro)
        [`vep`](http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html)

!!! hint "pVcfSplit"

    - **description**  
        Split multi-sample Vcf to single-sample Vcf files.

    - **input**  
        - `infile:file`: The input vcf file  
        - `samples`: The samples, if not provided, will extract all samples  

    - **output**  
        - `outdir:dir`: The output directory containing the extracted vcfs  

    - **args**  
        - `tool`: The tool used to do extraction. Default: vcftools (gatk, awk)  
        - `vcftools`: The path of vcftools' vcf-subset  
        - `bcftools`: The path of bcftools, used to extract the sample names from input vcf file.  
        - `gatk`: The path of gatk.  

!!! hint "pVcfMerge"

    - **description**  
        Merge single-sample Vcf files to multi-sample Vcf file.

    - **input**  
        - `infiles:files`: The input vcf files  
        - `outfile:dir`: The output multi-sample vcf.  

    - **args**  
        - `tool`: The tool used to do extraction. Default: vcftools  
        - `vcftools`: The path of vcftools' vcf-subset  
        - `bcftools`: The path of bcftools, used to extract the sample names from input vcf file.  
        - `gatk`: The path of gatk.  

!!! hint "pVcf2Maf"

    - **description**  
        Convert Vcf file to Maf file

    - **input**  
        - `infile:file` : The input vcf file  
        	- see `args.somatic`

    - **output**  
        - `outfile:file`: The output maf file  

    - **args**  
        - `tool`     : Which tool to use. Default: vcf2maf  
        - `vcf2maf`  : The path of vcf2maf.pl  
        - `vep`      : The path of vep  
        - `vepDb`    : The path of database for vep  
        - `filtervcf`: The filter vcf. Something like: ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz  
        - `ref`      : The reference genome  
        - `nthread`  : Number of threads used to extract samples. Default: 1  
        - `tumor1st` : Whether tumor sample comes first. Default: `True`  
        - `bcftools` : Path to bcftools used to extract sample names.  
        - `vcftools` : Path to vcftools used to split vcf.  
        - `samfunc`  : A lambda function used to deduce sample names from file name.  
        - `somatic`  : Whether input vcf file is a somatic mutation file. Default: False  
        	- somatic mutation vcf file can only have one sample TUMOR, or two samples, TUMOR and NORMAL, but will be considered as single sample.
        	- otherwise, multiple samples are supported in the input vcf file. Tumor id will be sample name for each sample, normal id will be NORMAL.

!!! hint "pVcf2Plink"

    - **description**  
        Convert vcf to plink binary files (.bed/.bim/.fam)

    - **input**  
        - `infile:file`: The input vcf file, needs to be tabix indexed.  

    - **output**  
        - `outdir:dir`: The output directory containing the plink binary files  

    - **args**  
        - `plink` : The path to plink  
        - `params`: Command arguments for `plink`. Some pre-settings:  
        	- `vcf-half-call`      : `m`
        	- `double-id`          : `True`
        	- `vcf-filter`         : `True`
        	- `vcf-idspace-to`     : `_`
        	- `set-missing-var-ids`: `@_#`    # make sure no duplicate vars
        	- `biallelic-only`     : `strict`

!!! hint "pVcfLiftover"

    - **description**  
        Lift over vcf files.

    - **input**  
        - `infile:file`: The input vcf file.  

    - **output**  
        - `outfile:file`: The output vcf file.  
        - `umfile:file`: The unmapped records  

    - **args**  
        - `tool`: Which tool to use. Default: `picard`  
        - `picard`: The path to picard  
        - `lochain`: The liftover chain file  
        - `ref`: The reference genome  
        - `mem`: The memory to use  
        - `tmpdir`: The temporary directory  
        - `params`: The extra params.  

!!! hint "pVcfAddChr"

    - **description**  
        Add `chr` to records of vcf files.

    - **args**  
        - `chr`: The prefix to add to each record.  

!!! hint "pVcfCleanup"

    - **description**  
        Remove configs from vcf file according to the given reference.

    - **input**  
        - `infile:file`: The input vcf file  

    - **output**  
        - `outfile:file`: The output vcf file. Default: `{{i.infile | fn}}.vcf`  

    - **args**  
        - `ref`: The reference file  

!!! hint "pVcf2GTMat"

    - **description**  
        Convert Vcf file to genotype matrix.
        Rownames are in the format of '<chr>_<pos>_<rs>_<ref>_<alt>'

    - **input**  
        - `infile:file`: The input vcf file  

    - **output**  
        - `outfile:file`: the output filename. Default: `{{i.infile | fn2}}.gtmat`  

    - **args**  
        - `novel`: The snp name used if not mapped to any rsid. Default: `NOVEL`  
        - `useid`: Use the id in vcf file is possible. Default: `True`  
        - `dbsnp`: The dbsnp vcf file used to get the rsid. If not provided, will use `novel`  
        - `na`   : The value to replace missing genotypes.  
        - `bialt`: bi-allelic snps only. Default: `True`  

    - **requires**  
        `pytabix`
        `pysam`

!!! hint "pVcfSort"

    - **description**  
        Sort the vcf records

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output file  

    - **args**  
        - `header`: Output header? Default: `True`  
        - `by`    : Sort by what, Coordinates (coord) or names (name)? Default: `coord`  
        - `tool`  : The tool used to do the sort. Default: `sort` (linux command)  

!!! hint "pVcfSubtract"

    - **description**  
        Subtract one vcf file from another

    - **input**  
        - `infile1:file`: The vcf file to be subtracted  
        - `infile2:file`: The background vcf file  

    - **output**  
        - `outfile:file`: The subtracted vcf file.  

    - **args**  
        - `header`  : Output header? Default: `True`  
        - `bychrom` : Split the vcf file by chromosomes, do subtraction and then merge them. Default: `False`  
        	- In case the vcf file is too big. 
        	- Requires both vcf files indexed (.tbi). If not they will be indexed there.
        - `nthread` : # threads to use, only when `bychrom` is True. Default: `1`  
        - `tool`    : The tool to be used. Default: `mem` (or pyvcf/bedtools)  
        - `bedtools`: The path to bedtools.  
        - `tabix`   : The path to tabix.  
        - `any`     : Remove record in `infile1` with any overlap in `infile2`. Default: `True`  

!!! hint "pVcfExtract"

    - **description**  
        Extract variants from a VCF file by given regions

    - **args**  
        - `tabix` : The path to tabix.  
        - `params`: Other parameters for `tabix`. Default: `Box(h = True, B = True)`  
        	- See `tabix --help`
## vcfnext

!!! hint "pVcfStatsPlot"

    - **description**  
        Convert csvstat file from snpEff to R-readable matrix and plot them.

    - **input**  
        - `indir:file`: The directory containing the csv stat files from `snpEff ann`  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `chroms`: The chromsome filter. Default: "" (all chroms)  
        - Note: snpEff csvstat file has no "chr" prefix

!!! hint "pGTMatAddRs"

    - **description**  
        Add rs id to a genotype matrix

    - **input**  
        - `infile:file`: The input genotype matrix, columns are samples, rows are mutations in format:  
        	- `<chr>_<pos>_<ref>_<alt>` or `<chr>_<pos>_<name>_<ref>_<alt>`
        	- has to be sorted (see `args.chrsort`)

    - **output**  
        - `outfile:file`: The output genotype matrix. Row names will turn into:  
        	- `<chr>_<pos>_<rs>_<ref>_<alt>`

    - **args**  
        - `dbsnp`: the dbsnp vcf file used to annotation the snps.  
        	- assume sorted by coordinates
        	- see `args.chrsort` for chromsome sorting
        - `notfound`: What to used if RS id not found. Default: `NOVEL`  
        - `chrsort`: How chromsome is sorted. Default: `version`  
        	- version sort: chr1, chr2, chr3, ..., chr10, chr11, ... chr20, ..., chrM(T), chrX, chrY
        	- natural sort: chr1, chr10, chr11, ..., chr19, chr2, chr21, chr22, chr3, ..., chrM(T), chrX, chrY
        	- or a list of chromsome order: e.g: `["chr1", "chr2", ..., "chrM", "chrX", "chrY"]`

!!! hint "pGTMat2Plink"

    - **description**  
        Convert a genotype matrix to plink binary files

    - **input**  
        - `infile:file`: The genotype matrix, probably generated by `pVcf2GTMat`  
        - `metafile:file`: The metadata file.  
        	- column names could be `['FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno']`, see plink's `ped` format
        	- row names are samples

    - **output**  
        - `outdir:dir`: The output directory. Default: `{{i.infile | fn}}.plink`  

    - **args**  
        - `plink`: The path to `plink`  
        - `keeptxt`: Keep the text files (.ped and .map) or not. Default: `False`  

    - **requires**  
        `plink 1.x`

!!! hint "pGTMat2Bed"

    - **description**  
        Convert a genotype matrix to a bed file containing the coordinates of the mutations

    - **input**  
        - `infile:file`: The genotype matrix. Row names must follow `<chr>_<pos>_<rsid>_<ref>_<alt>`  

    - **output**  
        - `outfile:file`: The output bed file. Default: `outfile:file:{{i.infile | fn}}.bed`  

    - **args**  
        - `ncol`: How many columns of bed to output. Default: `6`. Possible values: 3, 6 and 8  
        - `name`: Use the neat name (usually rsid) or full name (row names). Default: `neat`  

!!! hint "pCallRate"

    - **description**  
        Calculate sample/snp call rate from single sample vcfs

    - **input**  
        - `indir:file`: The dir containing the vcfs  

    - **output**  
        - `outsample:file`: The report of call rate for each sample  
        - `figsample:file`: The bar chat of sample call rates  
        - `outsnp:file`: The report of call rate for each snp  
        - `figsnp:file`: The bar chat of snp call rates  

!!! hint "pCepip"

    - **description**  
        Run CEPIP.

    - **input**  
        - `infile:file`: The input file (vcf or avinput)  

    - **output**  
        - `outfile:file`: The cepip result file  

    - **args**  
        - `cepip`: The path of cepip  
        - `cell` : The related cell line  
        - `params`: Other params for cepip  

    - **requires**  
        [`cepip`](http://jjwanglab.org/cepip/)

!!! hint "pMutSig"

    - **description**  
        MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.
        For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
        
        See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)

    - **input**  
        - `infile:file`: mutation table  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `mutsig` : The path to `run_MutSigCV.sh`, default: 'mutsig'  
        - `mcr`    : The Matlab MCR path  
        - `cvrg`   : coverage table  
        - `cvrt`   : covariates table  
        - `mutdict`: mutation_type_dictionary_file  
        - `chrdir` : chr_files_hg18 or chr_files_hg19  

    - **requires**  
        [MutSig](http://archive.broadinstitute.org/cancer/cga/mutsig_download)

!!! hint "pMafLiftover"

    - **description**  
        Liftover maf file from one assembly to another

    - **input**  
        - `infile:file`: The input maf file  

    - **output**  
        - `outfile:file`: The output maf file  

    - **args**  
        - `liftover`: The liftOver program.  
        - `lochain`: The liftOver chain file.  
        - `genome`: The target genome.  

    - **requires**  
        liftOver from UCSC

!!! hint "pMafMerge"

    - **description**  
        Merge maf files.

    - **input**  
        - `infiles:files`: The maf files  

    - **output**  
        - `outfile:file`: The merged maf file  

    - **args**  
        - `excols`: How to deal with extra columns other than 34 standard columns from TCGA.  
        	- merge(default): Merge the columns, if one not exists, fill with an empty string.
        	- discard: Just discard the extra columns, with only 34 columns left. So you can also put just one maf file in the indir with some columns missed to fill it with standard columns.

!!! hint "pMaf2Mat"

    - **description**  
        Convert maf file to a gene(row)-sample(column) matrix

    - **input**  
        - `infile:file`: The input file  

    - **output**  
        - `outfile:file`: The output matrix  

    - **args**  
        - `mutypes`: Provide manual list of variant classifications to be counted, only effective when `args.binary = False`. Default: `None` (all counted)  
        - `binary` : Just generate a binary matrix instead of a count matrix. Default: `False`  
        - `na`: What value to use for no mutations reported on a gene. Default: `0`  
        - `samfn`  : A function (in r) to transform the sample names. Default: `function(sample) sample`  

!!! hint "pMaftools"

    - **description**  
        Use maftools to draw plots.

    - **input**  
        - `indir:dir`: The input directory. Could contain:  
        	- `*.maf` or `*.maf.gz` file (required)
        	- `*.annot.tsv` or `*.annot.txt` file (see: https://github.com/PoisonAlien/maftools/blob/master/inst/extdata/tcga_laml_annot.tsv)
        	- `all_lesions.conf_*.txt`: Gistic cnv data
        	- `amp_genes.conf_*.txt`: Gistic cnv data
        	- `del_genes.conf_*.txt`: Gistic cnv data
        	- `scores.gistic`: Gistic cnv data
        	- `*.seg.txt`: CBS segments data
        	- `*sig_genes.txt` or `*sig_genes.txt.gz`: Mutsig results, to do pancancer somparison.

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `ngenes` : Top number of genes to plot for some plots. Default: `10`  
        - `mutypes`: Provide manual list of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants. Default: `["Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"]`  
        - `isTCGA`: If the maf file is from TCGA? Default: `False`  
        - `ref`   : The reference file for signature plot.  
        - `plot`  : Which plots to plot.   
        	- Default:
        	  ```python
        	  Box(
        	  	summary        = True,
        	  	oncoplot       = True,
        	  	oncostrip      = True,
        	  	titv           = True,
        	  	lollipop       = True,
        	  	cbsseg         = True,
        	  	rainfall       = True,
        	  	tcgacomp       = True,
        	  	vaf            = True,
        	  	genecloud      = True,
        	  	gisticGenome   = True,
        	  	gisticBubble   = True,
        	  	gisticOncoplot = True,
        	  	somInteraction = True,
        	  	oncodrive      = True,
        	  	pfam           = True,
        	  	pancan         = True,
        	  	survival       = True,
        	  	heterogeneity  = True,
        	  	signature      = True,
        	  )
        	  ```
        - `params`: The extra parameters for each plot function.  
        	- Default:
        	  ```python
        	  Box(
        	  	summary        = Box(rmOutlier = True, addStat = 'median', dashboard = True),
        	  	oncoplot       = Box(),
        	  	oncostrip      = Box(),
        	  	titv           = Box(),
        	  	lollipop       = Box(AACol = 'Protein_Change'),
        	  	cbsseg         = Box(labelAll = True),
        	  	rainfall       = Box(detectChangePoints = True),
        	  	tcgacomp       = Box(),
        	  	vaf            = Box(flip = True),
        	  	genecloud      = Box(minMut = 3),
        	  	gisticGenome   = Box(markBands = 'all'),
        	  	gisticBubble   = Box(),
        	  	gisticOncoplot = Box(),
        	  	somInteraction = Box(),
        	  	oncodrive      = Box(AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore', fdrCutOff = 0.1, useFraction = True),
        	  	pfam           = Box(AACol = 'Protein_Change'),
        	  	pancan         = Box(qval = 0.1, label = 1, normSampleSize = True),
        	  	survival       = Box(),
        	  	heterogeneity  = Box(),
        	  	signature      = Box(nTry = 6, plotBestFitRes = False),
        	  )
        	  ```
        - `devpars`: The parameters for plot device. Default: `Box(res = 300, height = 2000, width = 2000)`  
        - `nthread`: Number of threads used for multiple plot of one type. Default: `1`  

    - **requires**  
        [Maftools](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)

!!! hint "pMutationSigs"

    - **description**  
        Find similar COSMIC mutation signatures for MAF file 
        using https://github.com/pwwang/deconstruct_sigs_py

    - **input**  
        - `infile:file`: The input maf file.  

    - **output**  
        - `outdir:dir`: The output directory  

    - **args**  
        - `font_family`: Font family for plotting.   
        - `font_weight`: Font weight for plotting.   
        - `sig_cutoff` : Significance cutoff for signatures.   
        - `err_thres`  : The threshold to top the iteration.  
        - `ref`        : The reference genome.  

!!! hint "pSnpEff"

    - **description**  
        This is the default command. It is used for annotating variant filed (e.g. VCF files).

    - **input**  
        - `infile:file`: The input file   

    - **output**  
        - `outdir:file`: The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html  

    - **args**  
        - `snpEff`: The snpEff executable, default: "snpEff"  
        - `params`: Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"  
        - `genome`: The genome used for annotation, default: "hg19"  
        - `informat`: The format of input file [vcf or bed], default: "vcf"  
        - `outformat`: The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"  
        - `csvStats`: Whether to generate csv stats file, default: True.  
        - `htmlStats`: Whether to generate the html summary file, default: False.  
        - `javamem`: The memory to use. Default: '-Xms1g -Xmx8g'  

    - **requires**  
        [snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)
## web

!!! hint "pDownloadForm"

    - **description**  
        Download results by submitting a form, supporting pagination.

    - **input**  
        - `url`   : the URL contains the form  
        - `data`  : the data used to fill the form (JSON string or transformed from dict by json.dumps).  
        - `submit`: the submit button to submit the form (use Xpath).  
        - `next`  : the button for next page (use Xpath)  

    - **output**  
        - `outdir:file`: The directory saves the results  

    - **args**  
        - `interval`: seconds to wait between fetching each page. Default: 1  

    - **requires**  
        [`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
        [`Phantomjs`](http://phantomjs.org/)

!!! hint "pDownloadGet"

    - **description**  
        Download results by urls.

    - **input**  
        - `url`: the URLs to download  

    - **output**  
        - `outfile:file`: The output file  

!!! hint "pDownload"

    - **description**  
        Alias of `pDownloadGet`

!!! hint "pDownloadPost"

    - **description**  
        Download results by POST.

    - **input**  
        - `url` : the URLs to download  
        - `data`: the POST data.  

    - **output**  
        - `outfile:file`: The output file  
## xlsx

!!! hint "pTsvs2Xlsx"

    - **description**  
        Save tsv files to xlsx sheets.

    - **input**  
        - `infiles:files`: The input tsv files  

    - **output**  
        - `outfile:file`: The output xlsx file  

    - **args**  
        - `fn2sheet`: How to convert filename(without extension) to sheet name  

    - **requires**  
        python packages: `csv` and `openpyxl`
