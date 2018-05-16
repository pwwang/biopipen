# genomeplot
<!-- toc -->
{% raw %}

## pInteractionTrack

### description
   Gererate genomic interaction track for Gviz

### input
   `name`: The name of the track
   `infile:file`: The input file. 
       - See the `type` argument for `makeGenomicInteractionsFromFile` from `GenomicInteractions` r-package
   `region`: the region, just chromosome!

### output
   `outfile:file`: The dumped track data

## pGeneTrack

### description
   Generate gene track using ucsc data source

### input
   `name`:   The name of the track

### output
   `outfile:file`: The file to save the track

### args
   `genome`: The genome
   `params`: use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args

## pAnnoTrack

### description
   The annotation track of Gviz

### input
   `name`:        the name of the track
   `infile:file`: the file for the track. (wig, bigWig or bedGraph, bam, need to be indexed!)
   `chrom`:       the chrom

### output
   `outfile:file`:the rds file for the track

### args
   `genome`: The genome
   `params`:  See `displayPars(DataTrack())` for all available display params

## pDataTrack

### description
   The data track of Gviz

### input
   `name`:        the name of the track
   `infile:file`: the file for the track. (wig, bigWig or bedGraph, bam, need to be indexed!)
   `chrom`:       the chrom

### output
   `outfile:file`:the rds file for the track

### args
   `genome`: The genome
   `params`:  See `displayPars(DataTrack())` for all available display params

## pUcscTrack

### description
   Generate track from ucsc

### input
   `name`     : the name of the track
   `track`    : the UCSC track
   `trackType`: the Gviz track
   `region`   : the region

### output
   `outfile:file`:the dumped track

### args
   use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.

## pGenomePlot

### description
   plot the genomic features

### input
   `trkfiles:files`: the list of track dumped files
   `region`:         the region, in format of `chr1:1-1000`
   `highlight`:      the highlight regions, informat of start1-end1; start2-end2; ...

### output
   `outfile:file`:   the figure

### args
   `genome`  : The genome
   `showIdeo`: Show ideogram track? Default: True
   `showAxis`: Show axis? Default: True
   `showGenes`: Show geneTrack? Default: True
   `params`:   The params
       - `genneral`:  General params for plotTracks
       - `geneTrack`: The params for geneTrack
{% endraw %}
