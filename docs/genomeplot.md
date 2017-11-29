<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pGeneTrack](#pgenetrack)
- [pAnnoTrack](#pannotrack)
- [pDataTrack](#pdatatrack)
- [pUcscTrack](#pucsctrack)
- [pGenomePlot](#pgenomeplot)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pGeneTrack

### description
    Generate the gene track using ucsc data source

### input
    `name`:   The name of the track
    `genome`: The genome
    `chrom`:  The chromosome
    `from`:   The start
    `to`:     The end

### output
    `outfile:file`: The file to save the track

### args
    use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args

## pAnnoTrack

### description
    Generate annotation track

### input
    `infile:file`: the file for the track
    `name`:        the name of the track
    `genome`:      the genome
    `chrom`:       the chromosome
    `from`:        the start position to display
    `to`:          the end position to display

### output
    `outfile:file`:the dumped track

### args
    use `displayPars(AnnotationTrack())` to see all available args.

## pDataTrack

### description
    The data track of Gviz

### input
    `infile:file`: the file for the track
    `name`:        the name of the track
    `genome`:      the genome
    `chrom`:       the chromosome
    `from`:        the start position to display
    `to`:          the end position to display

### output
    `outfile:file`:the rds file for the track
    `gout`:        the genome
    `cout`:        the chromosome
    `fout`:        the start
    `tout`:        the end

### args
    See `displayPars(DataTrack())` for all available display params
    Quote all params!

## pUcscTrack

### description
    Generate track from ucsc

### input
    `ucscTrack`:   the track to fetch from ucsc. [Avialable tracks](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
    `table`:       the table from ucsc. [Available table](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
    `gvizTrack`:   the object track to generate. One of "AnnotationTrack", "GeneRegionTrack", "DataTrack", "GenomeAxisTrack"
    `name`:        the name of the track
    `genome`:      the genome
    `chrom`:       the chromosome
    `from`:        the start position to display
    `to`:          the end position to display

### output
    `outfile:file`:the dumped track

### args
    use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.

## pGenomePlot

### description
    plot the genomic features

### input
    `trkfiles:files`: the list of track dumped files
    `genome`:         the genome
    `chrom`:          the chromosome
    `from`:           the start position to display
    `to`:             the end position to display

### output
    `outfile:file`:   the figure
