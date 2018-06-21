
## Main command
```shell
> bioprocs

Usage: bioprocs <command> [options ...]

Method commands:
  params                        - Show bioprocs.params items.
  list                          - List the process
                                  Use "bioprocs list <MOD>" to list all processes of a module.
  help                          - Print this help information.
                                  Use "bioprocs help <CMD>" to show help for commands.

Script commands:
  genomeplot                    - Plot genomic elements.
  mutcall                       - Call mutations from sequencing data.
  pwmscan                       
  runcmds                       - Using PyPPL to distribute and run commands.

Process commands:
  <module.proc>                 - Use "{prog} list" to show all processes, or
                                  use "{prog} list <module>" to show processes of a module.
```

Run `bioprocs help <command>` to show the help information of the subcommand.

## Subcommand: params
```shell
> bioprocs help params

Show bioprocs.params items.

Usage:
  bioprocs params
  bioprocs params [PARAM1] [PARAM2 ...]
  bioprocs params [-param PARAM1 PARAM2 ...]

Options:
  param                         - The parameter names
```

## Subcommand: list
```shell
> bioprocs help list

List the process
Use "bioprocs list <MOD>" to list all processes of a module.

Usage:
  bioprocs list
  bioprocs list [MOD]
  bioprocs list [-mod MOD]

Options:
  mod                           - The name of the module.
```

## Script: genomeplot
```shell
> bioprocs help genomeplot
USAGE:
  bioprocs -inputs <LIST> -tracks <LIST> -outdir <STR> -region <STR> -names <LIST> [OPTIONS]

REQUIRED OPTIONS:
  -inputs     <LIST>                    - The input of the tracks.
                                          "<ucscTrack>:<gvizTrack>" for ucsc track;
                                          "<infile>:<intype>" for interaction tracks;
                                          files for other tracks.
  -tracks     <LIST>                    - The track types. Could be data, anno, interaction or ucsc, or multiple of them.
  -outdir     <STR>                     - Output directory.
  -region     <STR>                     - The region to plot. E.g. chr1:1000000-1200000
  -names      <LIST>                    - The corresponding names of the tracks.

OPTIONAL OPTIONS:
  -ideo       <STR>                     - Show ideogram track? DEFAULT: 'True'
  -plotparams <STR>                     - The params for pGenomePlot. DEFAULT: ''
  -highlights <LIST>                    - The highlight regions in format of "start-end"
                                          DEFAULT: []
  -genes      <STR>                     - Show gene track?
                                          DEFAULT: '/data2/junwenwang/shared/reference/hg19/hg19-gene.gtf'
  -devpars    <STR>                     - The device parameters for plotting.
                                          DEFAULT: '{"res": 300, "height": 300, "width": 2000}'
  -params     <LIST>                    - The params for each track DEFAULT: []
  -splitlen   <INT>                     - Split the plot into 2 if the region is longer then splitlen.
                                          DEFAULT: 1000000
  -genome     <STR>                     - Commonly used genome assembly. DEFAULT: 'hg19'
  -forks      <INT>                     - Number of cores used to plot if split. DEFAULT: 2
  -axis       (BOOL)                    - Show axis? DEFAULT: True
  -h, --help, -H, -?                    - Print this help information.
```

## Script: mutcall
```shell
> bioprocs help mutcall
USAGE:
  bioprocs -exdir <STR> -saminfo <STR> -indir <STR> [OPTIONS]

REQUIRED OPTIONS:
  -exdir    <STR>                       - Where to export the result files.
  -saminfo  <STR>                       - The sample information file:
                                          Column 1: the basename of the sample file in '-indir'
                                          Column 2: Group if '-muts' includes 'soma', basically, 'TUMOR' and 'NORMAL'
                                          Column 3: Patient if multiple samples belong to the same patient.
                                          Example:
                                          Sample	Group	Patient
                                          A1.bam	TUMOR	A
                                          A2.bam	NORMAL	A
                                          B1.bam	TUMOR	A
                                          B2.bam	NORMAL	A
  -indir    <STR>                       - The input directory containing input files.

OPTIONAL OPTIONS:
  -runner   <STR>                       - The runner to run the processes. DEFAULT: 'local'
  -intype   <STR>                       - The input file types. Either bam, ebam or fastq.
                                          Ebam means the bam files need to reformatted into fastq files.
                                          DEFAULT: 'bam'
  -compress (BOOL)                      - Use gzip and bam file to save space. DEFAULT: True
  -muts     <LIST>                      - What kind of mutations to call.
                                          Note: soma need paired information DEFAULT: ['germ']
  -forks    <INT>                       - How many jobs to run simultanuously. DEFAULT: 1
  -logfile  <STR>                       - Where to save the logs. DEFAULT: ''
  -h, --help, -H, -?                    - Print this help information.
```

## Script: runcmds
```shell
> bioprocs help runcmds
USAGE:
  bioprocs -cmds <STR> [OPTIONS]

REQUIRED OPTIONS:
  -cmds   <STR>                         - The cmd list. If not provided, STDIN will be used.

OPTIONAL OPTIONS:
  -runner <STR>                         - The runner. DEFAULT: 'local'
  -forks  <INT>                         - How many jobs to run simultaneously. DEFAULT: 1
  -h, --help, -H, -?                    - Print this help information.
```

## Process

Convert a `pyppl.Proc` to a command line interface.  
Try `bioprocs <module>.<process>` to show the help information of the process.
