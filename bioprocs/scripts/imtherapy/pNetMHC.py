from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell, logger
from bioprocs.utils.parallel import Parallel, distributeList

{%from os import path%}
{%from pyppl.utils import always_list%}
infile    = {{i.infile | quote}}
afile     = {{i.afile | ?path.isfile | =readlines | !always_list | repr}}
outfile   = Path({{o.outfile | quote}})
allfile   = {{o.outfile | prefix | @append: '.all' | @append: ext(o.outfile) | quote}}
netmhc    = {{args.netmhc | quote}}
isfa      = {{args.isfa | repr}}
nthread   = {{args.nthread | repr}}
params    = {{args.params | repr}}
tmpdir    = {{args.tmpdir | repr}}
lens      = {{args.lens | ?isinstance: list | =@join: ',' | quote}}

shell.load_config(netmhc = netmhc)

# support HLA-A*03:79 -> HLA-A0379
alleles = [allele.strip().replace('*', '').replace(':', '') for allele in afile if 'HLA-' in allele]
valid_alleles = shell.netmhc(listMHC = True).splitlines()

for i in range(nthread):
	shell.mkdir(p = outfile.parent.joinpath('threads', str(i+1)))

# split infile
if isfa:
	seqs = [line.strip() for line in shell.grep('>', infile).splitlines() if line.strip()]
	seqs_to_threads = distributeList(seqs, nthread)
	seqs = {}
	for i, tseqs in enumerate(seqs_to_threads):
		for tseq in tseqs:
			seqs[tseq] = i
	handlers = {}
	lastindex = None
	with open(infile) as fin:
		for line in fin:
			if line.startswith('>'):
				seq = line.strip()
				index = seqs[seq]
				if index not in handlers:
					handlers[index] = open(outfile.parent.joinpath('threads', str(index+1), 'peptides.txt'), 'w')
				handlers[index].write(line)
				lastindex = index
			elif lastindex is None:
				raise IndexError('Sequence tag not found!')
			else:
				handlers[lastindex].write(line)
	for handler in handlers.values():
		if not handler.closed:
			handler.close()
else:
	with open(infile) as fin:
		peptides = fin.readlines()
	pep_to_threads = distributeList(peptides, threads)
	for i, pep in enumerate(pep_to_threads):
		with open(outfile.parent.joinpath('threads', str(i+1), 'peptides.txt'), 'w') as fpep:
			fpep.write(''.join(pep))

"""
	PARAMETER            DEFAULT VALUE        DESCRIPTION
	[-a filename]        HLA-A0201            HLA allele name
	[-f filename]                             Input file (by default in FASTA format)
	[-p]                 0                    Switch on if input is a list of peptides (Peptide format)
	[-l string]          9                    Peptide length (multiple lengths separated by comma e.g. 8,9,10)
	[-s]                 0                    Sort output on decreasing affinity
	[-rth float]         0.500000             Threshold for high binding peptides (%Rank)
	[-rlt float]         2.000000             Threshold for low binding peptides (%Rank)
	[-listMHC]           0                    Print list of alleles included in netMHC
	[-xls]               0                    Save output to xls file
	[-xlsfile filename]  NetMHC_out.xls       File name for xls output
	[-t float]           -99.900002           Threshold for output
	[-thrfmt filename]   $NETMHC/data/threshold/%s.thr Format for threshold filenames
	[-hlalist filename]  $NETMHC/data/allelelist File with covered HLA names
	[-rdir filename]     $NETMHC              Home directory for NetMHC
	[-tdir filename]     $TMPDIR              Temporary directory (Default $$)
	[-syn filename]      $NETMHC/data/synlists/%s.synlist Format of synlist file
	[-v]                 0                    Verbose mode
	[-dirty]             0                    Dirty mode, leave tmp dir+files
	[-inptype int]       0                    Input type [0] FASTA [1] Peptide
	[-version filename]  $NETMHC/data/version File with version information
	[-w]                 0                    w option for webface
"""
# common options
params.tdir = tmpdir
params.l    = lens

def do_one(allele, ifile, ithread):
	ps        = params.copy()
	ps.p      = not isfa
	ps.f      = ifile
	ps.a      = allele
	ps._out   = outfile.parent.joinpath('threads', str(ithread+1), allele + '.out.txt')
	ps._debug = True
	shell.netmhc(**ps)

args = []
for allele in alleles:
	if allele not in valid_alleles:
		logger.warning('Not a valid allele: %s', allele)
	for i in range(nthread):
		if outfile.parent.joinpath('threads', str(i+1), 'peptides.txt').is_file():
			args.append((allele, outfile.parent.joinpath('threads', str(i+1), 'peptides.txt'), i))

if not args:
	raise ValueError('No valid alleles found.')

para = Parallel(nthread = nthread)
para.run(do_one, args)

# merge results
with open(outfile, 'w') as fout, open(allfile, 'w') as fall:
	header_written = False
	pos_written = False
	for i, ofile in enumerate(outfile.parent.joinpath('threads').glob('*/*.out.txt')):
		with open(ofile) as fo:
			for line in fo:
				line = line.strip()
				if not line or line.startswith('-'):
					continue
				if header_written and line.startswith('#'):
					continue
				if i == 0 and line.startswith('#'):
					fout.write(line + '\n')
					fall.write(line + '\n')
				else:
					header_written = True
					parts = line.split()
					if parts and parts[0] == 'pos' and i == 0 and not pos_written:
						fout.write('\t'.join(parts) + '\n')
						fall.write('\t'.join(parts) + '\n')
						pos_written = True
					elif not parts or parts[0] in ('pos', 'Protein'):
						continue
					elif len(parts) > 14:
						del parts[-2]
						fout.write('\t'.join(parts) + '\n')
						fall.write('\t'.join(parts) + '\n')
					else:
						fall.write('\t'.join(parts) + '\n')
