import random
from pathlib import Path
from diot import Diot, OrderedDiot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader
from bioprocs.utils.parallel import Parallel
from bioprocs.utils.reference import faIndex

infile   = Path({{i.infile | quote}})
ref      = {{args.ref | quote}}
outfile  = Path({{o.outfile | quote}})
nmuts    = {{args.nmuts | repr}}
params   = {{args.params | repr}}
vcfy     = {{args.vcfy | quote}}
nthread  = {{args.nthread | repr}}
bcftools = {{args.bcftools | repr}}
bedtools = {{args.bedtools | quote}}
tmpdir   = {{args.tmpdir | quote}}
samples  = {{args.samples | repr}}
gz       = {{args.gz | repr}}

shell.load_config(vcfy = vcfy, bcftools = bcftools, bedtools = bedtools)

mergedfile = outfile.parent.joinpath(infile.stem + '.merged.bed')
shell.bedtools.merge(i = infile).r > mergedfile

# get total length of regions
total = 0
reader = TsvReader(mergedfile, cnames = False)
for r in reader:
	total += int(r[2]) - int(r[1])

mutrate = float(nmuts) / float(total)

subfafile = outfile.parent.joinpath(infile.stem + '.sub.fa')
# get fasta
shell.bedtools.getfasta('-name+', fi = ref, fo = subfafile, bed = mergedfile).fg

# index it
faIndex(subfafile)
faindex = subfafile.with_suffix('.fa.fai')
reader  = TsvReader(faindex, cnames = False)
regions = reader.dump(0)

params.m = mutrate
params._ = subfafile

para = Parallel(nthread, raiseExc = True)

def one_run(region):
	ps = params.copy()
	ps.r = region
	ps.o = outfile.with_suffix('.' + region + '.unmerged.vcf')
	shell.vcfy(**ps).fg

para.run(one_run, [(reg, ) for reg in regions])

tmpfile = outfile.parent.joinpath(outfile.stem + '.tmp.vcf')
shell.bcftools.concat(_ = [vcf for vcf in outfile.parent.glob('*.unmerged.vcf')],
	o = tmpfile, threads = nthread).fg

shell.rm_rf(*[str(vcf) for vcf in outfile.parent.glob('*.unmerged.vcf')])

# get the original configs and coordinates
faIndex(ref)

refindex = ref + '.fai'
contig_len = OrderedDiot()
for r in TsvReader(refindex, cnames = False):
	contig_len[r[0]] = r[1]

if samples is True:
	samples = ['Sample']
if samples and not isinstance(samples, list):
	samples = [samples]

unsorted = outfile.parent.joinpath(outfile.stem + '.unsorted.vcf')
with open(tmpfile) as fin, open(unsorted, 'w') as fout:
	for line in fin:
		if line.startswith('##contig=<') and contig_len:
			for contig, clen in contig_len.items():
				fout.write(f'##contig=<ID={contig},length={clen}>\n')
			contig_len = None
			if samples:
				fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Possible genotype">\n')
		elif line.startswith('##contig=<') and not contig_len:
			continue
		elif line.startswith('##'):
			fout.write(line)
		elif line.startswith('#'):
			if samples:
				fout.write(line.rstrip() + '\tFORMAT\t{}\n'.format('\t'.join(samples)))
			else:
				fout.write(line)
		else:
			parts = line.split('\t')
			coord = parts[0].split('::', 1)[1]
			chrom, startend = coord.split(':', 1)
			start = startend.split('-', 1)[0]
			parts[0] = chrom
			parts[1] = str(int(start) + int(parts[1]))
			parts[6] = 'PASS'
			if samples:
				parts[-1] = parts[-1].rstrip()
				parts.append('GT')
				for sample in samples:
					rd = random.choice([0,0,1,1,1,1,2,2,2,2])
					if rd == 0:
						parts.append('0|0')
					elif rd == 1:
						parts.append('0|1')
					else:
						parts.append('1|1')
				fout.write('\t'.join(parts) + '\n')
			else:
				fout.write('\t'.join(parts))

shell.bcftools.sort(_ = unsorted, o = outfile, T = tmpdir, O = 'z' if gz else 'v').fg
