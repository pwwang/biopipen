from pyppl import proc

"""
@name:
	pPeakToRegPotential
@description:
	Convert peaks to regulatory potential score for each gene
	The formula is:
	```
		                 -(0.5 + 4*di/d0)
		PC = sum (pi * e                  )
	```
	Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489297/
@input:
	`peakfile:file`: The BED/peak file for peaks
	`genefile:file`: The BED file for gene coordinates
@output:
	`outfile:file`: The regulatory potential file for each gene
@args:
	`intensity`: `pi` in the formula. Boolean value, whether use the peak intensity or not, default: `True`,
	`geneformat`: The format for `genefile`, default: `ucsc+gz`. It could be:
		- ucsc or ucsc+gz: typically, you can download from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
		- bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 4th column required as gene identity.
	`peakformat`: The format for `peakfile`, default: `peak`. It could be:
		- peak or peak+gz: (either [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) or [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13), the 7th column will be used as intensity
		- bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 5th column will be used as intensity.
	`window`: `2 * d0` in the formula. The window where the peaks fall in will be consided, default: `100000`. 
```
		|--------- window ----------|
		|---- d0 -----|
		|--- 50K --- TSS --- 50K ---|
		     ^ (peak center)
		     |-- di --|
```
"""
pPeakToRegPotential = proc ()
pPeakToRegPotential.input     = "peakfile:file, genefile:file"
pPeakToRegPotential.output    = "outfile:file:{{peakfile.fn}}.rp"
pPeakToRegPotential.args      = {'intensity': True, 'geneformat': 'ucsc+gz', 'peakformat': 'peak', 'window': 100000}
pPeakToRegPotential.defaultSh = "python"
pPeakToRegPotential.script    = """
import math, gzip
peakfile = "{{peakfile}}"
genefile = "{{genefile}}"
arg_inst = {{proc.args.intensity}} == True
arg_gf   = "{{proc.args.geneformat}}"
arg_pf   = "{{proc.args.peakformat}}"
arg_wd   = int({{proc.args.window}})
d0       = arg_wd / 2

assert (isinstance(arg_inst, bool))
assert (arg_gf in ['ucsc', 'bed', 'ucsc+gz', 'bed+gz'])
assert (arg_pf in ['peak', 'bed', 'peak+gz', 'bed+gz'])

open_gf = open_pf = open
if arg_gf.endswith ('+gz'):
	arg_gf  = arg_gf[:-3]
	open_gf = gzip.open
if arg_pf.endswith ('+gz'):
	arg_pf  = arg_pf[:-3]
	open_pf = gzip.open

# read genes
genes = {}
if arg_gf == 'bed':
	with open_gf (genefile) as f:
		for line in f:
			line    = line.strip()
			if not line or line.startswith('track') or line.startswith('#'): continue
			items   = line.split("\\t")
			chr     = items[0]
			start   = int(items[1])
			end     = int(items[2])
			gene    = items[3]
			strand  = '-' if len(items)>5 and items[5] == '-' else '+'
			tss     = start if strand == '+' else end
			rstart  = tss - d0
			rend    = tss + d0
			genes[gene] = [chr, start, end, tss, rstart, rend]
else:
	with open_gf (genefile) as f:
		for line in f:
			line    = line.strip()
			if not line or line.startswith('track') or line.startswith('#'): continue
			items   = line.split("\\t")
			chr     = items[2]
			start   = int(items[4])
			end     = int(items[5])
			gene    = items[12]
			strand  = items[3]
			tss     = start if strand == '+' else end
			rstart  = tss - d0
			rend    = tss + d0
			genes[gene] = [chr, start, end, tss, rstart, rend]

# read peaks
peaks = {}
if arg_pf == 'peak':
	with open_pf (peakfile) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith('track') or line.startswith('#'): coninue
			items   = line.split("\\t")
			chr     = items[0]
			start   = int(items[1])
			end     = int(items[2])
			signal  = float(items[6])
			if peaks.has_key(chr):
				peaks[chr].append ([start, end, (start+end) / 2, signal])
			else:
				peaks[chr] = [[start, end, (start+end) / 2, signal]]
else:
	with open_pf (peakfile) as f:
		for line in f:
			line = line.strip()
			if not line or line.startswith('track') or line.startswith('#'): coninue
			items   = line.split("\\t")
			chr     = items[0]
			start   = int(items[1])
			end     = int(items[2])
			signal  = float(items[4])
			if peaks.has_key(chr):
				peaks[chr].append ([start, end, (start+end) / 2, signal])
			else:
				peaks[chr] = [[start, end, (start+end) / 2, signal]]

for key, val in peaks.iteritems():
	peaks[key] = sorted (val, cmp = lambda x, y: x[0] - y[0])

rp = {}
for gene, ginfo in genes.iteritems():
	(gchr, gstart, gend, gtss, grstart, grend) = ginfo
	rp[gene] = 0
	if not peaks.has_key(gchr): continue
	for pinfo in peaks[gchr]:
		(pstart, pend, pcenter, psignal) = pinfo
		if pcenter < grstart: continue
		if pcenter > grend: break
		score  = psignal if arg_inst else 1
		score *= math.exp (-(.5 + 4*abs(pcenter - tss)/d0))
		rp[gene] += score

with open ("{{outfile}}", 'w') as f:
	for key in sorted (rp, key=rp.get, reverse = True):
		f.write ("%s\\t%.3f\\n" % (key, rp[key]))
"""




