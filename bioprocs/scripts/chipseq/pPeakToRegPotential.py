import math, gzip
peakfile = "{{peakfile}}"
genefile = "{{genefile}}"
arg_inst = {{args.signal | repr}}
arg_gf   = "{{args.genefmt}}"
arg_pf   = "{{args.peakfmt}}"
arg_wd   = int({{args.window | repr}})
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
			items   = line.split("\t")
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
			items   = line.split("\t")
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
			items   = line.split("\t")
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
			items   = line.split("\t")
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
		f.write ("%s\t%.3f\n" % (key, rp[key]))