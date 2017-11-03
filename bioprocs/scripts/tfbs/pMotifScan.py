import re
from os import path
from shutil import rmtree
{{ runcmd }}
{{ params2CmdArgs }}
{{ parallel }}

params   = {{args.params}}
ucsclinkTpl = {{args.ucsclink | quote}}

motifs = {}
with open({{in.tffile | quote}}) as f:
	for line in f:
		line = line.strip("\n")
		if not line: continue
		parts = line.split("\t")
		if len(parts) == 1:
			motifs[parts[0]] = parts[0]
		else:
			motifs[parts[0]] = parts[1]

{% if args.tool | lambda x: x == 'meme' %}

cmds   = []

params['thresh']    = {{args.pval | quote}}
params['verbosity'] = 4

for mid in motifs.keys():
	params['o']         = path.join({{out.outdir | quote}}, re.sub(r'[^\w_]', '', mid))
	params['motif']     = mid

	if path.isdir(params['o']): rmtree(params['o'])
	cmd = '{{args.meme}} %s {{args.tfmotifs | quote}} {{in.sfile | quote}}' % params2CmdArgs(params, dash='--', equal=' ')

	cmds.append(cmd)

parallel(cmds, {{args.nthread}})

with open({{out.outfile | quote}}, 'w') as fout:
	cont = ['#Chr', 'Start', 'End', 'Name', 'NormedScore', 'Strand', 'Motif', 'Sequence', 'StartOnSeq', 'StopOnSeq', 'Score', 'PValue', 'QValue', 'MatchedSeq', 'UCSCLink']
	fout.write('\t'.join(cont) + '\n')
	minscore = 9999
	maxscore = -1
	for mid in motifs.keys():
		with open(path.join({{out.outdir | quote}}, re.sub(r'[^\w_]', '', mid), 'fimo.txt')) as f:
			for line in f:
				line   = line.strip()
				if not line or line.startswith('#'): continue
				parts  = line.split('\t')
				score  = float(parts[5])
				minscore = min(minscore, score)
				maxscore = max(maxscore, score)

	for mid in motifs.keys():
		with open(path.join({{out.outdir | quote}}, re.sub(r'[^\w_]', '', mid), 'fimo.txt')) as f:
			for line in f:
				line   = line.strip()
				if not line or line.startswith('#'): continue
				parts  = line.split('\t')
				motif  = parts[0]
				# GENE::chr1:111-222 or ::chr1:111-222 or chr1:111-222
				seq    = parts[1] 
				strt0  = int(parts[2])
				end0   = int(parts[3])
				strand = parts[4]
				score  = int(1000 * (float(parts[5]) - minscore) / (maxscore - minscore))
				seqs   = re.split(r'[:-]', seq)
				chrom  = seqs[-3]
				start  = int(seqs[-2])
				end    = int(seqs[-1])
				name   = motifs[motif] + '::' + (seqs[0] if len(seqs)>3 else chrom + ':' + seqs[-2] + '-' + seqs[-1])
				ucsclink = ucsclinkTpl.format(chrom + ':' + str(start + strt0) + '-' + str(start + end0))

				cont   = [chrom, start + strt0, start + end0, name, score, strand, motif, seq, strt0, end0, parts[5], parts[6], parts[7], parts[8], ucsclink]
				cont   = [str(c) for c in cont]
				fout.write('\t'.join(cont) + '\n')

{% endif %}



