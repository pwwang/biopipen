import re
from os import path
from shutil import rmtree
{{ runcmd }}
{{ params2CmdArgs }}
{{ parallel }}

{% if args.cleanmname %}
def cleanMotifName(mname, fn = {{in.mfile | fn | quote}}):
	if 'HOCOMOCO' in fn:
		return mname.split('_HUMAN')[0]
	return mname
{% else %}
def cleanMotifName(mname, fn = ''):
	return mname
{% endif %}


params   = {{args.params}}
ucsclinkTpl = {{args.ucsclink | quote}}
{% if args.tool | lambda x: x == 'meme' %}
cmds   = []

params['thresh']    = {{args.pval | quote}}
params['verbosity'] = 4
mids = {{args.mids | lambda x: list(map(lambda x: x.strip(), x.split(','))) if not isinstance(x, list) else x}}
for mid in mids:
	params['o']         = path.join({{out.outdir | quote}}, re.sub(r'[^\w_]', '', mid))
	params['motif']     = mid

	if path.isdir(params['o']): rmtree(params['o'])
	cmd = '{{args.meme}} %s {{in.mfile | quote}} {{in.sfile | quote}}' % params2CmdArgs(params, dash='--', equal=' ')

	cmds.append(cmd)

parallel(cmds, {{args.nthread}})

with open({{out.outfile | quote}}, 'w') as fout:
	cont = ['#Chr', 'Start', 'End', 'Name', 'NormedScore', 'Strand', 'Motif', 'Sequence', 'StartOnSeq', 'StopOnSeq', 'Score', 'PValue', 'QValue', 'MatchedSeq', 'UCSCLink']
	fout.write('\t'.join(cont) + '\n')
	minscore = 9999
	maxscore = -1
	for mid in mids:
		with open(path.join({{out.outdir | quote}}, re.sub(r'[^\w_]', '', mid), 'fimo.txt')) as f:
			for line in f:
				line   = line.strip()
				if not line or line.startswith('#'): continue
				parts  = line.split('\t')
				score  = float(parts[5])
				minscore = min(minscore, score)
				maxscore = max(maxscore, score)

	for mid in mids:
		with open(path.join({{out.outdir | quote}}, re.sub(r'[^\w_]', '', mid), 'fimo.txt')) as f:
			for line in f:
				line   = line.strip()
				if not line or line.startswith('#'): continue
				parts  = line.split('\t')
				motif  = parts[0]
				seq    = parts[1]
				strt0  = int(parts[2])
				end0   = int(parts[3])
				strand = parts[4]
				score  = int(1000 * (float(parts[5]) - minscore) / (maxscore - minscore))
				seqs   = re.split(r'[:-]', seq)
				chrom  = seqs[2]
				start  = int(seqs[3])
				end    = int(seqs[4])
				name   = cleanMotifName(motif) + '::' + seqs[0]
				ucsclink = ucsclinkTpl.format(chrom + ':' + str(start + strt0) + '-' + str(start + end0))

				cont   = [chrom, start + strt0, start + end0, name, score, strand, motif, seq, strt0, end0, parts[5], parts[6], parts[7], parts[8], ucsclink]
				cont   = [str(c) for c in cont]
				fout.write('\t'.join(cont) + '\n')

{% endif %}



