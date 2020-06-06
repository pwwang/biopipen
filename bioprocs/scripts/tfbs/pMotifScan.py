"""Script for tfbs.pMotifScan"""
# pylint: disable=undefined-variable,unused-import,invalid-name

import re
from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell, logger
from bioprocs.utils.meme import MemeReader
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord
from bioprocs.utils.parallel import Parallel

tffile = {{i.tffile | quote}}
sfile = {{i.sfile | quote}}
outfile = {{o.outfile | quote}}
outdir = {{o.outdir | quote}}
tool = {{args.tool | quote}}
meme = {{args.meme | quote}}
fimo = {{args.fimo | quote}}
params = {{args.params | repr}}
tfmotifs = {{args.tfmotifs | quote}}
cutoff = {{args.cutoff | ?isinstance: dict | !: dict(by='p', value=_) | $repr}}
ucsclink = {{args.ucsclink | repr}}
nthread = {{args.nthread | repr}}

shell.load_config(fimo=fimo or meme)

# get all motifs
logger.info('Loading motif names ...')
reader = TsvReader(tffile, cnames=False)
ncol = len(next(reader))
reader.rewind()
if ncol == 1:
    motifns = {r[0]: r[0] for r in reader}
else:
    motifns = {r[0]: r[1] for r in reader}
logger.info('%s motif names read.', len(motifns))

# match motifs
logger.info('Matching motifs in database ...')
mnames = [m.name for m in MemeReader(tfmotifs)]
motifs = {k: v for k, v in motifns.items() if k in mnames}
logger.info('%s motifs loaded', len(motifs))

if tool == 'meme':
    cmdparams = []
    if 'p' in cutoff:
        params.thresh = cutoff['p']
    else:
        params['qv-thresh'] = cutoff['q']
    params.verbosity = 4
    params._ = [tfmotifs, sfile]
    for motif, name in motifs.items():
        prms = params.copy()
        prms.oc = path.join(outdir, name + '.' +
                            re.sub(r'[^\w_]', '', motif))
        prms.motif = motif
        cmdparams.append((prms,))
    Parallel(nthread, raiseExc=True).run(shell.fimo, cmdparams)

    writer = TsvWriter(outfile)
    writer.cnames = [
        "CHR", "START", "END", "NAME", "SCORE", "STRAND",
        "MOTIF", "SEQ", "STARTONSEQ", "STOPONSEQ", "RAWSCORE", "PVAL",
        "QVAL", "MATCHEDSEQ"
    ]
    if ucsclink:
        writer.cnames.append("UCSCLINK")
    writer.writeHead(callback=lambda cnames: "#" + "\t".join(cnames))

    def rowfactory(r):
        """Row factory to refactorize the output"""
        if 'p-value' not in r:
            return None
        r.PVAL = float(r['p-value'])
        r.QVAL = float(r['q-value'])
        if 'p' in cutoff and r.PVAL >= cutoff['p']:
            return None
        elif 'q' in cutoff and r.QVAL >= cutoff['q']:
            return None
        r.RAWSCORE = r.score
        try:
            r.SCORE = int(float(r.score) * 10)
        except TypeError:
            r.SCORE = 0
        r.STRAND = r.strand
        r.MOTIF = r.motif_id
        # split motif_alt_id
        # GENE or GENE::chr1:111-222 or ::chr1:111-222 or chr1:111-222
        r.SEQ = r.sequence_name
        if '::' in r.SEQ:
            seqname, seqcoord = r.SEQ.split('::', 1)
            r.CHR, start, end = re.split(r'[:-]', seqcoord)
            start = int(start)
            end = int(end)
        elif ':' in r.SEQ:
            seqname, seqcoord = r.SEQ, r.SEQ
            r.CHR, start, end = re.split(r'[:-]', seqcoord)
            start = int(start)
            end = int(end)
        else:
            seqname, seqcoord = r.SEQ, ''
            r.CHR, start, end = '-', 0, 0
        r.NAME = seqname
        r.STARTONSEQ = r.start
        r.STOPONSEQ = r.stop
        r.START = start + int(r.start)
        r.END = start + int(r.stop)
        #r.QVAL = r['q-value']
        r.MATCHEDSEQ = r.matched_sequence
        if ucsclink:
            r.UCSCLINK = ucsclink.format('{chrom}:{start}:{end}'.format(
                chrom=r.CHR,
                start=start,
                end=end
            ))
        return r

    for motif, name in motifs.items():
        retfile = path.join(outdir, name + '.' +
                            re.sub(r'[^\w_]', '', motif), 'fimo.tsv')
        # motif_id	motif_alt_id	sequence_name	start	stop	strand
        # score	p-value	q-value	matched_sequence
        logger.info('- Merging %s ...', retfile)
        reader = TsvReader(
            retfile,
            #cnames  = lambda header: header[1:].strip().split('\t'),
            cnames=True,
            row=rowfactory,
            comment=None)
        for row in reader:
            if not row:
                continue
            row.NAME = name + '::' + row.NAME
            writer.write(row)
    writer.close()
else:
    raise ValueError('Unsupported tool: {}'.format(tool))
