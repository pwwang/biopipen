"""Script for seq.pSeqMutate"""
# pylint: disable=invalid-name,undefined-variable,unused-import,global-statement

import sys
from pathlib import Path
from bioprocs.utils.parallel import Parallel
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import logger, shell2 as shell

fafile = Path({{i.fafile | quote}})
mutfile = Path({{i.mutfile | quote}})
outfile = Path({{o.outfile | quote}})
nthread = {{args.nthread | repr}}
seqbase = {{args.seqbase | repr}}

# fasta index
indexes = {}
# # muts
nmuts = 0
muts_diffref = set()
# write chroms to different threading files
threaddir = outfile.parent.joinpath('seqs')
shell.mkdir(p=threaddir)

def split_fa():
    """Split fasta file into sequences"""
    if fafile.suffix == '.gz':
        import gzip
        ffa = gzip.open(fafile, 'rt')
    else:
        ffa = open(fafile, 'r')

    with ffa:
        seqfile = None
        lastseq = None
        for line in ffa:
            if line.startswith('>'):
                if seqfile:
                    seqfile.close()
                seqname0 = seqname = line[1:].split()[0].strip()
                shell.mkdir(p=threaddir / seqname)
                seqfile = open(threaddir.joinpath(seqname, 'seq.fa'), 'w')
                seqfile.write(line)

                indexes[seqname] = dict(file=threaddir.joinpath(
                    seqname, 'seq.fa'
                ))
                if '::' in seqname:
                    # adopt bedtools getfasta -name+
                    # DDX11L1::chr1:1-261869
                    seqname = seqname.split('::', 1)[1]
                if ':' not in seqname:
                    chrom = seqname
                    start = 1
                elif '-' not in seqname:
                    chrom, start = seqname.split(':')
                    start = int(start) + 1 - seqbase
                else:
                    chrom, startend = seqname.split(':')
                    start, _ = startend.split('-')
                    start = int(start) + 1 - seqbase
                indexes[seqname0]['chrom'] = chrom
                indexes[seqname0]['start'] = start
                indexes[seqname0]['seqlen'] = 0
                lastseq = seqname0
            else:
                indexes[lastseq]['seqlen'] += len(line.strip())
                seqfile.write(line)

def split_mutfile():
    """Split mutation file by chromosome"""
    mutfile_reseq = outfile.parent.joinpath(mutfile.name + '.reseq')
    mutreader = TsvReader(mutfile, cnames=False)
    mutwriter = TsvWriter(mutfile_reseq)
    global nmuts
    for mut in mutreader:
        # mut:
        # chr start end name score strand ref alt
        found = False
        for seq, seqinfo in indexes.items():
            if seqinfo['chrom'] != mut[0]:
                continue
            if (int(mut[2]) > seqinfo['start'] + seqinfo['seqlen'] - 1
                    or int(mut[2]) < seqinfo['start']):
                continue

            found = True
            #mfile = threaddir.joinpath(f'{mutfile.stem}.{seq}{mutfile.suffix}')
            #if seq not in mfiles:
            #    mfiles[seq] = open(mfile, 'a')

            rec = [seq,
                   str(int(mut[2]) - seqinfo['start'] + 1),
                   str(int(mut[2]) - seqinfo['start'] + 1),
                   mut[3],
                   mut[4],
                   mut[5],
                   mut[6],
                   mut[7]]
            mutwriter.write(rec)
        if not found:
            logger.warning('Mutation %s not found in sequence', list(mut))
        else:
            nmuts += 1
    mutreader.close()
    mutwriter.close()

    # split reseq
    mutfile_sorted = mutfile_reseq.with_suffix('.reseq.sorted')
    shell.sort(_=[mutfile_reseq], k=['1,1', '2,2n']).r > mutfile_sorted
    mutreader = TsvReader(mutfile_sorted, cnames=False)
    last_seq = None
    mutwriter = None
    for mut in mutreader:
        mfile = threaddir.joinpath(mut[0], 'mut.bed')
        if mut[0] != last_seq:
            last_seq = mut[0]
            mutwriter = TsvWriter(mfile)
        mutwriter.write(mut)
    mutwriter.close()
    mutreader.close()

def replace_pos(string, pos1, rep):
    """Replace a character by position in a string"""
    return string[:(pos1-1)] + rep + string[pos1:]

def run_chrom(chrom):
    """Run each chromosome"""
    mfile = threaddir.joinpath(chrom, 'mut.bed')
    seqfile = threaddir.joinpath(chrom, 'seq.fa')
    retfile = threaddir.joinpath(chrom, 'ret.fa')
    if not mfile.exists():
        shell.ln_s(seqfile, retfile, f=True)
        return

    mutreader = TsvReader(mfile, cnames=False)
    muts = mutreader.dump()
    mutreader.close()

    coocur = mutcur = 0
    global nmuts_diffref
    with open(seqfile, 'r') as fseq, open(retfile, 'w') as fret:
        fret.write(fseq.readline()) # seqname
        for line in fseq:
            line = line.strip()
            # see if any muts hit
            for i in range(mutcur, len(muts)):
                pos = int(muts[i][2])
                if coocur < pos <= coocur + len(line):
                    pos -= coocur
                    if line[pos-1].upper() != muts[i][6].upper():
                        muts_diffref.add(muts[i][3])
                        logger.warning('Mutation (%s) has a differnt reference '
                                       'allele (%s) than reported in the '
                                       'sequence (%s)', muts[i][3], muts[i][6],
                                       line[pos-1])
                    line = replace_pos(line, pos, muts[i][7])
                elif pos > coocur + len(line):
                    mutcur = i
                    coocur += len(line)
                    break
            fret.write(line + '\n')

def concat():
    """Concatenate all chrom output to outfile, in the order of chroms"""
    with open(outfile, 'w') as fout:
        for seqname in indexes:
            refile = threaddir.joinpath(seqname, 'ret.fa')
            fout.write(refile.read_text())

if __name__ == "__main__":
    split_fa()
    split_mutfile()

    para = Parallel(nthread, raiseExc=True)
    para.run(run_chrom, [(seqname, ) for seqname in indexes])
    logger.warning('%s/%s mutations have different reference alleles',
                   len(muts_diffref), nmuts)
    concat()
