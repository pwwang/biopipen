"""Script for seq.pSeqMutate"""
# pylint: disable=invalid-name,undefined-variable,unused-import

from pathlib import Path
from bioprocs.utils.parallel import Parallel
from bioprocs.utils.reference import fa_index
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import shell2 as shell, logger

fafile = {{i.fafile | quote}}
mutfile = Path({{i.mutfile | quote}})
outfile = Path({{o.outfile | quote}})
refcheck = {{args.refcheck | repr}}
nthread = {{args.nthread | repr}}
seqkit = {{args.seqkit | repr}}

shell.load_config(seqkit=seqkit)
fa_index(fafile)

fai = fafile + '.fai'
indexes = {}

# write chroms to different threading files
threaddir = outfile.parent.joinpath('chroms')
shell.mkdir(p=threaddir)

def read_indexes():
    """read the indexes, for different chromsomes"""
    chreader = TsvReader(fai, cnames=False)
    for irow in chreader:
        #chr1:1000       31      130     31      32
        if ':' not in irow[0]:
            chrom = irow[0]
            start = 1
            end = int(irow[1])
        elif '-' not in irow[0]:
            chrom, start = irow[0].split(':')
            start = int(start)
            end = start + int(irow[1]) - 1
        else:
            chrom, startend = irow[0].split(':')
            start, end = startend.split('-')
            start, end = int(start), int(end)
            if end - start + 1 != int(irow[1]):
                logger.warning('Sequence (%s) length (%s) is different as name '
                               'defined: %s', irow[0], irow[1], end - start + 1)
            end = start + int(irow[1]) - 1
        indexes[irow[0]] = dict(chrom=chrom, start=start, end=end)
    chreader.close()

def split_mutfile():
    """Split mutation file by chromosome"""
    mutreader = TsvReader(mutfile, cnames=False)
    mfiles = {}
    for mut in mutreader:
        # mut:
        # chr start end name score strand ref alt
        found = False
        for seq, seqinfo in indexes.items():
            if seqinfo['chrom'] != mut[0]:
                continue
            if int(mut[2]) > seqinfo['end'] or int(mut[2]) < seqinfo['start']:
                continue

            found = True
            mfile = threaddir.joinpath(f'{mutfile.stem}.{seq}{mutfile.suffix}')
            if seq not in mfiles:
                mfiles[seq] = open(mfile, 'a')

            rec = [seq,
                   str(int(mut[2]) - seqinfo['start'] + 1),
                   str(int(mut[2]) - seqinfo['start'] + 1),
                   mut[3],
                   mut[4],
                   mut[5],
                   mut[6],
                   mut[7]]
            mfiles[seq].write('\t'.join(rec) + '\n')
        if not found:
            logger.warning('Mutation %s not found in sequence', list(mut))
    mutreader.close()

    for mfile in mfiles.values():
        mfile.close()

def check_ref(chrom):
    """Check the reference allele"""
    # get the reference allele
    mfile = threaddir.joinpath(f'{mutfile.stem}.{chrom}{mutfile.suffix}')
    # sequence for the mutation
    msfile = threaddir.joinpath(f'{mutfile.stem}.{chrom}.ref')
    shell.seqkit.subseq(
        bed=mfile,
        chr=chrom,
        _out=msfile,
        _=[fafile]
    )
    # msfile and mutfile should be 1-1
    mreader = TsvReader(mfile)
    ms = open(msfile)
    for mrow in mreader:
        mutname = ms.readline().lstrip('>').rstrip()
        mutref = ms.readline().strip()
        if not mutname.startswith(mrow[0] + '_'):
            raise ValueError('Unmatched sequence name for '
                             'reference allele check')
        if mutref != mrow[6]:
            raise ValueError('Reference allele differs: '
                             f'{mrow[6]} in mutation file, '
                             f'whereas {mutref} in sequence file.')
    ms.close()
    mreader.close()

def run_chrom(chrom):
    """Run each chromosome"""
    if refcheck:
        check_ref(chrom)

    mfile = threaddir.joinpath(f'{mutfile.stem}.{chrom}{mutfile.suffix}')
    refile = threaddir.joinpath(f'mutated.{chrom}.fa')
    if not mfile.exists():
        shell.seqkit.grep(p=chrom, _debug=True, _out=refile,
                          _=[fafile])
        return

    # get position and alternate allele for each mutation
    positions = []
    mreader = TsvReader(mfile, cnames=False)
    for mrow in mreader:
        positions.append(f'{mrow[2]}:{mrow[7]}')
    mreader.close()

    refile.write_text('')
    params = dict(
        _dupkey=True,
        p=positions,
        s=chrom,
        _out=refile,
        _debug=True
    )
    shell.seqkit.grep( # pylint: disable=expression-not-assigned
        _pipe=True, p=chrom, _=[fafile]
    ) | shell.seqkit.mutate(**params)

def concat():
    """Concatenate all chrom output to outfile, in the order of chroms"""
    for chrom in indexes:
        refile = threaddir.joinpath(f'mutated.{chrom}.fa')
        shell.cat(refile, _out_=outfile)

if __name__ == "__main__":
    read_indexes()
    split_mutfile()
    para = Parallel(nthread)
    para.run(run_chrom, [(chrom, ) for chrom in indexes])

    concat()
