"""Reference file and index file managing"""
from os import path, readlink
from bioprocs.utils import shell2 as shell, gztype


def check(ref):
    """Check if ref exists"""
    if not ref or not path.exists(ref):
        raise Exception('Reference file not exists: %s' % ref)


def check_index(refindex):
    """Check if index exists"""
    if not refindex:
        refindex = [refindex]
    return all([path.exists(ri) for ri in refindex])

checkIndex = check_index # pylint: disable=invalid-name

def build_index(ref, cmd, ref2=None, cmd2=None):
    """Build index file"""
    try:
        shell.bash(c=cmd)
        return ref
    except Exception:
        try:
            if not path.exists(ref2):
                shell.ln_s(ref, ref2)
            shell.bash(c=cmd2)
            return ref2
        except Exception:
            return None
buildIndex = build_index # pylint: disable=invalid-name

def fa_index(fa, ext='.fai', samtools='samtools'): # pylint: disable=invalid-name
    """Index fasta file"""
    if not ext.startswith('.'):
        ext = '.' + ext

    shell.load_config(samtools=samtools)

    expected_index = path.join(path.dirname(fa), path.basename(fa) + ext)
    if path.isfile(expected_index):
        return
    # if fa is not a link, there is nowhere else to find index,
    # create it using samtools
    if not path.islink(fa):
        if samtools:
            shell.samtools.faidx(fa)
        else:
            raise ValueError('Index not found: {}'.format(expected_index))
        return
    # find the index in original directory
    origfa = readlink(fa)
    orig_index = origfa + ext
    if path.isfile(orig_index):
        shell.ln_s(orig_index, expected_index)
        return
    # find the index in realpath directory
    realfa = path.realpath(fa)
    real_index = path.splitext(realfa)[0] + ext
    if path.isfile(real_index):
        shell.ln_s(real_index, expected_index)
        return
    # if all failed, create it
    if samtools:
        shell.samtools.faidx(fa)
    else:
        raise ValueError('Index not found: {}'.format(fa))

faIndex = fa_index # pylint: disable=invalid-name

def bam_index(bam, ext='.bam.bai', samtools='samtools', nthread=1):
    """
    Index bam files
    If bam file is a link, try to find the index file in its orginal
    directory or its realpath directory
    If nothing found, try to create the index file using samtools
    @params:
        `ext`: The expected extension of index file. Default: `.bam.bai`
            - Some tools requird `XXX.bai` without `.bam`
        `samtools`: The path to samtools. Default: `samtools`
            - If it's None, then an exception will raised
              instead of creating the index file
        `nthread`: The # threads used to create the index file. Default: `1`
    """
    if not ext.startswith('.'):
        ext = '.' + ext
    # /path/to/some.bam -> some.bam
    bname = path.basename(bam)
    # /path/to/some.bam -> /path/to/
    dname = path.dirname(bam)
    # some.bam -> some
    fname = path.splitext(bname)[0]
    # some -> some
    # [1]some -> some
    #rname = fname.split(']', 1)[1] if fname.startswith('[') else fname

    shell.load_config(samtools=samtools)
    # /path/to/some.bam.bai
    expected_index = path.join(dname, fname + ext)
    if path.isfile(expected_index):
        return
    # if bam is not a link, there is nowhere else to find index,
    # create it using samtools
    if not path.islink(bam):
        if samtools:
            shell.samtools.index({'@': nthread}, b=True, _=[bam, expected_index])
        else:
            raise ValueError('Index not found: {}'.format(bam))
        return
    # find the index in original directory
    origbam = readlink(bam)
    orig_index = path.splitext(origbam)[0] + ext
    if path.isfile(orig_index):
        shell.ln_s(orig_index, expected_index)
        return
    # find the index in realpath directory
    realbam = path.realpath(bam)
    real_index = path.splitext(realbam)[0] + ext
    if path.isfile(real_index):
        shell.ln_s(real_index, expected_index)
        return
    # if all failed, create it
    if samtools:
        shell.samtools.index({'@': nthread}, b=True, _=[bam, expected_index])
    else:
        raise ValueError('Index not found: {}'.format(bam))

bamIndex = bam_index # pylint: disable=invalid-name

def tabix_index(filename, type, tabix='tabix'): # pylint: disable=redefined-builtin
    """Use tabix to index file"""
    # /path/to/some.vcf -> some.vcf
    # /path/to/some.vcf.gz -> some.vcf
    filename = str(filename)
    bname = path.basename(
        filename[:-3]) if filename.endswith('.gz') else path.basename(filename)
    # /path/to/some.bam -> /path/to/
    dname = path.dirname(filename)
    # some.vcf -> some
    # some.vcf.gz -> some
    fname = path.splitext(bname[:-3] if bname.endswith('.gz') else bname)[0]
    # some -> some
    # [1]some -> some
    #rname = fname.split(']', 1)[1] if fname.startswith('[') else fname

    expected_index = path.join(dname, fname + '.%s.gz.tbi' % type)
    if path.isfile(expected_index):
        return filename if filename.endswith('.gz') else filename + '.gz'

    shell.load_config(tabix=tabix)
    # type could bed3, bed6..
    ptype = 'bed' if type.startswith('bed') else type
    gt = gztype(filename) # pylint: disable=invalid-name
    if gt == 'bgzip':
        if path.islink(filename):
            linkfile = readlink(filename)
            if path.isfile(linkfile + '.tbi'):
                shell.ln_s(linkfile + '.tbi', expected_index)
                return filename
            realfile = path.realpath(filename)
            if path.isfile(realfile + '.tbi'):
                shell.ln_s(realfile + '.tbi', expected_index)
                return realfile
        shell.tabix(p=ptype, _=filename).fg
        return filename
    if gt == 'gzip':
        bgzfile = path.join(dname, bname + '.bgz.' + type)
        shell.gunzip(filename, c=True).r > bgzfile
        shell.bgzip(bgzfile)
        shell.tabix(p=ptype, _=bgzfile).fg
        return bgzfile + '.gz'
    shell.bgzip(filename, c=True).r > filename + '.gz'
    shell.tabix(p=ptype, _=filename + '.gz').fg
    return filename + '.gz'

tabixIndex = tabix_index # pylint: disable=invalid-name

def vcf_index(vcf, tabix='tabix'):
    """Index vcf file"""
    return tabix_index(vcf, 'vcf', tabix)

vcfIndex = vcf_index # pylint: disable=invalid-name

def bed_index(bed, tabix='tabix'):
    """Index bed file"""
    return tabix_index(bed, 'bed', tabix)

bedIndex = bed_index # pylint: disable=invalid-name
