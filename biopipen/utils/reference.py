import tempfile
from pathlib import Path

import cmdy
from ..core.config import config

def gztype(gzfile):
    import binascii
    with open(gzfile, 'rb') as f:
        flag = binascii.hexlify(f.read(4))

    if flag == b'1f8b0804':
        return 'bgzip'
    if flag == b'1f8b0808':
        return 'gzip'
    return 'flat'

def tabix_index(infile, preset, tmpdir=None, tabix=config.exe.tabix):
    """Index input file using tabix

    1. Try to check if there is an index file in the same directory where infile
       is.
    2. If so, return the infile
    3. Otherwise, check if infile is bgzipped, if not bgzip it and save
       it in tmpdir
    4. Index the bgzipped file and return the bgzipped file

    Returns:
        The infile itself or re-bgzipped infile. This file comes with the
        index file in the same directory
    """
    if tmpdir is None:
        tmpdir = Path(tempfile.mkdtemp(prefix="biopipen_tabix_index_"))
    else:
        tmpdir = Path(tmpdir)

    infile = Path(infile)
    gt = gztype(infile)

    if gt == "bgzip" and infile.with_suffix(infile.suffix + ".tbi").is_file():
        # only bgzipped file is possible to have index file
        return infile

    # /path/to/some.vcf -> some.vcf
    # /path/to/some.vcf.gz -> some.vcf
    basename = infile.stem if infile.name.endswith(".gz") else infile.name

    # try bgzip infile
    new_infile = tmpdir / (basename + ".gz")
    if gt == "gzip":
        # re-bgzip
        cmdy.gunzip(infile, c=True).r() > new_infile.with_suffix("")
        cmdy.bgzip(new_infile.with_suffix(""))
    elif gt == "flat":
        cmdy.bgzip(infile, c=True).r() > new_infile
    else:
        # directory of infile may not have write permission
        new_infile.symlink_to(infile)

    cmdy.tabix(p=preset, _=new_infile, _exe=tabix)
    return new_infile
