import tempfile
from pathlib import Path

from ..core.config import config
from biopipen.utils.misc import run_command


def gztype(gzfile):
    import binascii

    with open(gzfile, "rb") as f:
        flag = binascii.hexlify(f.read(4))

    if flag == b"1f8b0804":
        return "bgzip"
    if flag == b"1f8b0808":
        return "gzip"
    return "flat"


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
        run_command(
            ["gunzip", "-c", infile], stdout=new_infile.with_suffix(""),
        )
        run_command(["bgzip", new_infile.with_suffix("")], fg=True)
    elif gt == "flat":
        run_command(["bgzip", "-c", infile], stdout=new_infile)
    else:
        # directory of infile may not have write permission
        new_infile.symlink_to(infile)

    run_command([tabix, "-p", preset, new_infile], fg=True)
    return new_infile


def _run_bam_index(
    bam,
    idxfile=None,
    tool="samtools",
    samtools=config.exe.samtools,
    sambamba=config.exe.sambamba,
    ncores=1,
):
    if tool == "samtools":
        cmd = [samtools, "index", "-@", ncores, bam, idxfile]
    else:
        cmd = [sambamba, "index", "-t", ncores, bam, idxfile]
    run_command(cmd, fg=True)


def bam_index(
    bam,
    bamdir=tempfile.gettempdir(),
    tool="samtools",
    samtools=config.exe.samtools,
    sambamba=config.exe.sambamba,
    ncores=1,
    ext=".bam.bai",
    force=False,
):
    """Index a bam file

    First look for the index file in the same directory as the bam file,
    if found, return the bam file. Otherwise, generate a symbolic link of the
    bam file in bamdir, and generate a index there, return the path to the
    symbolic link

    Args:
        bam: The path to the bam file
        bamdir: If index file can't be found in the directory as the bam file,
            create a symbolic link to the bam file, and generate the index
            here
        tool: The tool used to generate the index file, either `samtools` or
            `sambamba`
        samtools: The path to samtools
        sambamba: The path to sambamba
        ncores: Number of cores (threads) used to generate the index file
        ext: The ext of the index file, default `.bam.bai`, in case, `.bai` is
            also treated as index file
        force: Force to generate the index file, with given bamfile.
            Don't check if the index file exists.

    Returns:
        The bam file if index exists in the directory as the bam file.
        Otherwise symbolic link to the bam file in bamdir.
    """
    bam = Path(bam)
    indexfile = bam.with_suffix(ext)
    if force:
        _run_bam_index(
            bam,
            indexfile,
            tool,
            samtools,
            sambamba,
            ncores,
        )
        return bam

    if indexfile.is_file():
        return str(bam)

    linkfile = Path(bamdir).joinpath(bam.name)
    indexfile = linkfile.with_suffix(ext)

    if linkfile.exists() and not linkfile.samefile(bam):
        linkfile.unlink()
        if indexfile.exists():
            indexfile.unlink()

    if not linkfile.exists():
        linkfile.symlink_to(bam)

    if indexfile.is_file():
        return linkfile

    _run_bam_index(
        linkfile,
        indexfile,
        tool,
        samtools,
        sambamba,
        ncores,
    )
    return linkfile
