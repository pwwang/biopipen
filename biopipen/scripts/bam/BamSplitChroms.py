from pathlib import Path
from biopipen.utils.misc import run_command
from biopipen.utils.reference import bam_index

bamfile: str = {{in.bamfile | quote}}  # pyright: ignore  # noqa
outdir: str = {{out.outdir | quote}}  # pyright: ignore
tool: str = {{envs.tool | quote}}  # pyright: ignore
samtools: str = {{envs.samtools | quote}}  # pyright: ignore
sambamba: str = {{envs.sambamba | quote}}  # pyright: ignore
ncores: int = {{envs.ncores | repr}}  # pyright: ignore
keep_other_sq = {{envs.keep_other_sq | repr}}  # pyright: ignore
chroms_to_keep = {{envs.chroms | repr}}  # pyright: ignore
should_index = {{envs.index | bool}}  # pyright: ignore


def _remove_other_sq(infile, chrom, outfile):
    exe = samtools if tool == "samtools" else sambamba
    print("\nRemoving other chromosomes in @SQ in header")
    header_cmd = [exe, "view", "-H", infile]
    header_p = run_command(  # type: ignore
        header_cmd,
        stdout=True,
        wait=False,
        print_command=True,
    )
    header = header_p.stdout.read().decode().strip().splitlines()  # type: ignore
    new_header = []
    for line in header:
        if line.startswith("@SQ"):
            if line.split()[1].split(":")[1] == chrom:
                new_header.append(line)
        else:
            new_header.append(line)

    newfile = Path(outdir) / f"{chrom}.sq_removed.sam"
    newfile.write_text("\n".join(new_header) + "\n")

    reheader_cmd = f"{exe} view {infile} >> {newfile}"
    run_command(reheader_cmd)

    # convert to bam
    dash_b = "-bS" if tool == "samtools" else "-S -f bam"
    convert_cmd = f"{exe} view {dash_b} -o {outfile} {newfile}"
    run_command(convert_cmd)

    # remove sam file
    newfile.unlink()


def use_samtools():
    print("\nUsing samtools")
    print("\nMaking sure bamfile is indexed")
    indexed_bamfile = bam_index(
        bamfile,
        tool="samtools",
        samtools=samtools,
        ncores=ncores,
    )
    if not chroms_to_keep:
        print("\nFetching chromosome names")
        cmd = (
            f"{samtools} view -H {indexed_bamfile} "
            "| grep '^@SQ' | cut -f 2 | cut -d ':' -f 2"
        )
        p = run_command(cmd, stdout=True, wait=False)
        chroms = p.stdout.read().decode().strip().splitlines()  # type: ignore
    else:
        print("\nUsing provided chromosomes")
        chroms = chroms_to_keep

    for chrom in chroms:
        print(f"Processing chromosome: {chrom}")
        outfile = (
            f"{outdir}/{chrom}.bam"
            if keep_other_sq
            else f"{outdir}/{chrom}.allsq.bam"
        )
        cmd = [
            samtools,
            "view",
            "-@",
            ncores,
            "-o",
            outfile,
            "-b",
            indexed_bamfile,
            chrom,
        ]
        run_command(cmd)

        if not keep_other_sq:
            _remove_other_sq(outfile, chrom, f"{outdir}/{chrom}.bam")
            Path(outfile).unlink()

        if should_index:
            print("\nIndexing bam file")
            bam_index(
                f"{outdir}/{chrom}.bam",
                tool="samtools",
                samtools=samtools,
                ncores=ncores,
                force=True,
            )

    print("\nDone")


def use_sambamba():
    print("\nUsing sambamba")
    print("\nMaking sure bamfile is indexed")
    indexed_bamfile = bam_index(
        bamfile,
        tool="sambamba",
        sambamba=sambamba,
        ncores=ncores,
    )
    if not chroms_to_keep:
        print("\nFetching chromosome names")
        cmd = (
            f"{sambamba} view -H {indexed_bamfile} "
            "| grep '^@SQ' | cut -f 2 | cut -d ':' -f 2"
        )
        p = run_command(cmd, stdout=True, wait=False)
        chroms = p.stdout.read().decode().splitlines()  # type: ignore
    else:
        print("\nUsing provided chromosomes")
        chroms = chroms_to_keep

    for chrom in chroms:
        print(f"Processing chromosome: {chrom}")
        outfile = (
            f"{outdir}/{chrom}.bam"
            if keep_other_sq
            else f"{outdir}/{chrom}.allSQ.bam"
        )
        cmd = [
            sambamba,
            "view",
            "-t",
            ncores,
            "-f",
            "bam",
            "-o",
            outfile,
            indexed_bamfile,
            chrom,
        ]
        run_command(cmd)

        if not keep_other_sq:
            _remove_other_sq(outfile, chrom, f"{outdir}/{chrom}.bam")
            Path(outfile).unlink()

        if should_index:
            print("\nIndexing bam file")
            bam_index(
                f"{outdir}/{chrom}.bam",
                tool="sambamba",
                sambamba=sambamba,
                ncores=ncores,
                force=True,
            )

    print("\nDone")


if __name__ == "__main__":
    if tool == "samtools":
        use_samtools()
    elif tool == "sambamba":
        use_sambamba()
    else:
        raise ValueError(f"Unknown tool: {tool}")
