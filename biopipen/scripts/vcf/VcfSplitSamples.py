import sys
import shlex
import concurrent.futures
from subprocess import Popen, check_output

infile: str = {{in.infile | quote}}  # pyright: ignore  # noqa
outdir: str = {{out.outdir | quote}}  # pyright: ignore
bcftools: str = {{envs.bcftools | repr}}  # pyright: ignore
gz = {{envs.gz | repr}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore
ncores: int = {{envs.ncores | int}}  # pyright: ignore
private = {{envs.private | repr}}  # pyright: ignore

if index:
    gz = True

# get sample list
o = check_output([bcftools, "query", "-l", infile])
samples = o.decode().strip().splitlines()


def do_sample(sample):
    """Do one sample"""
    print(f"- Processing sample {sample} ...")
    if gz:
        outfile = f"{outdir}/{sample}.vcf.gz"
    else:
        outfile = f"{outdir}/{sample}.vcf"

    # Get the index of the sample in the vcf file: sample in samples
    sample_idx = samples.index(sample)
    cmd = [
        bcftools,
        "view",
        "-e",
        f'GT[{sample_idx}]="RR" | GT[{sample_idx}]="mis"',
        "-Oz" if gz else "-Ov",
        "-s",
        sample,
        "-o",
        outfile,
        infile,
    ]
    print("  running:")
    print("  ", shlex.join(cmd))
    p = Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)
    p.wait()

    if index:
        cmd = [bcftools, "index", "-t", outfile]
        print("  running:")
        print("  ", shlex.join(cmd))
        Popen(cmd).wait()


with concurrent.futures.ProcessPoolExecutor(max_workers=ncores) as executor:
    executor.map(do_sample, samples)
