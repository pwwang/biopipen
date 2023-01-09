import concurrent.futures
import cmdy

infile = {{in.infile | repr}}  # pyright: ignore
outdir = {{out.outdir | repr}}  # pyright: ignore
bcftools = {{envs.bcftools | repr}}  # pyright: ignore
gz = {{envs.gz | repr}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore
private = {{envs.private | repr}}  # pyright: ignore

if index:
    gz = True

# get sample list
cmd = cmdy.bcftools.query(l=True, _=infile, _exe=bcftools)
samples = cmd.stdout.splitlines()


def do_sample(sample):
    """Do one sample"""
    print(f"Processing sample {sample} ...")
    if gz:
        outfile = f"{outdir}/{sample}.vcf.gz"
    else:
        outfile = f"{outdir}/{sample}.vcf"

    cmd = cmdy.bcftools.view(
        _=infile,
        _exe=bcftools,
        s=sample,
        O="z" if gz else "v",
        o=outfile,
        private=private,
    ).hold()
    print("  running:")
    print("  ", cmd.strcmd)
    cmd.run(wait=True)

    if index:
        cmd = cmdy.bcftools.index(_=outfile, t=True, _exe=bcftools).hold()
        print("  running:")
        print("  ", cmd.strcmd)
        cmd.run(wait=True)


with concurrent.futures.ProcessPoolExecutor(max_workers=ncores) as executor:
    executor.map(do_sample, samples)
