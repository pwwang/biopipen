import os
import glob
import shutil
from diot import Diot  # type: ignore
from biopipen.utils.misc import dict_to_cli_args, run_command

bamfile = {{ in.bamfile | quote }}  # pyright: ignore # noqa
snpfile = {{ in.snpfile | quote }}  # pyright: ignore
outdir = {{ out.outdir | quote }}  # pyright: ignore
freec: str = {{ envs.freec | repr }}  # pyright: ignore
ncores = {{ envs.ncores | repr }}  # pyright: ignore
bedtools = {{ envs.bedtools | repr }}  # pyright: ignore
sambamba = {{ envs.sambamba | repr }}  # pyright: ignore
samtools = {{ envs.samtools | repr }}  # pyright: ignore
tabix = {{ envs.tabix | repr }}  # pyright: ignore
rscript: str = {{ envs.rscript | repr }}  # pyright: ignore
ref = {{ envs.ref | quote }}  # pyright: ignore
refdir = {{ envs.refdir | quote }}  # pyright: ignore
binsize = {{ envs.binsize | repr }}  # pyright: ignore
args = {{ envs.args | dict }}  # pyright: ignore

chrLenFile = f"{ref}.fai"
if snpfile:
    chrLenFile2 = f"{outdir}/chrLenFile.fai.txt"
    # Filter chrs in chrLenFile based on snps
    # get seqs from snpfile
    seqs = run_command(
        dict_to_cli_args(
            {
                "": [tabix, snpfile],
                "l": True,
            }
        ),
        stdout="return",
    ).strip().splitlines()  # type: ignore

    kept_seqs = []
    with open(chrLenFile, "r") as fin, open(chrLenFile2, "w") as fout:
        for line in fin:
            chrom = line.split()[0]
            if chrom in seqs:
                kept_seqs.append(chrom)
                fout.write(line)
    if not kept_seqs:
        raise ValueError(
            "No contifs left in chrLenFile(.fai) after "
            "filtering with snp data, maybe wrong snp data? "
            "Also make sure the chromosome names match."
        )
    chrLenFile = chrLenFile2

configfile = f"{outdir}/config.ini"

config = Diot(args)
# See http://boevalab.inf.ethz.ch/FREEC/tutorial.html#CONFIG
# for full configuration
config.general |= Diot(
    bedtools=bedtools,
    chrFiles=refdir,
    chrLenFile=chrLenFile,
    outputDir=f"{outdir}/FREEC-output",
    sambamba=sambamba,
    maxThreads=ncores,
    samtools=samtools,
)
if "window" not in config.general:
    config.general.window = binsize

config.sample |= Diot(
    mateFile=bamfile,
    inputFormat="BAM",
)

config.BAF |= Diot(
    fastaFile=ref,
    makePileup=snpfile,
)

os.makedirs(f"{outdir}/FREEC-output", exist_ok=True)

config_ini = config.to_toml().replace('"', "")   # type: ignore

with open(configfile, "w") as fconf:
    fconf.write(config_ini)

# Does it terminate when freec is done?
run_command(
    dict_to_cli_args({"": freec, "conf": configfile}, prefix="-"),
    fg=True,
)

# plot cnvs
# get makeGraph.R
freec_path = os.path.realpath(shutil.which(freec).strip())  # type: ignore
mkgraph = os.path.join(os.path.dirname(freec_path), "makeGraph.R")
if not os.path.exists(mkgraph):
    raise RuntimeError("makeGraph.R not found")

try:
    ratiofile = glob.glob(f"{outdir}/FREEC-output/*_ratio.txt")[0]
except IndexError:
    raise RuntimeError("Control-FREEC failed to run") from None

rscript_path = os.path.realpath(shutil.which(rscript).strip())  # type: ignore
rpath = os.path.join(os.path.dirname(rscript_path), "R")

plotcmd = f"cat {mkgraph} | R --slave --args {config.general.ploidy} {ratiofile}"
print()
print(plotcmd)
os.system(plotcmd)
