import os
import glob
import rtoml
import cmdy
from diot import Diot

bamfile = {{ in.bamfile | repr }}
snpfile = {{ in.snpfile | repr }}
outdir = {{ out.outdir | repr }}
freec = {{ envs.freec | repr }}
ncores = {{ envs.ncores | repr }}
bedtools = {{ envs.bedtools | repr }}
sambamba = {{ envs.sambamba | repr }}
samtools = {{ envs.samtools | repr }}
tabix = {{ envs.tabix | repr }}
rscript = {{ envs.rscript | repr }}
ref = {{ envs.ref | repr }}
refdir = {{ envs.refdir | repr }}
binsize = {{ envs.binsize | repr }}
args = {{ envs.args | repr }}

chrLenFile = f"{ref}.fai"
if snpfile:
    chrLenFile2 = f"{outdir}/chrLenFile.fai.txt"
    # Filter chrs in chrLenFile based on snps
    # get seqs from snpfile
    seqs = cmdy.tabix(snpfile, _exe=tabix, l=True).stdout.strip().splitlines()
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

config_ini = rtoml.dumps(config).replace('"', "")

with open(configfile, "w") as fconf:
    fconf.write(config_ini)

cmdy_args = {"_exe": freec, "_prefix": "-"}

cmd = cmdy.freec(conf=configfile, **cmdy_args).hold()
print()
print(cmd.strcmd, flush=True)

# freec not terminate when done
os.system(cmd.strcmd)

# plot cnvs
# get makeGraph.R
freec_path = os.path.realpath(cmdy.which(freec).stdout.strip())
mkgraph = os.path.join(os.path.dirname(freec_path), "makeGraph.R")
if not os.path.exists(mkgraph):
    raise RuntimeError("makeGraph.R not found")

try:
    ratiofile = glob.glob(f"{outdir}/FREEC-output/*_ratio.txt")[0]
except IndexError:
    raise RuntimeError("Control-FREEC failed to run") from None

rscript_path = os.path.realpath(cmdy.which(rscript).stdout.strip())
rpath = os.path.join(os.path.dirname(rscript_path), "R")

plotcmd = f"cat {mkgraph} | R --slave --args {config.general.ploidy} {ratiofile}"
print()
print(plotcmd)
os.system(plotcmd)
