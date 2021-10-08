"""Script for bam.PatternCNV"""
# pylint: disable=undefined-variable,unused-import,bad-whitespace,invalid-name
# pylint: disable=not-a-mapping
from os.path import commonprefix
from pathlib import Path, PosixPath

from diot import Diot
from biopipen.utils import shell
from biopipen.utils.tsvio import TsvWriter

tumors = {{in.tumors | repr}}
normals = {{in.normals | repr}}
outdir = Path({{out.outdir | repr}})
patterncnv = {{args.patterncnv | repr}}
r = {{args.r | repr}}
samtools = {{args.samtools | repr}}
bedtools = {{args.bedtools | repr}}
targets = {{args.targets | repr}}
refflat = {{args.refflat | repr}}
gsize = {{args.gsize | repr}}
ref = {{args.ref | repr}}
params = {{args.params | repr}}

# check if tumor-normal samples are paired
if normals and len(tumors) != len(normals):
    raise ValueError("Tumor and normal samples are not paired.")

shell.load_config(patterncnv=patterncnv)
shell.mkdir(p=outdir)

samples = zip(tumors, normals or [None] * len(tumors))

# construct sample_info.txt file
sample_info_file = outdir / 'sample_info.txt'
sample_info_writer = TsvWriter(sample_info_file)
sample_info_writer.cnames = [
    'sample.name', 'subject.ID', 'sample.type', 'batch.ID', 'BAM.file'
]
sample_info_writer.write_head()

def get_bindir(bin):
    """Get the directory of the binary"""
    binpath = shell.which(bin) if Path(bin).parent == Path() else bin
    return Path(binpath).resolve().parent


subjects = set()
subject_index = 1
for tumor, normal in samples:
    tumor = Path(tumor)
    normal = normal and Path(normal)
    # get subjects
    if not normal:
        sample_info_writer.write([
            tumor.stem, tumor.stem, 'Germline', 1, tumor
        ])
    else:
        subject = commonprefix(tumor.stem, normal.stem)
        if not subject or subject in subjects:
            subject = f"Subject{subject_index}"
            subject_index += 1
        sample_info_writer.write([
            tumor.stem, subject, 'Somatic', 1, tumor
        ])
        sample_info_writer.write([
            normal.stem, subject, 'Germline', 1, normal
        ])
sample_info_writer.close()

# construct config.txt file
config_file = outdir / 'config.txt'
with config_file.open('w') as fconf:
    fconf.write(f"OUTPUT_DIR={outdir}\n")
    fconf.write(f"SAMPLE_INFO={sample_info_file}\n")
    fconf.write(f"CAPTUREKIT_BED={targets}\n")
    fconf.write(f"EXON_BED={refflat}\n")
    fconf.write(f"GENOME_SIZE={gsize}\n")
    fconf.write(f"GENOME_REF={ref}\n")

    fconf.write(f"PATTERNCNV={get_bindir(patterncnv)}\n")
    fconf.write(f"SAMTOOLS={samtools}\n")
    fconf.write(f"BEDTOOLS={get_bindir(bedtools)}\n")
    fconf.write(f"R={get_bindir(r)}\n")


shell.patterncnv(c=config_file, **params).fg()
