from datetime import date
from pathlib import Path

import cmdy
import numpy as np
from cyvcf2 import VCF, Writer
from pysam import FastaFile

inbed = {{in.inbed | quote}}  # pyright: ignore
outvcf = {{out.outvcf | quote}}  # pyright: ignore
tmpoutvcf = {{out.outvcf | append: ".tmp" | quote}}  # pyright: ignore
joboutdir = Path({{job.outdir | quote}})  # pyright: ignore
ref = {{envs.ref | quote}}  # pyright: ignore
headers = {{envs.headers | repr}}  # pyright: ignore
infos = {{envs.infos | repr}}  # pyright: ignore
base = {{envs.base | int}}  # pyright: ignore
formats = {{envs.formats | repr}}  # pyright: ignore
index = {{envs.index | repr}}  # pyright: ignore
bcftools = {{envs.bcftools | quote}}  # pyright: ignore
nonexisting_contigs = {{envs.nonexisting_contigs | quote}}  # pyright: ignore
genome = {{envs.genome | quote}}  # pyright: ignore
{{envs.helpers}}  # pyright: ignore
{% if envs.sample.startswith("lambda") %}  # pyright: ignore
instem = {{in.inbed | stem | quote}}  # pyright: ignore
sample = ({{envs.sample}})(instem)  # pyright: ignore
{% else %}  # pyright: ignore
sample = {{envs.sample | quote}}  # pyright: ignore
{% endif %}  # pyright: ignore
converters = {}
{% for key, val in envs.converters | items: %}  # pyright: ignore
converters[{{key | quote}}] = {{val}}  # pyright: ignore
{% endfor %}  # pyright: ignore

if ref == "":
    raise ValueError("Please specify the reference fasta file.")
fai = f"{ref}.fai"
if not Path(fai).exists():
    raise ValueError(f"{fai} does not exist.")

TEMPLATE_VCF = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate={date.today().strftime("%Y%m%d")}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0|0
"""

TEMPLATE_VCF_FILE = joboutdir / "template.vcf"
TEMPLATE_VCF_FILE.write_text(TEMPLATE_VCF)

vcf = VCF(TEMPLATE_VCF_FILE)

# Add source
vcf.add_to_header(f"##source=biopipen.ns.bed.Bed2Vcf")

# Add genome assembly
if genome:
    vcf.add_to_header(f"##reference={genome}")

vcf.add_info_to_header(
    {
        "ID": "END",
        "Number": "1",
        "Type": "Integer",
        "Description": "End position of the variant described in this record"
    }
)

vcf.add_format_to_header(
    {
        "ID": "GT",
        "Number": "1",
        "Type": "String",
        "Description": "Genotype",
    }
)

# Add contigs
contigs = set()
with open(fai) as f:
    for line in f:
        contig, length, *_ = line.strip().split("\t")
        contigs.add(contig)
        vcf.add_to_header(f"##contig=<ID={contig},length={length}>")

for header in headers:
    vcf.add_to_header(header)

for info in infos:
    vcf.add_info_to_header(info)

for fmt in formats:
    vcf.add_format_to_header(fmt)

header_types = {}
for header in vcf.header_iter():
    try:
        header_types[header["ID"]] = header["HeaderType"]
    except KeyError:
        pass

refseq = FastaFile(ref)
variant = next(vcf)
writer = Writer(tmpoutvcf, vcf)
try:
    with open(inbed) as f:
        for line in f:
            # chr,start,end,name,...
            items = line.rstrip("\n\r").split("\t")
            if nonexisting_contigs == "drop" and items[0] not in contigs:
                continue
            variant.CHROM = items[0]
            start = int(items[1])
            end = int(items[2])
            # The pos will be incremented by 1 with .set_pos()
            variant.set_pos(start - base)
            # If it is not a SNP
            if end - start > base:
                variant.INFO["END"] = end + 1 - base

            skip = False
            for key, converter in converters.items():
                val = converter(items)
                if val is None:
                    skip = True
                    continue
                if key == "ID":
                    variant.ID = val
                elif key == "REF":
                    variant.REF = val
                elif key == "ALT":
                    if isinstance(val, str):
                        variant.ALT = [val]
                    else:
                        variant.ALT = val
                elif key == "QUAL":
                    variant.QUAL = val
                elif key == "FILTER":
                    variant.FILTER = val

                elif header_types[key] == "FORMAT":
                    variant.set_format(key, np.array([val]))

                elif header_types[key] == "INFO":
                    variant.INFO[key] = val

            if skip:
                continue

            if "REF" not in converters:
                variant.REF = refseq.fetch(
                    variant.CHROM,
                    variant.POS - 1,
                    variant.POS
                )

            writer.write_record(variant)
finally:
    vcf.close()
    writer.close()

if index:
    cmdy.bcftools.sort(_=tmpoutvcf, O='z', o=outvcf, _exe=bcftools)
    cmdy.bcftools.index(_=outvcf, t=True, _exe=bcftools)

else:
    Path(tmpoutvcf).replace(outvcf)
