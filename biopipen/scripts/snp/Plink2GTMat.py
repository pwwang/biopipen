
from os import path
from glob import glob
from biopipen.utils.misc import run_command, logger

indir: str = {{in.indir | quote}}  # noqa: E999 # pyright: ignore
outfile: str = {{out.outfile | quote}}  # pyright: ignore
plink: str = {{envs.plink | quote}}  # pyright: ignore
ncores: int = {{envs.ncores | repr}}  # pyright: ignore
transpose: bool = {{envs.transpose | repr}}  # pyright: ignore
samid: str = {{envs.samid | repr}}  # pyright: ignore
varid: str = {{envs.varid | repr}}  # pyright: ignore
trans_chr: dict = {{envs.trans_chr | repr}}  # pyright: ignore
missing_id: str = {{envs.missing_id | repr}}  # pyright: ignore
gtcoding: str = {{envs.gtcoding | repr}}  # pyright: ignore
trans_chr = trans_chr or {}

bedfile = glob(path.join(indir, '*.bed'))
if len(bedfile) == 0:
    raise FileNotFoundError(f"No .bed file found in `in.indir`")
elif len(bedfile) > 1:
    logger.warning(f"Multiple .bed files found in `in.indir`, using the first one.")

bedfile = bedfile[0]
input   = path.splitext(bedfile)[0]
output  = path.splitext(outfile)[0]

cmd = [
    plink,
    "--bfile", input,
    "--out", output,
    "--threads", ncores,
    "--keep-allele-order",
    "--recode", "A-transpose" if not transpose else "A",
]
# if transpose:
#     cmd += ["tabx"]

run_command(cmd, fg=True, env={"cwd": path.dirname(outfile)})


def _vcf_gtcoding(gt):
    try:
        return str(2 - int(gt))
    except (ValueError, TypeError):
        return "NA"


if not transpose:  # rows are variants, columns are samples
    # .traw file is created, tab-separated, with the following columns:
    trawfile = output + ".traw"
    # CHR     Chromosome code
    # SNP     Variant identifier
    # (C)M	  Position in morgans or centimorgans
    # POS     Base-pair coordinate
    # COUNTED Counted allele (defaults to A1), the actual alternative allele
    #           with --keep-allele-order
    # ALT     Other allele(s), comma-separated, the actual reference allele
    # <FID>_<IID>... Allelic dosages
    #   (0/1/2/'NA' for diploid variants, 0/2/'NA' for haploid)
    with open(trawfile, 'r') as fin:
        with open(outfile, 'w') as fout:
            samples = fin.readline().strip().split('\t')[6:]
            header = ["Variant"]
            for sam in samples:
                try:
                    fid, iid = sam.split('_')
                except ValueError:
                    raise ValueError(
                        f"Can't determine FID and IID from sample ID: {sam}, "
                        f"extra underscore (_) detected."
                    ) from None
                sam = samid.replace('{fid}', fid).replace('{iid}', iid)
                header.append(sam)
            fout.write('\t'.join(header) + '\n')

            for line in fin:
                line = line.strip().split('\t')
                chrom = trans_chr.get(line[0], line[0])
                var = line[1]
                if var == "." or var == "":
                    var = missing_id
                pos = line[3]
                ref = line[5]
                alt = line[4]
                variant = (
                    varid
                    .replace('{chr}', chrom)
                    .replace('{varid}', var)
                    .replace('{pos}', pos)
                    .replace('{ref}', ref)
                    .replace('{alt}', alt)
                )
                if gtcoding == "plink":
                    record = [variant] + line[6:]
                else:  # vcf
                    record = [variant] + [_vcf_gtcoding(x) for x in line[6:]]
                fout.write('\t'.join(record) + '\n')

else:
    # .raw file is created, tab-separated, with the following columns:
    rawfile = output + ".raw"
    # FID       Family ID
    # IID       Individual ID
    # PAT       Paternal ID
    # MAT       Maternal ID
    # SEX       Sex (1 = male, 2 = female, 0 = unknown)
    # PHENOTYPE Main phenotype value
    # <VariantID>... Allelic dosage (0/1/2/NA for diploid variants, 0/2/NA for haploid)
    #
    # Variant information may not be included in <VariantID>
    # We use the .bim file to get the variant information
    bimfile = input + ".bim"
    with open(rawfile, 'r') as fin:
        with open(outfile, 'w') as fout:
            header = ["Sample"]
            with open(bimfile, 'r') as fbim:
                for line in fbim:
                    line = line.strip().split('\t')
                    chrom = trans_chr.get(line[0], line[0])
                    var = line[1]
                    if var == "." or var == "":
                        var = missing_id
                    pos = line[3]
                    ref = line[5]
                    alt = line[4]
                    variant = (
                        varid
                        .replace('{chr}', chrom)
                        .replace('{varid}', var)
                        .replace('{pos}', pos)
                        .replace('{ref}', ref)
                        .replace('{alt}', alt)
                    )
                    header.append(variant)
            fout.write('\t'.join(header) + '\n')

            next(fin)  # skip header
            for line in fin:
                line = line.strip().split('\t')
                fid = line[0]
                iid = line[1]
                sam = samid.replace('{fid}', fid).replace('{iid}', iid)
                if gtcoding == "plink":
                    record = [sam] + line[6:]
                else:  # vcf
                    record = [sam] + [_vcf_gtcoding(x) for x in line[6:]]
                fout.write('\t'.join(record) + '\n')
