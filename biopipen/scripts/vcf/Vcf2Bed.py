from cyvcf2 import VCF, Variant

infile: str = {{in.infile | quote}}  # pyright: ignore  # noqa: E999
outfile: str = {{out.outfile | quote}}  # pyright: ignore
# vcf, default 1
inbase = {{envs.inbase | int}}  # pyright: ignore
# bed, default 0
outbase = {{envs.outbase | int}}  # pyright: ignore


def variant2bedrecord(variant: Variant):
    start = variant.start + outbase - inbase
    end = variant.end or (start + 1 + outbase - inbase)
    return (
        variant.CHROM,
        str(start),
        str(end),
        variant.ID or f"{variant.CHROM}:{start}-{end}",
        str(variant.QUAL or 0),
        variant.INFO.get("Strand", "+"),
    )


def main():
    with open(outfile, "w") as fbed:
        for variant in VCF(infile):
            bedrecord = variant2bedrecord(variant)
            fbed.write('\t'.join(bedrecord) + '\n')


if __name__ == "__main__":
    main()
