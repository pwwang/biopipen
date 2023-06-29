from pathlib import Path

import pandas
from biopipen.scripts.vcf.VcfFix_utils import HeaderContig, fix_vcffile
from biopipen.utils.reference import bam_index
from biopipen.utils.misc import run_command, dict_to_cli_args

bamfile = {{in.bamfile | quote}}  # pyright: ignore
snpfile = {{in.snpfile | repr}}  # pyright: ignore
outdir = Path({{out.outdir | quote}})  # pyright: ignore
cnvpytor = {{envs.cnvpytor | quote}}  # pyright: ignore
cnvnator2vcf = {{envs.cnvnator2vcf | quote}}  # pyright: ignore
samtools = {{envs.samtools | quote}}  # pyright: ignore
ncores = {{envs.ncores | int}}  # pyright: ignore
refdir = {{envs.refdir | quote}}  # pyright: ignore
genome = {{envs.genome | quote}}  # pyright: ignore
chrsize = {{envs.chrsize | quote}}  # pyright: ignore
filters = {{envs.filters | repr}}  # pyright: ignore
args = {{envs | repr}}  # pyright: ignore

del args['cnvpytor']
del args['ncores']
del args['cnvnator2vcf']
del args['samtools']
del args['refdir']
del args['genome']
del args['chrsize']
del args['filters']


bamfile = bam_index(bamfile, outdir, samtools, ncores)

NOSNP_COLS = [
    "CNVtype",
    "CNVregion",
    "CNVsize",
    "CNVlevel",
    "eval1",
    "eval2",
    "eval3",
    "eval4",
    "q0",
    "pN",
    "dG",
]

SNP_COLS = [
    "CNVtype",
    "CNVregion",
    "CNVsize",
    "CNVlevel",
    "eval1",
    "eval2",
    "eval3",
    "eval4",
    "q0",
    "pN",
    "dNS",
    "pP",
    "bin_size",
    "n",
    "delta_BA",
    "e",
    "baf_eval",
    "hets",
    "homs",
    "cn_1",
    "genotype",
    "likeliho",
    "cf_1",
    "cn_2",
    "genotype",
    "likeliho",
    "cf_2",
]

## Without snp
#0 CNV type: "deletion" or "duplication",
#1 CNV region (chr:start-end),
#2 CNV size,
#3 CNV level - read depth normalized to 1,
#4 e-val1 -- e-value (p-value multiplied by genome size divided by bin size) calculated using t-test statistics between RD statistics in the region and global,
#5 e-val2 -- e-value (p-value multiplied by genome size divided by bin size) from the probability of RD values within the region to be in the tails of a gaussian distribution of binned RD,
#6 e-val3 -- same as e-val1 but for the middle of CNV,
#7 e-val4 -- same as e-val2 but for the middle of CNV,
#8 q0 -- fraction of reads mapped with q0 quality in call region,
#9 pN -- fraction of reference genome gaps (Ns) in call region,
#10 dG -- distance from closest large (>100bp) gap in reference genome.

## With snp
#0 CNV type: "deletion", "duplication", or ”cnnloh",
#1 CNV region (chr:start-end),
#2 CNV size,
#3 CNV level - read depth normalized to 1,
#4 e-val1 -- e-value (p-value multiplied by genome size divided by bin size) calculated using t-test statistics between RD statistics in the region and global,
#5 e-val2 -- e-value (p-value multiplied by genome size divided by bin size) from the probability of RD values within the region to be in the tails of a gaussian distribution of binned RD,
#6 e-val3 -- same as e-val1 but for the middle of CNV,
#7 e-val4 -- same as e-val2 but for the middle of CNV,
#8 q0 -- fraction of reads mapped with q0 quality in call segments,
#9 pN -- fraction of reference genome gaps (Ns) within call region,
#10 dNS -- fraction of reference genome gaps (Ns) within call segments,
#11 pP -- fraction of P bases (1kGP strict mask) within call segments,
#12 bin_size – size of bins
#13 n – number of bins within call segments,
#14 delta_BAF – change in BAF from ½,
#15 e-val1 -- e-value RD based (repeted, reserved for future upgrades),
#16 baf_eval – e-value BAF based,
#17 hets – number of HETs,
#18 homs – number of HOMs,
#19 cn_1 – most likely model copy number,
#20 genotype_1 - most likely model genotype,
#21 likelihood_1 – most likely model likelihood,
#22 cf_1 -- most likely model cell fraction,
#23 cn_2 – the second most likely model copy number,
#24 genotype_2 - the second most likely model genotype,
#25 likelihood_2 – the second most likely model likelihood,
#26 cf_2 -- the second most likely model cell fraction.


def parse_cnvpytor_withoutsnp_line(line):
    parts = line.strip("\n").split("\t")
    cnvtype = parts[0]
    chrom = parts[1].split(":")[0]
    start, end = parts[1][len(chrom)+1:].split("-", 1)
    size = parts[2]
    level = parts[3]
    e1, e2, e3, e4 = parts[4:8]
    return dict(
        chrom=chrom,
        start=start,
        end=end,
        size=size,
        level=level,
        cnvtype=cnvtype,
        e1=e1,
        e2=e2,
        e3=e3,
        e4=e4,
        q0=parts[8],
        pn=parts[9],
        dg=parts[10],
    )

def parse_cnvpytor_withsnp_line(line):
    parts = line.strip("\n").split("\t")
    cnvtype = parts[0]
    chrom = parts[1].split(":")[0]
    start, end = parts[1][len(chrom)+1:].split("-", 1)
    size = parts[2]
    level = parts[3]
    e1, e2, e3, e4 = parts[4:8]

    return dict(
        chrom=chrom,
        start=start,
        end=end,
        size=size,
        level=level,
        cnvtype=cnvtype,
        e1=e1,
        e2=e2,
        e3=e3,
        e4=e4,
        q0 = parts[8],
        pn = parts[9],
        dns = parts[10],
        pp = parts[11],
        binsize = parts[12],
        n = parts[13],
        delta_baf = parts[14],
        e_val1 = parts[15],
        baf_eval = parts[16],
        hets = parts[17],
        homs = parts[18],
        cn1 = parts[19],
        genotype1 = parts[20],
        likelihood1 = parts[21],
        cf1 = parts[22],
        cn2 = parts[23],
        genotype2 = parts[24],
        likelihood2 = parts[25],
        cf2 = parts[26],
    )


def parsed_to_gff(i, parsed):
    outs = [
        parsed["chrom"],
        "CNVpytor",
        parsed["cnvtype"],
        parsed["start"],
        parsed["end"],
        parsed["level"],
        "+",
        ".",
        "; ".join(f"{name}={value}" for name, value in sorted(parsed.items()))
    ]
    return "\t".join(outs) + "\n"


def parsed_to_bed(i, parsed):
    outs = [
        parsed["chrom"],
        parsed["start"],
        parsed["end"],
        f"{parsed['cnvtype']}_{i}",
        parsed["level"],
        "+",
    ]
    return "\t".join(outs) + "\n"


def cnvpytor2other(infile, snp, out):
    outfile = Path(infile).with_suffix(f".{out}")
    with open(infile, "rt") as fin, open(outfile, "wt") as fout:
        for i, line in enumerate(fin):
            if snp:
                parsed = parse_cnvpytor_withsnp_line(line)
            else:
                parsed = parse_cnvpytor_withoutsnp_line(line)

            if out == "gff":
                outline = parsed_to_gff(i, parsed)
            else:
                outline = parsed_to_bed(i, parsed)

            fout.write(outline)


def load_chrsize():
    with open(chrsize) as f:
        for line in f:
            line = line.strip()
            if line:
                chrom, size = line.split()
                yield chrom, int(size)


def cnvpytor2vcf(infile, snp, fix=True):
    unfixedfile = Path(infile).with_suffix(f".unfixed.vcf")
    outfile = Path(infile).with_suffix(f".vcf")
    stdout = run_command(
        dict_to_cli_args(
            {
                "": cnvpytor2vcf,
                "reference": genome,
                "_": [infile, refdir],
            },
            prefix="-",
        ),
        stdout="return",
    )
    if fix:
        unfixedfile.write_text(stdout)

        fixes = [
            {
                "kind": "format",
                "id": "PE",
                "fix": lambda obj: setattr(obj, 'Type', 'String')
            },
            {
                "kind": "fields",
                "fix": lambda items: items.__setitem__(-1, Path(bamfile).stem)
            }
        ]
        for chrom, size in load_chrsize():
            fixes.append({
                "kind": "contig",
                "append": True,
                "fix": (
                    lambda obj, chrom=chrom, size=size:
                    HeaderContig(ID=chrom, length=size)
                )
            })

        fix_vcffile(unfixedfile, outfile, fixes)
    else:
        outfile.write_text(stdout)


def do_case():
    case = args.copy()

    chrom = case.pop("chrom", [])
    binsizes = case.pop("binsizes", [10000, 100000])
    snp = case.pop("snp", False)
    if snpfile is None:
        snp = False
    rootfile = outdir / "file.pytor"
    case["j"] = case.get("j", ncores)

    # read depth signal
    run_command(
        dict_to_cli_args(
            {
                "": cnvpytor,
                "root": rootfile,
                "rd": bamfile,
                "chrom": chrom,
                **case,
            },
            prefix="-",
        ),
        fg=True,
    )

    # predicting cnv
    run_command(
        dict_to_cli_args(
            {
                "": cnvpytor,
                "root": rootfile,
                "his": binsizes,
            },
            prefix="-",
        ),
        fg=True,
    )

    run_command(
        dict_to_cli_args(
            {
                "": cnvpytor,
                "root": rootfile,
                "partition": binsizes,
            },
            prefix="-",
        ),
        fg=True,
    )

    if snp is not False:
        snp["chrom"] = snp.get("chrom", chrom)
        snp["sample"] = False if not snp.get("sample", False) else snp["sample"]
        mask_snps = snp.pop("mask_snps", True)
        baf_nomask = snp.pop("baf_nomask", False)

        run_command(
            dict_to_cli_args(
                {
                    "": cnvpytor,
                    "root": rootfile,
                    "snp": snpfile,
                    **snp,
                },
                prefix="-",
            ),
            fg=True,
        )

        if mask_snps:
            run_command(
                dict_to_cli_args(
                    {
                        "": cnvpytor,
                        "root": rootfile,
                        "mask_snps": True,
                    },
                    prefix="-",
                ),
                fg=True,
            )

        run_command(
            dict_to_cli_args(
                {
                    "": cnvpytor,
                    "root": rootfile,
                    "baf": binsizes,
                    "nomask": baf_nomask,
                },
                prefix="-",
            ),
            fg=True,
        )

    # call
    for binsize in binsizes:
        outfile = outdir / f"calls{'.combined' if snp is not False else ''}.{binsize}.tsv"
        outfile_filtered = outdir / f"calls{'.combined' if snp is not False else ''}.{binsize}.filtered.tsv"
        run_command(
            dict_to_cli_args(
                {
                    "": cnvpytor,
                    "root": rootfile,
                    "call": "combined" if snp is not False else True,
                    "_": binsize,
                },
                prefix="-",
            ),
            stdout=outfile,
        )

        cnvpytor2other(outfile, bool(snp), "gff")
        cnvpytor2other(outfile, bool(snp), "bed")
        cnvpytor2vcf(outfile, bool(snp))

        # filter
        df = pandas.read_csv(outfile, sep="\t", header=None)
        df.columns = SNP_COLS if snp else NOSNP_COLS
        for key in (
            "CNVsize",
            "eval1",
            "eval2",
            "eval3",
            "eval4",
            "q0",
            "pN",
            "dG",
        ):
            # no `dG` column with snp data
            if snp and key == "dG":
                continue
            if key in filters and filters[key]:
                q1, q2 = filters[key]
                q1 = float(q1)
                q2 = float(q2)
                df = df[(df[key] >= q1) & (df[key] <= q2)]

        df.to_csv(outfile_filtered, sep="\t", index=False, header=False)
        cnvpytor2other(outfile_filtered, bool(snp), "gff")
        cnvpytor2other(outfile_filtered, bool(snp), "bed")
        cnvpytor2vcf(outfile_filtered, bool(snp))

        # plots
        manplot = outdir / f"manhattan.{binsize}.png"
        run_command(
            dict_to_cli_args(
                {
                    "": cnvpytor,
                    "root": rootfile,
                    "plot": ["manhattan", binsize],
                    "chrom": chrom,
                    "o": manplot,
                },
                prefix="-",
            ),
            fg=True,
        )

        circplot = outdir / f"circular.{binsize}.png"
        try:
            run_command(
                dict_to_cli_args(
                    {
                        "": cnvpytor,
                        "root": rootfile,
                        "plot": ["circular", binsize],
                        "chrom": chrom,
                        "o": circplot,
                    },
                    prefix="-",
                ),
                fg=True,
            )
        except RuntimeError as e:
            pass


do_case()
