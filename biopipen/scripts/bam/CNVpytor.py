from pathlib import Path

import warnings
import pandas
from datetime import datetime
from diot import Diot  # pyright: ignore
from biopipen.utils.reference import bam_index
from biopipen.utils.misc import run_command, dict_to_cli_args, logger

bamfile: str = {{in.bamfile | quote}}  # pyright: ignore # noqa
snpfile: str = {{in.snpfile | quote}}  # pyright: ignore
outdir = Path({{out.outdir | quote}})  # pyright: ignore
cnvpytor: str = {{envs.cnvpytor | quote}}  # pyright: ignore
samtools: str = {{envs.samtools | quote}}  # pyright: ignore
ncores: int = {{envs.ncores | int}}  # pyright: ignore
refdir = {{envs.refdir | quote}}  # pyright: ignore
genome = {{envs.genome | quote}}  # pyright: ignore
chrsize: str = {{envs.chrsize | quote}}  # pyright: ignore
filters: dict = {{envs.filters | repr}}  # pyright: ignore
args: Diot = {{envs | repr}}  # pyright: ignore

del args['cnvpytor']
del args['ncores']
del args['samtools']
del args['refdir']
del args['genome']
del args['chrsize']
del args['filters']


bamfile: Path = bam_index(bamfile, str(outdir), samtools, ncores=ncores)

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


def parse_chrom(chrom, chromdir):
    file = Path(chromdir) / f"{chrom}.fa"
    if not file.exists():
        warnings.warn(f"Chromosome file not found in refdir: {chrom}")
        return ""

    seq = ""
    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seq = ""
            else:
                seq += line
    return seq


def cnvpytor2vcf(infile, snp):
    # snp: in case to be used in the future
    outfile = Path(infile).with_suffix(f".vcf")
    # stdout = run_command(
    #     dict_to_cli_args(
    #         {
    #             "": cnvnator2vcf,
    #             "reference": genome,
    #             "_": [infile, refdir],
    #         },
    #         prefix="-",
    #     ),
    #     stdout="return",
    # )
    ## command hangs
    with open(infile) as fin, open(outfile, "w") as fout:
        fout.write("##fileformat=VCFv4.2\n")
        fout.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        fout.write(f"##reference={genome}\n")
        fout.write(f"##source=CNVpytor\n")
        for chrom, size in load_chrsize():
            fout.write(f"##contig=<ID={chrom},length={size}>\n")
        fout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
        fout.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
        fout.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
        fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        fout.write('##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">\n')
        fout.write('##INFO=<ID=natorP1,Number=1,Type=Float,Description="e-val by t-test">\n')
        fout.write('##INFO=<ID=natorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">\n')
        fout.write('##INFO=<ID=natorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">\n')
        fout.write('##INFO=<ID=natorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">\n')
        fout.write('##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">\n')
        fout.write('##INFO=<ID=natorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">\n')
        fout.write('##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">\n')
        fout.write('##ALT=<ID=DEL,Description="Deletion">\n')
        fout.write('##ALT=<ID=DUP,Description="Duplication">\n')
        fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fout.write('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n')
        fout.write('##FORMAT=<ID=PE,Number=1,Type=String,Description="Number of paired-ends that support the event">\n')
        fout.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{bamfile.stem}\n")
        prev_chrom, chrom_seq, count = "", "", 0
        for line in fin:
            # type, coor, length, rd, p1, p2, p3, p4, q0, pe = line.strip("\n").split()
            items = line.strip("\n").split()
            type, coor, length = items[:3]
            rd = float(items[3]) if len(items) > 3 else False
            p1 = items[4] if len(items) > 4 else ""
            p2 = items[5] if len(items) > 5 else ""
            p3 = items[6] if len(items) > 6 else ""
            p4 = items[7] if len(items) > 7 else ""
            q0 = items[8] if len(items) > 8 else ""
            pe = items[9] if len(items) > 9 else ""
            chrom, pos = coor.split(":")
            start, end = pos.split("-")
            start, end = int(start), int(end)
            is_del = type == "deletion"
            is_dup = type == "duplication"

            if not is_del and not is_dup:
                warnings.warn(f"Skipping unrecognized CNV type: {type}")
                continue

            if chrom != prev_chrom:
                chrom_seq = parse_chrom(chrom, refdir)
                prev_chrom = chrom

            count += 1
            info = f"END={end}"
            info += f";SVTYPE=DEL;SVLEN=-{length}" if is_del else f";SVTYPE=DUP;SVLEN={length}"
            info += ";IMPRECISE"
            info += f";natorRD={rd}" if rd is not False else ""
            info += f";natorP1={p1}" if p1 else ""
            info += f";natorP2={p2}" if p2 else ""
            info += f";natorP3={p3}" if p3 else ""
            info += f";natorP4={p4}" if p4 else ""
            info += f";natorQ0={q0}" if q0 else ""
            info += f";natorPE={pe}" if pe else ""

            gt = "GT"
            if rd is not False:
                gt += ":CN"
                gt += ":PE" if pe else ""
                gt += "\t"
                if is_del and rd < 0.25:
                    gt += "1/1:0"
                elif is_del and rd >= 0.25:
                    gt += "0/1:1"
                elif rd <= 1.75:
                    gt += "0/1:2"
                elif rd > 1.75 and rd <= 2.25:
                    gt += "1/1:2"
                elif rd > 2.25:
                    gt += f"./2:{rd:.0f}"
                else:
                    gt = "GT:PE\t./." if pe else "GT\t./."

                gt += f":{pe}" if pe else ""
            else:
                gt += "\t./."

            fout.write("\t".join(
                [
                    chrom,
                    str(start),
                    f"CNVpytor_{'del_' if is_del else 'dup_'}{count}",
                    chrom_seq[start - 1] if start < len(chrom_seq) else "N",
                    "<DEL>" if is_del else "<DUP>",
                    ".",
                    "PASS",
                    info,
                    gt,
                ]
            ) + "\n")


def do_case():
    case = args.copy()

    chrom = case.pop("chrom", [])
    binsizes = case.pop("binsizes", [10000, 100000])
    snp = case.pop("snp", False)
    if snpfile is None:
        snp = False
    rootfile = outdir / "file.pytor"
    case["j"] = case.get("j", ncores)

    logger.info("Reading depth signals ...")
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

    logger.info("Predicting CNVs ...")
    run_command(
        dict_to_cli_args(
            {
                "": cnvpytor,
                "root": rootfile,
                "his": binsizes,
            },
            prefix="-",
            dup_key=False,
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
            dup_key=False,
        ),
        fg=True,
    )

    if snp is not False:
        snp["chrom"] = snp.get("chrom", chrom)
        snp["sample"] = False if not snp.get("sample", False) else snp["sample"]
        mask_snps = snp.pop("mask_snps", True)
        baf_nomask = snp.pop("baf_nomask", False)

        logger.info("Importing SNP data ...")
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
            logger.info("Masking 1000 Genome SNPs ...")
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

        logger.info("Calculating BAF histograms ...")
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

    logger.info("Predicting CNV regions using joint caller ...")
    for binsize in binsizes:
        logger.info(f"- binsize: {binsize}")
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

        logger.info("  Converting to other formats ...")
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
        logger.info("  Plotting ...")
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
