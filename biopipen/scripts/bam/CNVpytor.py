from pathlib import Path

import cmdy
from biopipen.scripts.vcf.VcfFix_utils import HeaderContig, fix_vcffile

bamfile = {{in.bamfile | quote}}  # pyright: ignore
snpfile = {{in.snpfile | repr}}  # pyright: ignore
outdir = Path({{out.outdir | quote}})  # pyright: ignore
cnvpytor = {{envs.cnvpytor | quote}}  # pyright: ignore
cnvnator2vcf = {{envs.cnvnator2vcf | quote}}  # pyright: ignore
ncores = {{envs.ncores | int}}  # pyright: ignore
refdir = {{envs.refdir | quote}}  # pyright: ignore
genome = {{envs.genome | quote}}  # pyright: ignore
chrsize = {{envs.chrsize | quote}}  # pyright: ignore
args = {{envs | repr}}  # pyright: ignore

del args['cnvpytor']
del args['ncores']
del args['cnvnator2vcf']
del args['refdir']
del args['genome']
del args['chrsize']

cmdy_args = {"_exe": cnvpytor, "_prefix": "-", "_deform": None}

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


def cnvpytor2vcf(infile, snp):
    unfixedfile = Path(infile).with_suffix(f".unfixed.vcf")
    outfile = Path(infile).with_suffix(f".vcf")
    cnvnator2vcf_cmd = cmdy.cnvnator2vcf(
        reference=genome,
        _=[infile, refdir],
        _exe=cnvnator2vcf,
        _prefix="-",
    ).h()
    print()
    print(cnvnator2vcf_cmd.strcmd, flush=True)
    out = cnvnator2vcf_cmd.run()
    unfixedfile.write_text(out.stdout)

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
            "fix": lambda obj: HeaderContig(ID=chrom, length=size)
        })

    fix_vcffile(unfixedfile, outfile, fixes)


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
    cmd = cmdy.cnvpytor(
        root=rootfile,
        rd=bamfile,
        chrom=False if not chrom else chrom,
        **case,
        **cmdy_args,
    ).h()
    print()
    print(cmd.strcmd, flush=True)
    cmd.fg().run()

    # predicting cnv
    cmd = cmdy.cnvpytor(
        root=rootfile,
        his=binsizes,
        **cmdy_args,
    ).h()
    print()
    print(cmd.strcmd, flush=True)
    cmd.fg().run()

    cmd = cmdy.cnvpytor(
        root=rootfile,
        partition=binsizes,
        **cmdy_args,
    ).h()
    print()
    print(cmd.strcmd, flush=True)
    cmd.fg().run()

    if snp is not False:
        snp["chrom"] = snp.get("chrom", chrom)
        snp["sample"] = False if not snp.get("sample", False) else snp["sample"]
        mask_snps = snp.pop("mask_snps", True)
        baf_nomask = snp.pop("baf_nomask", False)
        cmd = cmdy.cnvpytor(
            root=rootfile,
            snp=snpfile,
            **snp,
            **cmdy_args,
        ).h()
        print()
        print(cmd.strcmd, flush=True)
        cmd.fg().run()

        if mask_snps:
            cmd = cmdy.cnvpytor(
                root=rootfile,
                mask_snps=True,
                **cmdy_args,
            ).h()
            print()
            print(cmd.strcmd)
            cmd.fg().run()

        cmd = cmdy.cnvpytor(
            root=rootfile,
            baf=binsizes,
            nomask=baf_nomask,
            **cmdy_args,
        ).h()
        print()
        print(cmd.strcmd, flush=True)
        cmd.fg().run()

    # call
    for binsize in binsizes:
        outfile = outdir / f"calls{'.combined' if snp is not False else ''}.{binsize}.tsv"
        print()
        print(cmdy.cnvpytor(
            root=rootfile,
            call="combined" if snp is not False else True,
            _=binsize,
            **cmdy_args,
        ).h.strcmd, flush=True)
        cmdy.cnvpytor(
            root=rootfile,
            call="combined" if snp is not False else True,
            _=binsize,
            **cmdy_args,
        ).r > outfile
        cnvpytor2other(outfile, bool(snp), "gff")
        cnvpytor2other(outfile, bool(snp), "bed")
        cnvpytor2vcf(outfile, bool(snp))

        # plots
        manplot = outdir / f"manhattan.{binsize}.png"
        cmd = cmdy.cnvpytor(
            root=rootfile,
            plot=["manhattan", binsize],
            chrom=False if not chrom else chrom,
            o=manplot,
            **cmdy_args,
        ).h()
        print()
        print(cmd.strcmd, flush=True)
        cmd.fg().run()

        circplot = outdir / f"circular.{binsize}.png"
        cmd = cmdy.cnvpytor(
            root=rootfile,
            plot=["circular", binsize],
            chrom=False if not chrom else chrom,
            o=circplot,
            _raise=False, # ZeroDivisionError
            **cmdy_args,
        ).h()
        print()
        print(cmd.strcmd, flush=True)
        cmd.fg().run()


do_case()
