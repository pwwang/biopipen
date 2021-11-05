from pathlib import Path
import cmdy

bamfile = {{in.bamfile | quote}}
snpfile = {{in.snpfile | repr}}
outdir = Path({{out.outdir | quote}})
cnvpytor = {{envs.cnvpytor | quote}}
cases = {{envs.cases | repr}}
ncores = {{envs.ncores | int}}

cmdy_args = {"_exe": cnvpytor, "_prefix": "-", "_deform": None}
def do_case(casename, case):
    casedir = outdir / casename
    casedir.mkdir(exist_ok=True)

    chrom = case.pop("chrom", [])
    binsizes = case.pop("binsizes", [10000, 100000])
    snp = case.pop("snp", False)
    rootfile = casedir / "file.pytor"
    case["j"] = case.get("j", ncores)

    # read depth signal
    cmdy.cnvpytor(
        root=rootfile,
        rd=bamfile,
        chrom=False if not chrom else chrom,
        **case,
        **cmdy_args,
    )

    # predicting cnv
    cmdy.cnvpytor(
        root=rootfile,
        his=binsizes,
        **cmdy_args,
    ).fg
    cmdy.cnvpytor(
        root=rootfile,
        partition=binsizes,
        **cmdy_args,
    ).fg

    if snp is not False:
        snp["chrom"] = snp.get("chrom", chrom)
        snp["sample"] = False if not snp.get("sample", False) else snp["sample"]
        mask_snps = snp.pop("mask_snps", True)
        baf_nomask = snp.pop("baf_nomask", False)
        cmdy.cnvpytor(
            root=rootfile,
            snp=snpfile,
            **snp,
            **cmdy_args,
        ).fg

        if mask_snps:
            cmdy.cnvpytor(
                root=rootfile,
                mask_snps=True,
                **cmdy_args,
            ).fg

        cmdy.cnvpytor(
            root=rootfile,
            baf=binsizes,
            nomask=baf_nomask,
            **cmdy_args,
        ).fg

    # call
    for binsize in binsizes:
        outfile = casedir / f"calls{'.combined' if snp is not False else ''}.{binsize}.tsv"
        cmdy.cnvpytor(
            root=rootfile,
            call="combined" if snp is not False else True,
            _=binsize,
            **cmdy_args,
        ).r > outfile

        # plots
        manplot = casedir / f"manhattan.{binsize}.png"
        cmdy.cnvpytor(
            root=rootfile,
            plot=["manhattan", binsize],
            chrom=False if not chrom else chrom,
            o=manplot,
            **cmdy_args,
        ).fg

        circplot = casedir / f"circular.{binsize}.png"
        cmdy.cnvpytor(
            root=rootfile,
            plot=["circular", binsize],
            chrom=False if not chrom else chrom,
            o=circplot,
            _raise=False, # ZeroDivisionError
            **cmdy_args,
        ).fg


for casename, case in cases.items():
    do_case(casename, case)
