from __future__ import annotations

from os import path, PathLike
from biopipen.core.filters import dict_to_cli_args
from biopipen.utils.reference import tabix_index
from biopipen.utils.misc import run_command

invcf: str | PathLike = {{in.invcf | quote}}  # noqa: E999 # pyright: ignore
outprefix: str = {{in.invcf | stem0 | quote}} # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
args: dict = {{envs | dict}}  # pyright: ignore

plink = args.pop("plink")
tabix = args.pop("tabix")
ncores = args.pop("ncores")

# normalize vcf-filter
args.setdefault("vcf_filter", True)
if isinstance(args["vcf_filter"], str):
    args["vcf_filter"] = args["vcf_filter"].split()

# normalize biallelic-only
args.setdefault("max_alleles", 2)

# This makes it possible to keep the allele order in the output
# no need for plink2
# args["keep_allele_order"] = True
args.setdefault("keep_allele_order", True)

# resolve plink 1.x --set-missing-var-ids doesn't distinguish $1, $2,...
# for ref and alts
# if (
#     "set_missing_var_ids" in args
#     and args["set_missing_var_ids"]
#     and ("$" in args["set_missing_var_ids"] or "%" in args["set_missing_var_ids"])
# ):
#     tmpfile = path.join(outdir, 'with_var_ids.vcf')
#     set_missing_var_ids = args.pop("set_missing_var_ids")
#     set_missing_var_ids = (
#         set_missing_var_ids
#         .replace("@", "%CHROM")
#         .replace("#", "%POS")
#         .replace("$1", "%REF")
#         .replace("$2", "%ALT{0}")
#         .replace("$3", "%ALT{1}")
#         .replace("$4", "%ALT{2}")
#         .replace("$5", "%ALT{3}")
#         .replace("$6", "%ALT{4}")
#         .replace("%CHROM_", "%CHROM\\_")
#         .replace("%POS_", "%POS\\_")
#         .replace("%REF_", "%REF\\_")
#     )
#     set_vid_cmd = [
#         bcftools,
#         "annotate",
#         "--set-id",
#         f"+{set_missing_var_ids}",
#         "--output-type",
#         "z",
#         "--output",
#         tmpfile,
#         "--threads",
#         ncores,
#         invcf,
#     ]

#     run_command(set_vid_cmd, fg=True, env={"cwd": outdir})
#     invcf = tmpfile

invcf = tabix_index(invcf, "vcf", tabix=tabix)
args["vcf"] = invcf
args["out"] = path.join(outdir, outprefix)
args["threads"] = ncores

cmd = [
    plink,
    "--make-bed",
    *dict_to_cli_args(args, dup_key=False, dashify = True),
]

run_command(cmd, fg=True, env={"cwd": outdir})
