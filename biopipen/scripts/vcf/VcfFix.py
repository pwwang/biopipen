# for fixes to evaluate
from biopipen.scripts.vcf.VcfFix_utils import (  # noqa: F401
    HeaderItem,
    HeaderInfo,
    HeaderFormat,
    HeaderFilter,
    HeaderContig,
    HeaderGeneral,
    Fields,
    Info,
    Format,
    Alt,
    Filter,
    Sample,
    Samples,
    Variant,
)
from biopipen.scripts.vcf.VcfFix_utils import fix_vcffile

infile = {{in.infile | quote}}  # pyright: ignore
instem = {{in.infile | stem | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore

{{envs.helpers}}  # pyright: ignore
fixes = []
{%- for fix_item in envs.fixes %}  # pyright: ignore

# New fix item
fix = {}
{%-     for key, value in fix_item.items() %}  # pyright: ignore
{%-          if key == "fix" %}  # pyright: ignore
fix[{{key | quote}}] = {{value}}  # pyright: ignore
{%-          else %}  # pyright: ignore
fix[{{key | quote}}] = {{value | repr}}  # pyright: ignore
{%-          endif %}  # pyright: ignore
{%-     endfor %}  # pyright: ignore
fixes.append(fix)
{%- endfor %}  # pyright: ignore

fix_vcffile(infile, outfile, fixes)
