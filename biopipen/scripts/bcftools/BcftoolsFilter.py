import tempfile
from pathlib import Path
from hashlib import md5
import cmdy

infile = {{in.infile | repr}}
outfile = {{out.outfile | repr}}
bcftools = {{envs.bcftools | repr}}
keep = {{envs.keep | repr}}
args = {{envs.args | repr}}
ncores = {{envs.ncores | repr}}
tmpdir = {{envs.tmpdir | repr}}
includes = {{envs.includes | repr}}
excludes = {{envs.excludes | repr}}

args["_exe"] = bcftools
args["_"] = infile
args["o"] = outfile
args["threads"] = ncores
if "O" not in args and "output-type" not in args:
    args["O"] = "z" if infile.endswith(".gz") else "v"
if "m" not in args and "mode" not in args:
    args["m"] = "+"

tmpdir = (
    Path(tmpdir) / f"biopipen-bcftoolsfilter-{md5(infile.encode()).hexdigest()}"
)
tmpdir.mkdir(parents=True, exist_ok=True)
# a.vcf.gz -> a
# a.vcf -> a
stem = Path(infile).stem
if stem.endswith(".vcf"):
    stem = stem[:-4]
# .vcf.gz
# .gz
ext = Path(infile).name[len(stem):]

FILTER_INDEX = [1]

def handle_filter(vcf, fname, filt, flag):
    print("- Handling filter ", fname, ": ", filt, " ...")

    arguments = args.copy()
    arguments[flag] = filt
    arguments["_"] = vcf
    arguments["o"] = tmpdir / f"{stem}.{fname}{ext}"
    if keep:
        arguments["s"] = fname

    cmd = cmdy.bcftools.filter(**arguments).hold()
    print("  running:")
    print("  ", cmd.strcmd)
    cmd.run(wait=True)
    return arguments["o"]


def normalize_expr(expr, flag):
    out = {}
    if not expr:
        return out
    if isinstance(expr, list):
        for ex in expr:
            out[f"FILTER{FILTER_INDEX[0]}"] = (ex, flag)
            FILTER_INDEX[0] += 1
    elif isinstance(expr, dict):
        for name, ex in expr.items():
            out[name] = (ex, flag)
    else: # str
        out[f"FILTER{FILTER_INDEX[0]}"] = (expr, flag)
        FILTER_INDEX[0] += 1
    return out

includes = normalize_expr(includes, "include")
excludes = normalize_expr(excludes, "exclude")
includes.update(excludes)

# bcftools can be only done once at one filter
for fname, (filt, flag) in includes.items():
    infile = handle_filter(infile, fname, filt, flag)

cmdy.cp(infile, outfile)
