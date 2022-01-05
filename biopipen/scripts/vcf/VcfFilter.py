from cyvcf2 import VCF, Writer, Variant

infile = {{in.invcf | repr}}
outfile = {{out.outfile | repr}}

{{envs.helper}}

keep = {{envs.keep | repr}}
filters = {{envs.filters | repr}}
filter_descs = {{envs.filter_descs | repr}}

# builtin filters
BUILTIN_FILTERS = {}

def builtin_filters(func):
    BUILTIN_FILTERS[func.__name__] = func
    return func

@builtin_filters
def SNPONLY(variant: Variant, nonrev: bool = True):
    """Keep or remove SNPs only"""
    ret = (
        len(variant.REF) == 1 and
        all(len(alt) == 1 for alt in variant.ALT)
    )
    return nonrev and ret

@builtin_filters
def QUAL(variant: Variant, cutoff, nonrev: bool = True):
    """Filter variants with QUAL above or below cutoff"""
    ret = variant.QUAL >= cutoff
    return nonrev and ret

for name, filt in filters.items():
    if name in BUILTIN_FILTERS:
        if not isinstance(filt, tuple):
            filt = (filt, )
        filters[name] = lambda variant: BUILTIN_FILTERS[name](variant, *filt)
        filters[name].__doc__ = BUILTIN_FILTERS[name].__doc__
    else:
        filters[name] = eval(filt)
        filters[name].__doc__ = filter_descs.get(name, filt)


invcf = VCF(infile)
for name, filt in filters.items():
    invcf.add_filter_to_header({
        'ID': name,
        'Description': filt.__doc__,
    })

if outfile.endswith(".gz"):
    outvcf = Writer(outfile, invcf, "wz")
else:
    outvcf = Writer(outfile, invcf)

for variant in invcf:
    for name, filt in filters.items():
        if not filt(variant):
            if not variant.FILTER:
                variant.FILTER = name
            else:
                variant.FILTER = f"{variant.FILTER};{name}"
    if variant.FILTER and not keep:
        continue
    outvcf.write_record(variant)

invcf.close()
outvcf.close()
