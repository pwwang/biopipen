import pandas
from datar.all import c, right_join, select, relocate
from biopipen.utils.gene import gene_name_conversion

infile = {{in.infile | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
inopts = {{envs.inopts | repr}}  # pyright: ignore
outopts = {{envs.outopts | repr}}  # pyright: ignore
notfound = {{envs.notfound | repr}}  # pyright: ignore
genecol = {{envs.genecol | repr}}  # pyright: ignore
output = {{envs.output | repr}}  # pyright: ignore
infmt = {{envs.infmt | repr}}  # pyright: ignore
outfmt = {{envs.outfmt | repr}}  # pyright: ignore
species = {{envs.species | quote}}  # pyright: ignore

df = pandas.read_csv(infile, **inopts)

if isinstance(genecol, int):
    genes = df.iloc[:, genecol]
else:
    genes = df.loc[:, genecol]

colname = genes.name
genes = genes.tolist()

#        query  `outfmt`
#     <object> <object>
# 0  1255_g_at   GUCA1A
# 1    1316_at     THRA
# 2    1320_at   PTPN21
# 3    1294_at  MIR5193
converted = gene_name_conversion(
    genes=genes,
    species=species,
    infmt=infmt,
    outfmt=outfmt,
    notfound=notfound,
)
converted.columns = [colname] + converted.columns[1:].tolist()

if output == "only":
    out = converted

elif output == "keep":
    out = df >> right_join(converted, by=colname, suffix=["", "_converted"])

elif output == "drop":
    out = df >> right_join(
        converted,
        by=colname, suffix=["", "_converted"]
    ) >> select(~c(colname))

elif output == "replace":
    out = df >> right_join(
        converted, by=colname, suffix=["", "_converted"]
    )
    converted_cols = out.columns[-len(converted.columns)+1:].tolist()
    pos = df.columns.get_indexer([colname])[0]
    out = out >> relocate(
        converted_cols, _after=pos+1
    ) >> select(~c(colname))

else:
    raise ValueError(f"Unknown output mode: {output}.")

out.to_csv(outfile, **outopts)
