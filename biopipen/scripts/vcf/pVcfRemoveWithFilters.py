from cyvcf2 import VCF, Writer
from pyppl.utils import always_list

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
filters = {{args.filters | repr}}
reverse = {{args.reverse | bool}}
filters = always_list(filters)

reader = VCF(infile)
writer = Writer(outfile, reader)

for rec in reader:
    filt = rec.FILTER or 'PASS'
    filt = filt.split(';')
    hit = any(fil in filters for fil in filt)
    if (not reverse and hit) or (reverse and not hit):
        pass
    else:
        writer.write_record(rec)

reader.close()
writer.close()
