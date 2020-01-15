from diot import Diot
from bioprocs.utils.tsvio2 import TsvReader

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
inopts  = {{args.inopts | repr}}
col     = {{args.col | repr}}
unique  = {{args.unique | repr}}

reader = TsvReader(infile, **inopts)

data = reader.dump(col = col)

if unique:
    data2 = []
    for dat in data:
        if dat not in data2:
            data2.append(dat)
    data = data2

with open(outfile, 'w') as fout:
    for dat in data:
        fout.write(dat + '\n')
