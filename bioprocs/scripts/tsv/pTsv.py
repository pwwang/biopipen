"""Script for tsv.pTsv"""
# pylint: disable=unused-import,exec-used,undefined-variable,not-callable
# pylint: disable=unsupported-assignment-operation,not-a-mapping,invalid-name
from diot import Diot
from bioprocs.utils.tsvio2 import TsvWriter, TsvReader, TsvRecord

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts = {{args.inopts | repr}}
outopts = {{args.outopts | repr}}
helper = {{args.helper | repr}}
if not isinstance(helper, list):
    helper = [helper]

helper = [line for line in helper if line]
inopts['row'] = {{args.inopts.get('row', None)}}

reader = TsvReader(infile, **inopts)
exec('\n'.join(helper), globals())
row_func = {{args.row | render}}

writer = TsvWriter(outfile, delimit=outopts.get('delimit', "\t"))
outcnames = outopts.get('cnames', True)
head_callback = {{args.outopts.get('headCallback', True)}}
if outcnames is True:
    writer.cnames = reader.cnames
    writer.writeHead(head_callback)
elif isinstance(outcnames, (list, tuple)):
    writer.cnames = list(outcnames)
    writer.writeHead(head_callback)

for record in reader:
    rec = row_func(record)
    if rec is False:
        continue
    if rec is None or rec is True:
        rec = record
    writer.write(rec)
writer.close()
