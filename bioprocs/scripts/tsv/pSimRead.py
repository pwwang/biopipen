from os import path
from pyppl import Box
from bioprocs.utils import regionOverlap
from bioprocs.utils.tsvio import SimRead, TsvWriter

infiles = {{in.infiles}}
inopts  = {{args.inopts}}
usemeta = {{args.usemeta}} # int only
outfile = {{out.outfile | quote}}

outopts = Box(delimit = '\t', headPrefix = '', headDelimit = '\t', headTransform = None, head = False, ftype = '', cnames = [])
outopts.update({{args.outopts}})
head          = outopts['head']
headPrefix    = outopts['headPrefix']
headDelimit   = outopts['headDelimit']
headTransform = outopts['headTransform']
del outopts['head']
del outopts['headPrefix']
del outopts['headDelimit']
del outopts['headTransform']
#head = head if usemeta is None else True

writer  = TsvWriter(outfile, **outopts)
sm = SimRead(*infiles, **inopts)

if usemeta is not None:
	writer.meta.prepend(*sm.readers[usemeta - 1].meta.items())

{{args.helper}}
do      = {{args.do}}
match   = {{args.match}}

sm.do   = do
if callable(match):
	sm.match = match

if head:
	writer.writeHead(delimit = headDelimit, prefix = headPrefix, transform = headTransform)
sm.run()