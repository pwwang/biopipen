#!/usr/bin/env python
"""Plot genomic elements."""
import json, sys
from os import path
from pyppl import PyPPL, utils, Box
from bioprocs import params

params._prefix = '-'

params.ideo       = 'True'
params.ideo.desc  = 'Show ideogram track?'
params.genes      = params.refgene.value
params.genes.desc = 'Show gene track?'
params.axis       = True
params.axis.desc  = 'Show axis?'

params.genome.show     = True
params.outdir.required = True
params.outdir.desc     = 'Output directory.'

params.region.required = True
params.region.desc     = 'The region to plot. E.g. chr1:1000000-1200000'
params.tracks.required = True
params.tracks.type     = list
params.tracks.desc     = 'The track types. Could be data, anno, interaction or ucsc, or multiple of them.'
params.names.required  = True
params.names.type      = list
params.names.desc      = 'The corresponding names of the tracks.'
params.inputs.required = True
params.inputs.type     = list
params.inputs.desc     = 'The input of the tracks.\n"<ucscTrack>:<gvizTrack>" for ucsc track;\n"<infile>:<intype>" for interaction tracks;\nfiles for other tracks.'
params.params.type     = list
params.params.desc     = 'The params for each track'
params.plotparams.desc = 'The params for pGenomePlot.'
params.devpars         = '{"res": 300, "height": 300, "width": 2000}'
params.devpars.desc    = 'The device parameters for plotting.'
params.splitlen        = 5000000
params.splitlen.desc   = 'Split the plot into 2 if the region is longer then splitlen.'
params.forks           = 2
params.forks.desc      = 'Number of cores used to plot if split.'

# highlist
params.highlights      = []
params.highlights.desc = 'The highlight regions in format of "start-end"'

def main():
	opts = params._parse(dict_wrapper = Box)
	from bioprocs.genomeplot import pInteractionTrack, pAnnoTrack, pDataTrack, \
		pUcscTrack, pGenomePlot

	chrom, startend  = opts.region.split(':')
	start, end       = startend.split('-')
	start      = int(start)
	end        = int(end)
	trackProcs = []
	#uuid       = utils.uid(str(sys.argv))
	for i, tt in enumerate(opts.tracks):
		if tt == 'data':
			datatrackproc = pDataTrack.copy(tag = opts.names[i])
			datatrackproc.input = (opts.names[i], opts.inputs[i], chrom)
			if opts.params:
				datatrackproc.args.opts.update(json.loads(opts.params[i]))
			trackProcs.append(datatrackproc)
		elif tt == 'anno':
			annotrackproc = pAnnoTrack.copy(tag = opts.names[i])
			annotrackproc.input = (opts.names[i], opts.inputs[i], chrom)
			if opts.params:
				annotrackproc.args.opts.update(json.loads(opts.params[i]))
			trackProcs.append(annotrackproc)
		elif tt == 'ucsc':
			ucsctrackproc = pUcscTrack.copy(tag = opts.names[i])
			ucsctrack, gviztrack = opts.inputs[i].split(':')
			ucsctrackproc.input = (opts.names[i], ucsctrack, gviztrack, opts.region)
			if opts.params:
				ucsctrackproc.args.opts.update(json.loads(opts.params[i]))
			trackProcs.append(ucsctrackproc)
		else:
			intertrackproc = pInteractionTrack.copy(tag = opts.names[i])
			infile, intype = opts.inputs[i].split(':')
			intertrackproc.input = (opts.names[i], infile, opts.region)
			intertrackproc.args.intype = intype
			if opts.params:
				intertrackproc.args.opts.update(json.loads(opts.params[i]))
			trackProcs.append(intertrackproc)

	if end - start > opts.splitlen:
		pGenomePlot.depends        = trackProcs
		#pGenomePlot.tag            = uuid
		pGenomePlot.forks          = opts.forks
		pGenomePlot.exdir          = opts.outdir
		pGenomePlot.args.ideoTrack = opts.ideo
		pGenomePlot.args.axisTrack = opts.axis
		pGenomePlot.args.geneTrack = opts.genes
		if opts.devpars:
			pGenomePlot.args.devpars.update(json.loads(opts.devpars))
		if opts.plotparams:
			pGenomePlot.args.opts.update(json.loads(opts.plotparams))
		if len(opts.highlights) == 2 and '-' not in opts.highlist[0]:
			h1 = opts.highlights[0]
			h2 = opts.highlights[1]
		else:
			h1 = ';'.join(opts.highlights)
			h2 = h1
		pGenomePlot.input   = lambda *chs: [
			([ch.get() for ch in chs], "%s:%s-%s" % (chrom, start, start + 10000), h1),
			([ch.get() for ch in chs], "%s:%s-%s" % (chrom, end - 100000, end), h2),
		]

	else:
		pGenomePlot.depends        = trackProcs
		#pGenomePlot.tag            = uuid
		pGenomePlot.exdir          = opts.outdir
		pGenomePlot.args.ideoTrack = opts.ideo
		pGenomePlot.args.axisTrack = opts.axis
		pGenomePlot.args.geneTrack = opts.genes
		if opts.devpars:
			pGenomePlot.args.devpars.update(json.loads(opts.devpars))
		if opts.plotparams:
			pGenomePlot.args.opts.update(json.loads(opts.plotparams))
		pGenomePlot.input   = lambda *chs: [([ch.get() for ch in chs], opts.region, ';'.join(opts.highlights))]

	config = {'proc': {'file': None}}
	PyPPL(config).start(trackProcs).run()

if __name__ == "__main__":
	main()
