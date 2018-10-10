from pyppl import Box
from os import path, rename, rmdir, system
from bioprocs.utils import runcmd, cmdargs, logger

infile     = {{i.manifile | quote}}
outdir     = {{o.outdir | quote}}
params     = {{args.params}}
nthread    = {{args.nthread | int}}
token      = {{args.token | repr}}
gdc_client = {{args.gdc_client | repr}}

# run gdc-client to download the data
gdc  = '{} download '.format(gdc_client)
args = Box({
	'm': infile,
	'n': nthread,
	'd': outdir,
	'retry-amount': '3',
	'debug': True,
	'log-file': path.join(outdir, 'gdc-client.log')
})
if token:
	args.t = token

args.update(params)
cmd2run = gdc + cmdargs(args, equal = ' ')
system(cmd2run)

# check if all the data sucessfully downloaded
with open(infile) as fin:
	ids = [line.split()[0] for line in fin if line.strip() and not line.startswith('id')]

del args['m']
for i in ids:
	if not path.isdir(path.join(outdir, i)):
		logger.warning('File failed to download: {}'.format(i))
		cmd2run = gdc + cmdargs(args, equal = ' ') + ' ' + i
		runcmd(cmd2run)
		rename(path.join(outdir, 'https:', 'api.gdc.cancer.gov', 'data', i), path.join(outdir, i))

if path.isdir(path.join(outdir, 'https:')):
	rmdir(path.join(outdir, 'https:', 'api.gdc.cancer.gov', 'data'))
	rmdir(path.join(outdir, 'https:', 'api.gdc.cancer.gov'))
	rmdir(path.join(outdir, 'https:'))
