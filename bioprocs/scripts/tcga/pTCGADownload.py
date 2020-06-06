"""Script for tcga.pTCGADownload"""
# pylint: disable=undefined-variable,unused-import,invalid-name,used-before-assignment

from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell, logger

infile = {{i.manifile | quote}}
outdir = {{o.outdir | quote}}
params = {{args.params}}
nthread = {{args.nthread | int}}
token = {{args.token | repr}}
gdc_client = {{args.gdc_client | repr}}

shell.load_config(gdc=gdc_client)

dftargs = Diot({
    'retry-amount': '3',
    'debug': False,
    'log-file': path.join(outdir, 'gdc-client.log')
})
dftargs.update(params)
dftargs.m = infile
dftargs.n = nthread
dftargs.d = outdir
if token:
    dftargs.t = token
shell.gdc.download(**dftargs).fg

# check if all the data sucessfully downloaded
logger.warning(
    'Checking if all files in manifest has been sucessfully downloaded ...')
with open(infile) as fin:
    ids = [line.split()[0] for line in fin if line.strip()
           and not line.startswith('id')]
logger.warning('Got %s file ids', len(ids))
logger.warning('%s files failed', sum(
    int(not path.isdir(path.join(outdir, i))) for i in ids))

del dftargs['m']
for i in ids:
    if not path.isdir(path.join(outdir, i)):
        logger.warning('FAILED: %s', i)
        logger.warning('- Trying to find it in https directory ...')
        destfile = path.join(outdir, 'https:', 'api.gdc.cancer.gov', 'data', i)
        if path.exists(destfile):
            logger.warning('- Found, save it.')
            shell.mv(destfile, path.join(outdir, i))
        else:
            logger.warning('- Not found, try to download it.')
            dftargs['_'] = [i, 'shraw:1>&2']
            gdc.download(**dftargs).run()
            shell.mv(destfile, path.join(outdir, i))

if path.isdir(path.join(outdir, 'https:')):
    shell.rm_rf(path.join(outdir, 'https:'))
