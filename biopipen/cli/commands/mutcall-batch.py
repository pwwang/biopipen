"""
Call mutations in batches
We can use `biopipen mutcall` to do it, but since we have limited space,
we need to split the data into batches and remove the intermediate
files for each batch.
"""
# pylint: disable=invalid-name

import math
from pathlib import Path
import toml
import cmdy
from diot import Diot
from pyparam import Params
from biopipen import params as bp_params
from biopipen.utils import logger
from biopipen.utils.sampleinfo import SampleInfo2
from biopipen.utils.parallel import distribute_list
from biopipen.utils.tsvio2 import TsvWriter#, TsvReader


HERE = Path(__file__).parent.resolve()
ARGS_FILE = HERE / 'mutcall-batch.args.toml'

help_group = 'SCRIPTS'
params = Params(desc=__doc__)
params.from_file(ARGS_FILE)

params.get_param('metadir').callback = lambda val, allvals: (
    val or allvals.ppldir.joinpath(allvals.saminfo.stem + '.meta')
)
params.get_param('forks').callback = lambda val, allvals: (
    allvals.batchsize * 2 if '<batchsize>' in val else None
)

def parse_and_print_args(opts):
    """Parse and print the command-line arguments"""
    logger.info('Arguments parsed:')
    for key in opts:
        if params.get_param(key).show:
            logger.info('- %-10s: %s', key, opts[key])

    if not opts.metadir.is_dir():
        opts.metadir.mkdir(parents=True, exist_ok=True)
    opts.cachefile = opts.metadir / (('%s-cache.txt') % opts.saminfo.stem)
    return bp_params.parse() | opts

def batchdir(ppldir, saminfo, batch):
    """Get batch directory with its workdirs for processes"""
    return ppldir.joinpath('%s-batch%s' % (saminfo.stem, batch))

def clean_im(wdir, nthread):
    """Clean intermediate files (process workdirs) in wdir"""
    cmdy.pyppl.clean(
        cmdy_raise=False,
        nthread=nthread,
        nocheck=True,
        force=True,
        wdir=wdir
    ).fg()
    cmdy.rm(wdir, r=True, f=True)

def split_batches(saminfo, batchsize):
    """Split sample(pair)s into batches"""
    saminfo = SampleInfo2(saminfo)
    ret = []
    if saminfo.is_paired():
        # make sure pairs in same batch
        # pair-mates must be defined with the same patient
        # different pairs must be defined with different patients,
        # even they are from the same patients
        for patient in saminfo.all_patients():
            ret.append(tuple(saminfo.get_samples(by='Patient',
                                                 value=patient,
                                                 return_all=True)))
    else:
        # samples should be unique then
        ret = saminfo.get_samples(return_all=True)
    total = len(ret)
    ret = list(distribute_list(ret, math.ceil(total / batchsize)))
    logger.info('Splitting into %d batchs with batchsize %s for total %s',
                len(ret), batchsize, total)
    return ret

def get_current_batch(saminfo, cachefile, metadir, batchsize):
    """Get batch number and sample info file for current batch"""
    logger.info('Fetching current batch ...')
    cache = False
    if cachefile.is_file():
        with cachefile.open() as fcache:
            cache = toml.load(fcache, Diot)

        if (not Path(cache.saminfo).samefile(saminfo) or
                cache.batchsize != batchsize):
            cache = False
        else:
            #curr_batch = batches[cache.batch]
            batch_saminfo = metadir.joinpath('%s-batch%s.saminfo' %
                                             (saminfo.stem, cache.batch))
            # check if samples are the same
            if not batch_saminfo.is_file():
                cache = False
            # else:
            #     cached_samples = TsvReader(batch_saminfo).dump("Sample")
            #     batch_samples = (sum(([sample[0].Sample, sample[1].Sample]
            #                            for sample in curr_batch), [])
            #                      if isinstance(batches[0][0], tuple)
            #                      else [sample.Sample
            #                      for sample in curr_batch])
            #     if cached_samples != batch_samples:
            #         cache = False
    else:
        cache0 = Diot(saminfo=str(saminfo), batchsize=batchsize, batch=0)
        with cachefile.open('w') as fcache:
            toml.dump(cache0, fcache)
    return cache.batch if cache else 0

def save_batch_saminfo(samrows, batch, saminfo, metadir):
    """Save sample info the current batch"""
    batch_saminfo = metadir.joinpath('%s-batch%s.saminfo' %
                                     (saminfo.stem, batch))

    if batch_saminfo.is_file():
        return batch_saminfo
    writer = TsvWriter(batch_saminfo)
    if isinstance(samrows[0], tuple): #paired
        writer.cnames = ['Sample', 'Patient', 'Group']
    else:
        writer.cnames = ['Sample', 'Group']
    writer.write_head()

    for samples in samrows:
        if isinstance(samples, tuple):
            writer.write([samples[0].Sample,
                          samples[0].Patient, samples[0].Group])
            writer.write([samples[1].Sample,
                          samples[1].Patient, samples[1].Group])
        else:
            writer.write([samples.Sample, samples.Patient, samples.Group])
    writer.close()
    return batch_saminfo

def update_cache(batch, cachefile):
    """Update current batch to cache."""
    with cachefile.open() as fcache:
        cache = toml.load(fcache, Diot)
    cache.batch = batch
    with cachefile.open('w') as fcache:
        toml.dump(cache.as_dict(), fcache)

def main(opts):
    """Main function"""
    opts = parse_and_print_args(opts)

    batches = split_batches(opts.saminfo, opts.batchsize)
    nbatches = len(batches)
    curr = get_current_batch(opts.saminfo,
                             opts.cachefile,
                             opts.metadir,
                             opts.batchsize)
    for i in range(curr):
        logger.info('### Batch %03d/%03d finished, %s',
                    i, nbatches - 1,
                    'trying to clean intermediate files and skip'
                    if opts.dropim
                    else 'skip.')
        if opts.dropim:
            clean_im(batchdir(opts.ppldir, opts.saminfo, i), opts.forks)

    for i in range(curr, len(batches)):
        logger.info('')
        logger.info('### Batch %03d/%03d started:', i, nbatches - 1)
        logfile = opts.metadir.joinpath('%s-batch%s.log' %
                                        (opts.saminfo.stem, i))
        logfile.write_text('')

        batch_saminfo = save_batch_saminfo(batches[i], i,
                                           opts.saminfo, opts.metadir)

        cmd = cmdy.biopipen.mutcall(
            _raise=True,
            intype=opts.intype,
            muts=opts.muts,
            indir=opts.indir,
            aligner=opts.aligner,
            trimmer=opts.trimmer,
            errhow=opts.errhow,
            saminfo=batch_saminfo,
            nthread=opts.nthread,
            outdir=opts.outdir,
            runner=opts.runner,
            forks=opts.forks,
            logfile=logfile,
            compress=opts.compress,
            forcecache=opts.forcecache,
            ppldir=batchdir(opts.ppldir, opts.saminfo, i)
        ).h.fg
        logger.info(cmd.cmd)
        cmd.run()
        logger.info('')
        if opts.dropim:
            clean_im(batchdir(opts.ppldir, opts.saminfo, i), opts.forks)
        update_cache(i + 1, opts.cachefile)

if __name__ == '__main__':
    main(params.parse())
