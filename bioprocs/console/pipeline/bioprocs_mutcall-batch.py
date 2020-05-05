"""
Call mutations in batches
We can use `bioprocs mutcall` to do it, but since we have limited space,
we need to split the data into batches and remove the intermediate
files for each batch.
"""
import math
from pathlib import Path
import toml
import cmdy
from diot import Diot
from bioprocs import params
from bioprocs.utils import logger
from bioprocs.utils.sampleinfo import SampleInfo2
from bioprocs.utils.parallel import distribute_list
from bioprocs.utils.tsvio2 import TsvWriter#, TsvReader

params._desc = __doc__
params.metadir = './<ppldir>/<saminfo>.meta'
params.metadir.desc = 'The meta directory to save cache and log files'
params.metadir.callback = (
    lambda opt, ps:
    opt.set_value(str(Path(ps.ppldir.value).joinpath(
        '%s.meta' % Path(ps.saminfo.value).stem
    )))
    if '<ppldir>' in opt.value
    else None
)
params.indir.required = True
params.indir.desc = 'The data directory with raw data'
params.saminfo.required = True
params.saminfo.desc = 'The sample information'
params.batchsize = 10
params.batchsize.desc = ['The batch size (how many instances to run per batch)',
                         'Note if samples are paired (tumor-normal), '
                         'this will be number of pairs']
params.outdir.required = True
params.outdir.desc = 'The output directory'
params.runner = 'local'
params.runner.desc = 'The runner for the pipeline'
params.ppldir = './workdir'
params.ppldir.desc = ('The directory for PyPPL pipeline where the intermediate '
                      'files are saved. This directory will be clearned up if '
                      '`dropim` is `True`')
params.intype = 'bam' # ebam, fastq
params.intype.desc = ['The type of input files. Could be ebam, bam or fastq',
                      '- ebam: External bam files, need to be remapped to '
                      'current ref genome',
                      '- bam: Bam files already mapped to current ref genome',
                      '- fastq: Raw read files']
params.muts = ['germ'] # germ, soma, scnv, gcnv
params.muts.desc = [
    'What kind of mutations to call.',
    'Note: soma need paired information. scnv needs `batchsize` > 1.'
]
params.aligner = 'bwa'
params.aligner.desc = 'The alignment tool'
params.trimmer = 'trimmomatic'
params.trimmer.desc = 'The trimming tool'
params.compress = True
params.compress.desc = 'Use gzip and bam file to save space.'
params.dropim = False
params.dropim.desc = 'Should we drop the intermediate files?'
params.nthread = 1
params.nthread.desc = 'Threads used for each job'
params.forks = '<batchsize>*2'
params.forks.desc = 'How many instances to run simultaneously for each batch.'
params.forks.callback = (lambda opt, ps: opt.set_value(ps.batchsize.value * 2)
                         if '<batchsize>' in opt.value else None)
params.forcecache = True
params.forcecache.desc = ('Force processes to cache if they had successful '
                          'run before without rcfiles being generated.')

def parse_and_print_args():
    """Parse and print the command-line arguments"""
    opts = params._parse(dict_wrapper=Diot)
    logger.info('Arguments parsed:')
    for key in opts:
        if params[key].show and key not in params._hopts:
            logger.info('- %-10s: %s', key, opts[key])

    opts.metadir = Path(opts.metadir)
    if not opts.metadir.is_dir():
        opts.metadir.mkdir(parents=True, exist_ok=True)
    opts.saminfo = Path(opts.saminfo)
    opts.ppldir = Path(opts.ppldir)
    opts.cachefile = opts.metadir / (('%s-cache.txt') % opts.saminfo.stem)
    return opts

def batchdir(ppldir, saminfo, batch):
    """Get batch directory with its workdirs for processes"""
    return ppldir.joinpath('%s-batch%s' % (saminfo.stem, batch))

def clean_im(wdir, nthread):
    """Clean intermediate files (process workdirs) in wdir"""
    cmdy.pyppl.clean(
        _raise=False,
        _fg=True,
        nthread=nthread,
        nocheck=True,
        force=True,
        wdir=wdir
    )
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

def main():
    """Main function"""
    opts = parse_and_print_args()
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

        cmd = cmdy.bioprocs.mutcall(
            _hold=True,
            _fg=True,
            _raise=True,
            intype=opts.intype,
            muts=opts.muts,
            indir=opts.indir,
            aligner=opts.aligner,
            trimmer=opts.trimmer,
            saminfo=batch_saminfo,
            nthread=opts.nthread,
            outdir=opts.outdir,
            runner=opts.runner,
            forks=opts.forks,
            logfile=logfile,
            compress=opts.compress,
            forcecache=opts.forcecache,
            ppldir=batchdir(opts.ppldir, opts.saminfo, i)
        )
        logger.info(cmd.cmd)
        cmd.run()
        logger.info('')
        if opts.dropim:
            clean_im(batchdir(opts.ppldir, opts.saminfo, i), opts.forks)
        update_cache(i + 1, opts.cachefile)

if __name__ == '__main__':
    main()
