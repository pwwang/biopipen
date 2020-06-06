"""Shell utilities using cmdy"""
from modkit import modkit

import cmdy

DEFAULT_CONFIG = dict(
    default=dict(_raise=True),

    bbmap_repair=dict(_prefix='', _sep='='),
    bedtools=dict(_prefix='-'),
    biobambam=dict(_sep='=', _prefix=''),
    bowtie2=dict(_dupkey=True),
    dtoxog=dict(_prefix='-'),
    sort=dict(_sep='', _dupkey=True),
    gatk3=dict(_dupkey=True),
    hla_la=dict(_deform=None),
    liftover=dict(_prefix='-', _sep='='),
    oncotator=dict(_sep='auto'),
    optitype=dict(_dupkey=False),
    maf2vcf=dict(_sep=' '),
    netmhc=dict(_prefix='-'),

    # As of picard 2.20.5-SNAPSHOT
    # it's changing in the futher.
    # pylint: disable=line-too-long
    # See: https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    # Future one should be:
    # picard = dict(_sep = ' ', _prefix = '-')
    picard=dict(_sep='=', _prefix=''),
    plink=dict(_dupkey=False),
    pyclone=dict(_deform=None),
    razers3=dict(_prefix='-'),
    snpeff=dict(_prefix='-'),
    vcfanno=dict(_prefix='-'),
    vep=dict(_dupkey=True, _deform=None),
)
cmdy.CMDY_CONFIG._load(DEFAULT_CONFIG)

# aliases
rm_rf = cmdy.rm.bake(r=True, f=True)
ln_s = cmdy.ln.bake(s=True)
kill_9 = cmdy.kill.bake(s=9)
wc_l = cmdy.wc.bake(l=True)
cp = copy = cmdy.cp
mv = move = cmdy.mv
which = lambda x: cmdy.which(x).strip()
runcmd = lambda cmd, **kwargs: cmdy.bash(c=cmd, **kwargs)

def load_config(conf=None, **kwargs):
    """Load the configurations for the module"""
    conf = conf or {}
    conf.update(kwargs)
    conf2load = {'default': DEFAULT_CONFIG['default']}
    for key, val in conf.items():
        conf2load[key] = DEFAULT_CONFIG.get(key, {}).copy()
        conf2load[key].update(val if isinstance(val, dict) else {'_exe': val})

    cmdy.CMDY_CONFIG._load(conf2load)

@modkit.delegate
def _modkit_delegate(module, name):
    return getattr(cmdy, name)
