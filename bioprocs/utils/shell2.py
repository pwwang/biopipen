"""Shell utilities using cmdy"""
from modkit import modkit
import cmdy

# pylint: disable=invalid-name

DEFAULT_CONFIG = dict(
    default={'raise': True},

    bbmap_repair=dict(prefix='', sep='='),
    bedtools=dict(prefix='-'),
    biobambam=dict(sep='=', prefix=''),
    bowtie2=dict(dupkey=True),
    dtoxog=dict(prefix='-'),
    sort=dict(sep='', dupkey=True),
    gatk3=dict(dupkey=True),
    hla_la=dict(deform=None),
    liftover=dict(prefix='-', sep='='),
    oncotator=dict(sep='auto'),
    optitype=dict(dupkey=False),
    maf2vcf=dict(sep=' '),
    netmhc=dict(prefix='-'),

    # As of picard 2.20.5-SNAPSHOT
    # it's changing in the futher.
    # pylint: disable=line-too-long
    # See: https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
    # Future one should be:
    # picard = dict(sep = ' ', prefix = '-')
    picard=dict(sep='=', prefix=''),
    plink=dict(dupkey=False),
    pyclone=dict(deform=None),
    razers3=dict(prefix='-'),
    snpeff=dict(prefix='-'),
    vcfanno=dict(prefix='-'),
    vep=dict(dupkey=True, deform=None),
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
        conf2load[key].update(val if isinstance(val, dict) else {'exe': val})
    cmdy.CMDY_CONFIG._load(conf2load)

@modkit.delegate
def _modkit_delegate(module, name): # pylint: disable=unused-argument
    return getattr(cmdy, name)
