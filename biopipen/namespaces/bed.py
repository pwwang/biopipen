"""Tools to handle BED files"""
from ..core.proc import Proc
from ..core.config import config

class BedLiftOver(Proc):
    """Liftover a BED file using liftOver

    Input:
        inbed: The input BED file

    Output:
        outbed: The output BED file

    Envs:
        liftover: The path to liftOver
        chain: The map chain file for liftover
    """
    input = "inbed:file"
    output = "outbed:file:{{in.inbed | basename}}"
    envs = {
        "liftover": config.exe.liftover,
        "chain": config.path.liftover_chain,
    }
    lang = config.lang.bash
    script = "file://../scripts/bed/BedLiftOver.sh"
