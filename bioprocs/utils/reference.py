from os import path, symlink
from bioprocs.utils import runcmd

def check(ref):
    if not ref or not path.exists(ref):
        raise Exception('Reference file not exists: %s' % ref)

def checkIndex(refindex):
    if not refindex:
        refindex = [refindex]
    return all([path.exists(ri) for ri in refindex])

def buildIndex(ref, cmd, ref2 = None, cmd2 = None):
    try:
        runcmd(cmd)
        return ref
    except Exception:
        try:
            if not path.exists(ref2):
                symlink(ref, ref2)
            runcmd(cmd2)
            return ref2
        except Exception:
            return None
