import json
from os import path
from pathlib import Path
from glob import glob
from textwrap import dedent
from .defaults import ROOT_DIR

def _text_dedent(text):
    splits = text.split('\n', 1)
    if len(splits) > 1:
        first_line, rest = splits
        return f"{first_line}\n{dedent(rest)}"
    return splits[0]

def always_list(value, trim=True):
    if isinstance(value, (tuple, list)):
        return list(value)
    return [val if not trim else val.strip() for val in value.split(',')]

def readlines(file):
    with open(file) as ff:
        return list(ff.readlines())

basename = lambda x: Path(x).name
stem = lambda x: Path(x).stem
stem2 = lambda x: Path(x).stem.split('.')[0]
abspath = lambda x: Path(x).resolve()
joinpath = lambda x, *args: path.join(x, *args)
ext = lambda x: Path(x).suffix
is_dir = lambda x: path.exists(x) and path.isdir(x)
quote = lambda x: json.dumps(str(x))

def r(var, ignoreintkey=True):
    """Convert a value into R values"""
    if var is True or var is False:
        # TRUE, FALSE
        return str(var).upper()
    if var is None:
        # NULL
        return 'NULL'
    if isinstance(var, str):
        # str representations
        if var.upper() in ['+INF', 'INF']:
            return 'Inf'
        if var.upper() == '-INF':
            return '-Inf'
        if var.upper() == 'TRUE':
            return 'TRUE'
        if var.upper() == 'FALSE':
            return 'FALSE'
        if var.upper() == 'NA' or var.upper() == 'NULL':
            return var.upper()
        if var.startswith('r:') or var.startswith('R:'):
            return str(var)[2:]
        return repr(str(var))
    if isinstance(var, Path):
        # stringify paths
        return repr(str(var))
    if isinstance(var, (list, tuple, set)):
        # All iterables to c()
        return 'c({})'.format(','.join([r(i) for i in var]))
    if isinstance(var, dict):
        return 'list({})'.format(','.join([
            '`{0}`={1}'.format(k, r(v)) if isinstance(k, int) and not ignoreintkey
            else r(v) if isinstance(k, int) and ignoreintkey
            # list allow repeated names
            else '`{0}`={1}'.format(str(k).split('#')[0], r(v))
            for k, v in sorted(var.items())
        ]))
    return repr(var)

def rimport(*scripts):
    """Load R scripts from script directory in the template"""
    # convert scripts into a dict which then can be transformed into R list
    scripts = {i: script for i, script in enumerate(scripts)}
    rsource = f'''
    if (!exists('.rimport') || !is.function(.rimport)) {{
        .rimport = function(...) {{
            for (rfile in list(...)) {{
                source(file.path('{ROOT_DIR}', 'utils', rfile))
            }}
        }}
    }}
    do.call(.rimport, {r(scripts)})
    '''
    return _text_dedent(rsource)

envs = dict(
    r=r,
    R=r,
    basename=basename,
    bn=basename,
    rimport=rimport,
    is_dir=is_dir,
    always_list=always_list,
    readlines=readlines,
    stem=stem,
    stem2=stem2,
    ext=ext,
    fn=stem,
    fn2=stem2,
    quote=quote,
    abspath=abspath,
    joinpath=joinpath,
    path_exists=path.exists,
    glob=glob
)
