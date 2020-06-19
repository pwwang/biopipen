"""Environment to render templates"""
import json
from pathlib import Path
from pyppl.template import DEFAULT_ENVS

__all__ = []

HERE = Path(__file__).parent.resolve()

def _get_paths(paths):
    ret = []
    for path in paths:
        if isinstance(path, (list, tuple)):
            ret.extend(list(path))
        else:
            ret.append(path)
    return (repr(str(path)) for path in ret)


def rimport(*paths):
    """Import an R source file"""
    rimport_rfunc = f"""
    if (!exists('..rimport..') || !is.function(..rimport..)) {{
        ..rimport.. = function(...) {{
            for (rfile in list(...)) {{
                source(file.path({str(HERE)!r}, 'utils', rfile))
            }}
        }}
    }}
    """
    pathstr = ', '.join(_get_paths(paths))
    return f"""
    {rimport_rfunc}
    ..rimport..({pathstr})
    """


def bashimport(*paths):
    """Import a bash file"""
    bashimport_bashfunc = f"""
    type __bashimport__ 1>&2 2>/dev/null
    if [ $? -ne 0 ]; then
        function __bashimport__() {{
            for src in "$@"; do
                source {HERE}/utils/$src
            done
        }}
    fi
    """
    pathstr = ', '.join(_get_paths(paths))
    return f"""
    {bashimport_bashfunc}
    __bashimport__ {pathstr}
    """

def commonprefix(*paths, basename_only=True):
    """Get the common prefix of a set of paths"""
    from os import path
    paths = [path.basename(pth) if basename_only else pth
             for pth in _get_paths(paths)]
    return path.commonprefix(paths)

def read(var):
    """Read the contents from a file"""
    with open(var) as fvar:
        return fvar.read()


def readlines(var, skip_empty_lines=True):
    """Read the lines from a file"""
    ret = []
    with open(var) as fvar:
        for line in fvar:
            line = line.rstrip('\n\r')
            if not line and skip_empty_lines:
                continue
            ret.append(line)
    return ret


def basename(var, orig=False):
    """Get the basename of a path"""
    bname = Path(var).name
    if orig or not bname.startswith('['):
        return bname

    return bname[bname.find(']') + 1:]


def filename(var, orig=False, dot=-1):
    """
    Return the stem of the basename (stripping extension(s))
    @params:
        `var`: The path
        `orig`: If the path is a renamed file (like: `origin[1].txt`),
            - whether return its original filename or
                the parsed filename (`origin.txt`)
        `dot`: Strip to which dot.
            - `-1`: the last one
            - `-2`: the 2nd last one ...
            - `1` : remove all dots.
    """
    bname = basename(var, orig)
    if '.' not in bname:
        return bname
    return '.'.join(bname.split('.')[0:dot])


def prefix(var, orig=False, dot=-1):
    """Get the prefix part of a path"""
    return str(Path(var).parent.joinpath(filename(var, orig, dot)))


def R(var, ignoreintkey=True): # pylint: disable=too-many-return-statements,invalid-name
    """Convert a value into R values"""
    if var is True or var is False:
        return str(var).upper()
    if var is None:
        return 'NULL'
    if isinstance(var, str):
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
        return repr(str(var))
    if isinstance(var, (list, tuple, set)):
        return 'c({})'.format(','.join([R(i) for i in var]))
    if isinstance(var, dict):
        # list allow repeated names
        return 'list({})'.format(','.join([
            '`{0}`={1}'.format(k, R(v))
            if isinstance(k, int) and not ignoreintkey
            else R(v)
            if isinstance(k, int) and ignoreintkey
            else '`{0}`={1}'.format(str(k).split('#')[0], R(v))
            for k, v in sorted(var.items())
        ]))
    return repr(var)


def Rlist(var, ignoreintkey=True):  # pylint: disable=invalid-name
    """Convert a dict into an R list"""
    assert isinstance(var, (list, tuple, set, dict))
    if isinstance(var, dict):
        return R(var, ignoreintkey)
    return 'as.list({})'.format(R(var, ignoreintkey))


def render(var, data=None):
    """
    Render a template variable, using the shared environment
    """
    if not isinstance(var, str):
        return var

    import inspect
    from pyppl.template import TemplateJinja2, TemplateLiquid

    frames = inspect.getouterframes(inspect.currentframe())
    data = data or {}
    for frame in frames:
        lvars = frame[0].f_locals
        if lvars.get('__engine') == 'liquid':
            evars = lvars.get('_liquid_context', {})
            if 'true' in evars:
                del evars['true']
            if 'false' in evars:
                del evars['false']
            if 'nil' in evars:
                del evars['nil']
            if '_liquid_liquid_filters' in evars:
                del evars['_liquid_liquid_filters']
            break
        if '_Context__self' in lvars:
            evars = dict(lvars['_Context__self'])
            break

    engine = evars.get('__engine')
    if not engine:
        raise RuntimeError(
            "I don't know which template engine to use to render {}...".format(
                var[:10]))

    engine = TemplateJinja2 if engine == 'jinja2' else TemplateLiquid
    return engine(var, **evars).render(data)


def glob1(*paths, first=True):
    """
    Return the paths matches the paths
    """
    assert len(paths) >= 2, "glob1 needs at least two parts for path."
    assert not (paths[-1] is True
                or paths[-1] is False), "`first` should be a keyword argument."
    paths = list(paths)
    path0 = paths.pop(0)
    pattern = paths.pop(-1)
    ret = list(Path(path0).joinpath(*paths).glob(pattern))
    if ret and first:
        return ret[0]  # Path object
    if not ret and first:
        return '__NoNeXiStFiLe__'
    return ret


def array_join(var, element_quote=None, all_quote=None, separator=' '):
    """Join an array"""
    var = (repr(str(element))
           if element_quote in ("'", 'single')
           else json.dumps(str(element))
           if element_quote in ('"', 'double')
           else element
           for element in var)
    var = separator.join(var)
    if all_quote in ("'", 'single'):
        return repr(var)
    if all_quote in ('"', 'double'):
        return json.dumps(var)
    return var


TEMPLATE_ENVS = dict(
    R=R,
    #Rvec     = R,            # will be deprecated!
    Rlist=Rlist,
    realpath=lambda var: Path(var).resolve().as_posix(),
    readlink=__import__('os').readlink,
    dirname=lambda var: Path(var).parent.as_posix(),
    # /a/b/c[1].txt => c.txt
    basename=basename,
    commonprefix=commonprefix,
    stem=filename,
    # /a/b/c.d.e.txt => c
    stem2=lambda var, orig=False, dot=1: filename(var, orig, dot),
    # /a/b/c.txt => .txt
    ext=lambda var: Path(var).suffix,
    glob1=glob1,
    # /a/b/c[1].txt => /a/b/c
    prefix=prefix,
    # /a/b/c.d.e.txt => /a/b/c
    prefix2=lambda var, orig=False, dot=1: prefix(var, orig, dot),
    # double quote string
    quote=lambda var: json.dumps(str(var)),
    squote=lambda var: repr(str(var)),
    json=json.dumps,
    read=read,
    readlines=readlines,
    render=render,
    array_join=array_join,
    rimport=rimport,
    bashimport=bashimport,
)

# aliases or reuses
# TEMPLATE_ENVS['readlink'] = TEMPLATE_ENVS['realpath']
TEMPLATE_ENVS['parent'] = TEMPLATE_ENVS['dirname']
TEMPLATE_ENVS['bn'] = TEMPLATE_ENVS['basename']
TEMPLATE_ENVS['filename'] = TEMPLATE_ENVS['stem']
TEMPLATE_ENVS['fn'] = TEMPLATE_ENVS['stem']
TEMPLATE_ENVS['filename2'] = TEMPLATE_ENVS['stem2']
TEMPLATE_ENVS['fn2'] = TEMPLATE_ENVS['stem2']
TEMPLATE_ENVS['ext2'] = lambda var: TEMPLATE_ENVS['ext'](var).lstrip('.')

DEFAULT_ENVS.update(TEMPLATE_ENVS)
