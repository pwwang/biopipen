"""A set of bioinformatics for PyPPL"""

import inspect
from pathlib import Path
import inflection
from diot import Diot
from pyppl import Proc
from varname import varname

# pylint: disable=invalid-name
# open to R (reticulate) to get the path of r util scripts
HERE = Path(__file__).resolve().parent

EXT_MAP = {
    'Rscript': 'R',
    'python': 'py',
    'python2': 'py',
    'python3': 'py',
}


def _findscript(script: str, callerdir: str) -> str:
    """Try to find the relative script"""
    if not script or not script.startswith('file:'):
        return script
    scriptfile = Path(script[5:])
    if scriptfile.is_absolute():
        return script
    scriptfile = callerdir.joinpath(scriptfile)
    return 'file:{}'.format(scriptfile)

# pylint: disable=redefined-builtin,invalid-name
def proc_factory(id=None, tag='notag', desc='No description.', **kwargs):
    """A factory to produce processes with default script,
    envs and report_template"""
    # in case if we have too long description
    id = id or varname()
    proc = Proc(id, tag=tag, desc=desc, **kwargs)
    lang = Path(proc.lang).name
    ext = '.' + EXT_MAP.get(lang, lang)
    callerfile = Path(
        inspect.getframeinfo(inspect.currentframe().f_back).filename)

    script = proc._script or 'file:scripts/{}/{}{}'.format(
        callerfile.stem, id, ext)
    proc.script = _findscript(script, callerfile.parent)
    report_template = proc.config.report_template or \
     'file:reports/{}/{}.md'.format(callerfile.stem, id)
    report_template = _findscript(report_template, callerfile.parent)
    if Path(report_template[5:]).is_file() or (report_template and
                                               report_template[:5] != 'file:'):
        # pylint: disable=assigning-non-slot
        proc.config.report_template = report_template
    return proc

def module_delegator(module: "ModuleType", name: str) -> Proc:
    """A delegator for submodule to do:
    >>> from submodule import pProcess
    >>> # or
    >>> import submodule
    >>> pProcess = submodule.pProcess
    So that we don't have to instantiate each `Proc` instance each time
    we import the submodule

    We will fit the name of the process with the `name`, and try to find
    the relative script and report file.

    Args:
        module (ModuleType): The submodule
        name (str): The name to be imported from the submodule

    Returns:
        Proc: Initialized `Proc` object
    """
    if name not in module._PROC_FACTORY:
        raise AttributeError(f"No such process {name} in "
                             f"module {module.__name__}")

    pfac = module._PROC_FACTORY[name].factory
    pfac_type = module._PROC_FACTORY[name].factype
    aliasof = module._PROC_FACTORY[name].aliasof

    if aliasof:
        module[name] = module[aliasof].copy(id=name)
        return module[name]

    if pfac_type == 'proc':
        module[name] = pfac.copy(id=name)
        return module[name]

    if pfac_type == 'class':
        proc_meta: dict = {key: getattr(pfac, key) for key in dir(pfac)
                           if not key.startswith('__')}
    else:
        proc_meta: dict = pfac()

    proc_meta.setdefault("id", name)
    proc = Proc(**proc_meta)

    if not proc.config.get('annotate') or not proc.config.annotate.sections:
        # pylint: disable=assigning-non-slot
        proc.config.annotate = pfac.__doc__

    # Automatically map the script and report file
    lang = Path(proc.lang).name
    ext = '.' + EXT_MAP.get(lang, lang)
    callerfile = Path(module.__file__)

    script = proc._script or f'file:scripts/{callerfile.stem}/{proc.id}{ext}'
    proc.script = _findscript(script, callerfile.parent)

    # Do the same for report file as well
    report_template = (
        proc.config.get('report_template') or
        'file:reports/{}/{}.md'.format(callerfile.stem, proc.id)
    )
    report_template = _findscript(report_template, callerfile.parent)

    if (Path(report_template[5:]).is_file() or
            (report_template and report_template[:5] != 'file:')):
         # pylint: disable=assigning-non-slot
        proc.config.report_template = report_template

    module[name] = proc
    return proc

def factory_type(name):
    """Detect the type of factory from a name
    The factories in submodule have to be written in such way

    Examples:

        >>> factory_type("pBcftoolsQuery")
        >>> # proc
        >>> factory_type("p_bcftools_query")
        >>> # function
        >>> factory_type("PBcftoolsQuery")
        >>> # class
    """
    if not name:
        return None
    if name.startswith('p_'):
        return 'function'
    if name.startswith('P') and name[1:2].isupper():
        return 'class'
    if name.startswith('p') and name[1:2].isupper():
        return 'proc'
    return None

def name_convert(name: str, to: str) -> str:
    """Convert one name of factory type to another"""
    frm = factory_type(name)
    if not frm:
        raise ValueError(f"Not any type of proc factory: {name}")

    assert (to in ('function', 'class', 'proc')), (
        "Expecting function, class or proc"
    )
    if frm == to:
        return name

    if frm == 'function':
        if to == 'class':
            # convert p_bcftools_query to PBcftoolsQuery
            return inflection.camelize(name, True)
        else: # to == proc
            # convert p_bcftools_query to pBcftoolsQuery
            return inflection.camelize(name, False)
    elif frm == 'class':
        if to == 'function':
            # convert PBcftoolsQuery to p_bcftools_query
            return inflection.underscore(name)
        else: # to == proc
            # convert PBcftoolsQuery to pBcftoolsQuery
            return 'p' + name[1:]
    else: # proc
        if to == 'class':
            # convert pBcftoolsQuery to PBcftoolsQuery
            return 'P' + name[1:]
        else: # to == function
            # convert pBcftoolsQuery to p_bcftools_query
            return inflection.underscore(name)

def module_postinit(module):
    """Calculate the real process ids for the submodule
    from the proxy function or class

    We should compose a meta directory with the proc ids as keys and the meta
    informtion for the process to be generated with those proc ids as values.

    The meta information should contain:
    - The original factory of this process
    - The type of the factory (function, class, or Proc)
    - The name of the original proc if I am an alias

    Aliases specified by modkit.aliases could be arbitrary.
    All of the following are legal:
    1. Aliasing a factory to a factory
        >>> modkit.aliases(p_bcftools_query=p_query)
        We only need to get the proc id from the original factory, and use the
        factory for the new proc.
    2. Aliasing a proc to a factory
        >>> modkit.aliases(pBcftoolsQuery=p_query)
        In such case, we need to remove the alias `pBcftoolsQuery` from
        `modkit.__modkit_meta__['aliases']`, and let modkit.delegate handle it.
        At meantime, we also need to clean __all__ to force recalulation.
        For the rest part, we have to do the same thing as listed in 1)
    3. Aliasing a factory to a proc
        >>> modkit.aliases(p_bcftools_query=pQuery)
        In this case, the factory should be generated for the new proc to
        make a copy of original proc.

    In either way, we need the ablity to detect the type
    (factory(function, class), Proc) from a name.
    """
    # factype, factory, aliasof
    factory = module._PROC_FACTORY = Diot()

    for attr in dir(module):
        typeof = factory_type(attr)
        if not typeof:
            # I am not a proc definition
            continue

        proc_id = name_convert(attr, 'proc')
        if proc_id in factory:
            raise ValueError(f"A factory for proc {proc_id} already exists")

        aliasof = module.__modkit_meta__['aliases'].get(attr)
        atypeof = factory_type(aliasof)

        if not atypeof:
            # We are not alias
            # Let's see what I am
            factory[proc_id] = dict(
                factype=typeof,
                factory=module[attr],
                aliasof=None
            )

            if typeof == 'proc':
                # I am not an alias, but a proc directly
                # Since I am defined directly, I will
                # not go through delegator
                factory[proc_id].factype = 'self'
                assert attr == module[attr].id, (
                    "Cannot pass in id for a process"
                )
        else: # I am an alias
            # we should remove this from __all__ so that it goes through
            # delegator, and factory can handle the proc generation
            del module.__modkit_meta__['aliases'][attr]
            module.__modkit_meta__['all'] = set()

            factory[proc_id] = dict(
                factype=None,
                factory=None,
                aliasof=name_convert(aliasof, 'proc')
            )
            if aliasof in module.__modkit_meta__['envs']:
                # We don't want to initialize the process here,
                # So let's check if the origin exists in envs
                # In other words, if it is defined in the submodule
                factory[proc_id].factype = atypeof
                factory[proc_id].factory = module[aliasof]
            elif atypeof != 'proc':
                # Otherwise, it has to be a process that is generated
                # by module delegator
                # We can't have undefined class or function factory deleagted
                # It has to be a process
                raise ValueError("Alias to an undefined "
                                 "function or class factory")
            elif module.__modkit_meta__['aliases'].get(aliasof):
                raise ValueError("Alias to an origin that is not "
                                 "directly defined.")
            else:
                origin_class = name_convert(aliasof, 'class')
                origin_func = name_convert(aliasof, 'function')
                if origin_class in module.__modkit_meta__['envs']:
                    factory[proc_id].factype = 'class'
                    factory[proc_id].factory = module[origin_class]
                elif origin_func in module.__modkit_meta__['envs']:
                    factory[proc_id].factype = 'function'
                    factory[proc_id].factory = module[origin_func]
                else:
                    raise ValueError("Alias to undefined origin.")

    module.__modkit_meta__['delegate'] = module_delegator
