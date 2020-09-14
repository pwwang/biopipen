"""Proxy for biopipen-xxx"""
from pathlib import Path
from functools import lru_cache
from pyparam import Params
from diot import OrderedDiot
from pyppl import Proc, Channel, PyPPL
from pyppl._proc import IN_FILESTYPE
from modkit import modkit, Module as ModkitModule
from tempfile import gettempdir
from .utils import import_from_path, ROOT_MODULE

MODULE_DIR = Path(__file__).parent.parent.resolve()

PROCATTR_DESC = {
    'cache': 'Whether we should use the cached process.',
    'dirsig': 'Whether check the dir signature, slow for large directory.',
    'envs': 'The environment variables for template.',
    'errhow': 'What to do if errors happens.',
    'errntry': 'How many times to retry when errhow is `retry`.',
    'forks': 'How many jobs to run simutaneously.',
    'nthread': 'Number of threads to use to build jobs.',
    'ppldir': 'Pipeline directory',
    'runner': 'The running profile, or runner.',
    'tag': 'The tag of the process.'
}

class Module:
    """A module of processes"""

    def __init__(self, path: Path):
        self.name = path.stem
        # load the module from path
        module = import_from_path(path)
        self.doc = (module.__doc__.strip()
                    if module.__doc__
                    else '[ Undocumented ]')
        self.procs = Process.collections(module)
        self.params = Params(prog=f"{ROOT_MODULE}-{self.name}", desc=self.doc,
                             help_on_void=False)
        self.load_params()

    def load_params(self):
        """Load the parameters for each process"""
        for procname, proc in self.procs.items():
            command = self.params.add_command(
                procname,
                desc=proc.doc,
                help_on_void=True,
                help_callback=proc._help_callback(f"{ROOT_MODULE}-{self.name}")
            )
            proc.to_params(command)

    def main(self):
        """Main entry point"""
        parsed = self.params.parse()
        self.procs[parsed.__command__].run(parsed[parsed.__command__])

class Process:
    """A process"""

    @classmethod
    @lru_cache()
    def collections(cls, module):
        """The collections of processes of given module."""
        processes = {}
        if isinstance(module, ModkitModule):
            for proc_id in module._PROC_FACTORY:
                proc = getattr(module, proc_id)
                processes[proc_id] = cls(
                    proc,
                    module._PROC_FACTORY[proc_id].aliasof
                )

            return processes

        for attr in dir(module):
            proc = getattr(module, attr)
            if not isinstance(proc, Proc):
                continue
            processes[attr] = cls(proc)

        return processes

    def __init__(self, proc, aliasof=None):
        self.proc = proc
        self.doc = proc.config.annotate
        self.aliasof = aliasof

    def to_params(self, params):
        """Load to params"""
        self._desc_to_params(params)
        self._input_to_params(params)
        self._output_to_params(params)
        self._args_to_params(params)
        self._procattr_to_params(params)
        self._plugincfg_to_params(params)

    def _help_callback(self, root_prog):
        def help_callback(helps):
            langline = (
                f"{root_prog} {self.proc.id} "
                f"%s[lang = {self.proc.lang}]" % (
                    "" if not self.proc.origin else f"({self.proc.origin}) "
                )
            )
            helps.DESCRIPTION.insert(0, langline)
            helps.DESCRIPTION.insert(1, "")
            if self.aliasof:
                helps.DESCRIPTION.append("")
                if isinstance(self.doc.description, list):
                    helps.DESCRIPTION.extend(
                        self.doc.description or [self.proc.desc]
                    )
                else:
                    helps.DESCRIPTION.append(
                        self.doc.description or self.proc.desc
                    )

            # helps['PLUGIN CONFIGIRATIONS'] = HelpSectionOption([
            #     (['--config{.subconf}{.subconf}'],
            #      ['Configurations for plugins. For example: '
            #       '`--config.export_dir`'])
            # ])
        return help_callback

    def _desc_to_params(self, params):
        """Modify the description of the process in help page"""
        params.desc = [f'Alias of {self.aliasof}.'] if self.aliasof else []

        if not self.aliasof:
            if isinstance(self.doc.description, list):
                params.desc.extend(self.doc.description or [self.proc.desc])
            else:
                params.desc.append(self.doc.description or self.proc.desc)

    def _input_to_params(self, params):
        params.add_param('i', type="ns", desc="Input options", show=False)
        for inname in (self.doc.input or []):
            params.add_param(
                f'i.{inname}',
                desc=self.doc.input[inname]['desc'],
                argname_shorten=False,
                group='INPUT OPTIONS',
                type=("list:path"
                      if self.doc.input[inname]['type'] in ('files', 'paths')
                      else "list:dir"
                      if self.doc.input[inname]['type'] == 'dirs'
                      else "list:path"
                      if self.doc.input[inname]['type'] in ('file', 'path')
                      else "list:dir"
                      if self.doc.input[inname]['type'] == 'dir'
                      else "list:str")
            )

    def _output_to_params(self, params):
        params.add_param('o', type="ns", desc="Output options", show=False)
        for outname in (self.doc.output or []):
            params.add_param(
                f'o.{outname}',
                desc=self.doc.output[outname]['desc'],
                argname_shorten=False,
                group=('OUTPUT OPTIONS (--config.export_dir '
                       'implied if path spedified)'),
                type="str",
                default=str(self.doc.output[outname]['default'])
            )

    def _args_to_params(self, params):
        params.add_param('args', type="ns", desc="Output options", show=False)
        for aname in (self.doc.args or []):
            params.add_param(
                f'args.{aname}',
                desc=self.doc.args[aname]['desc'],
                argname_shorten=False,
                group=('ARGS OPTIONS'),
                type=("ns" if isinstance(self.doc.args[aname]['default'], dict)
                      else None),
                default=self.doc.args[aname]['default']
            )

    def _procattr_to_params(self, params):
        for attr in dir(self.proc):
            if attr not in PROCATTR_DESC:
                continue
            default = (
                getattr(self.proc, attr) if attr != 'runner'
                else self.proc._runner if isinstance(self.proc._runner, str)
                else self.proc._runner.runner
            )
            params.add_param(
                attr,
                group='OTHER PROCESS ATTRIBUTES',
                desc=PROCATTR_DESC[attr],
                type=("ns" if isinstance(default, dict) else None),
                default=default
            )

    def _plugincfg_to_params(self, params):
        params.add_param('config', type="ns", desc="Plugin configurations",
                         show=False)
        for cname in self.proc.config:
            if cname == 'annotate':
                continue
            default = self.proc.config[cname]
            params.add_param(
                f'config.{cname}',
                argname_shorten=False,
                group=('PLUGIN CONFIGIRATIONS'),
                type=("ns" if isinstance(default, dict) else None),
                default=default
            )

    def run(self, opts):
        """Run the process as a command line tool"""
        opts = self.params.parse()

        indata = OrderedDiot()
        for inkey in self.doc.input:
            # We should allow some input options to be empty
            # to be smart on "files"
            intype = self.doc.input[inkey]['type']
            if (intype in IN_FILESTYPE
                    and inkey in opts.i
                    and not isinstance(opts.i[inkey], list)):
                indata[inkey + ':' + intype] = Channel.create([[opts.i[inkey]]
                                                               ])
            else:
                indata[inkey + ':' + intype] = Channel.create(opts.i.get(inkey))
        self.proc.input = indata

        outdata = OrderedDiot()
        for outkey in self.doc.output:
            outype, outdeft = self.doc.output[outkey][
                'type'], self.doc.output[outkey]['default']
            if not opts.o.get(outkey):
                outdata[outkey + ':' + outype] = outdeft
                continue
            if outype not in ('file', 'dir'):
                outdata[outkey + ':' + outype] = opts.o[outkey]
                continue
            # try to extract exdir from output
            if '/' in opts.o[outkey]:
                out = Path(opts.o[outkey])
                if (self.proc.config.export_dir
                        and out.parent != self.proc.config.export_dir):
                    raise ValueError('Cannot have output files/dirs '
                                        'with different parents as exdir.')
                self.proc.config.export_dir = str(out.parent)
                outdata[outkey + ':' + outype] = out.name
            else:
                outdata[outkey + ':' + outype] = opts.o[outkey]
        self.proc.output = outdata
        self.proc.args.update(vars(opts.args))
        self.proc.config.update(vars(opts.config))
        # config = {
        # 	'logger'  : {'file': None },
        # 	'ppldir': Path(gettempdir()) / 'bioprocs.workdir'
        # }
        # config.update(Process._update_args(opts.get('config', {})))

        for key, val in opts.items():
            if key in ('i', 'o', 'args', 'config'):
                continue
            setattr(self.proc, key, val)

        PyPPL(logger_file=False,
              ppldir=Path(gettempdir()) / 'bioprocs.workdir').start(
                  self.proc
              ).run()

@modkit.delegate
def _console_script_delegator(module, name): # pylint: disable=unused-argument
    """The delegator for console script

    For example:
    biopipen-vcf = biopipen.cli.modules.vcf
    """
    return Module(MODULE_DIR.joinpath(name + '.py')).main
