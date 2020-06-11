"""Utilities for bioprocs script"""
import re
import sys
from pathlib import Path
from tempfile import gettempdir
from contextlib import contextmanager
import inflection
from modkit import Module as ModkitModule
from colorama import Back
from pyparam.help import Helps, HelpAssembler
from diot import Diot, OrderedDiot
from pyppl import Channel, PyPPL, Proc
from pyppl._proc import IN_FILESTYPE
from liquid.stream import LiquidStream
import bioprocs


def _split(string, delimit, trim=True):
    return LiquidStream.from_string(string).split(delimit, trim=trim)


def substr_replace(string, starts, lengths, replace):
    """Replace substrings"""
    if not isinstance(starts, (tuple, list)):
        starts = [starts]
    if not isinstance(lengths, (tuple, list)):
        lengths = [lengths]
    assert len(starts) == len(lengths)
    if not isinstance(replace, (tuple, list)):
        replace = [replace] * len(starts)

    delta = 0
    for i, start in enumerate(starts):
        # adjust starts
        string = string[:start +
                        delta] + replace[i] + string[start + lengths[i] +
                                                     delta:]
        delta += len(replace[i]) - lengths[i]
    return string


def highlight(origin, query, incase=True, hicolor=Back.LIGHTRED_EX):
    """Highlight string with query string"""
    # get all occurrences of q
    if incase:
        occurs = [
            m.start() for m in re.finditer(query.lower(), origin.lower())
        ]
    else:
        occurs = [m.start() for m in re.finditer(query, origin)]
    lengths = [len(query)] * len(occurs)
    return substr_replace(origin, occurs, lengths, [
        '{}{}{}'.format(hicolor, origin[occur:occur + length], Back.RESET)
        for occur, length in zip(occurs, lengths)
    ])


def highlight_multi(line, queries):
    """Highlight a string with multiple queries"""
    for query in queries:
        line = highlight(line, query)
    return line


def subtract_dict(bigger, smaller, prefix=''):
    """Subtract a dict from another"""
    ret = bigger.copy()
    for key, val in smaller.items():
        if key not in ret:
            continue
        if isinstance(ret[key], dict) and isinstance(val,
                                                     dict) and ret[key] != val:
            ret[key] = subtract_dict(ret[key], val, prefix + '  ')
        elif ret[key] == val:
            del ret[key]
    return ret


class Module:
    """A module of bioprocs"""
    @staticmethod
    def modules():
        """Get all modules"""
        return [
            module.stem
            for module in Path(bioprocs.__file__).parent.glob('*.py')
            if not module.stem.startswith('_')
        ]

    def __init__(self, name):
        self.name = name
        try:
            self.module = getattr(__import__('bioprocs', fromlist=[name]),
                                  name)
        except AttributeError as ex:
            if "has no attribute '%s'" % name in str(ex):
                raise AttributeError('No such module: %s' % name) from None
            raise
        self.desc = (self.module.__doc__.strip()
                     if self.module.__doc__
                     else '[ Not documented ]')
        self._procs = {}

    def procs(self):
        """Get the processes of the module"""
        if self._procs:
            return self._procs

        if isinstance(self.module, ModkitModule):
            for proc_id in self.module._PROC_FACTORY:
                proc = getattr(self.module, proc_id)
                self._procs[proc_id] = Process(
                    proc, self.module,
                    proc.config.annotate,
                    self.module._PROC_FACTORY[proc_id].aliasof
                )

            return self._procs

        for attr in dir(self.module):
            proc = getattr(self.module, attr)
            if not isinstance(proc, Proc):
                continue
            self._procs[attr] = Process(proc, self.module,
                                        proc.config.annotate)

        return self._procs

    def to_helps(self, helpsec):
        """Send me to a Helps section"""
        helpsec.prefix = ''
        helpsec.add((self.name, '', self.desc))

    @contextmanager
    def load_procs(self, index=None):
        """Loading procs with indicators"""
        info = 'Collecting module%s: %s ...' % (index and
                                                (' ' + index) or '', self.name)
        print('\r' + info, end='')
        yield self.procs()
        print(end='\r')
        print(' ' * len(info), end='')
        if index:
            idx, total = index.split('/')
            if idx == total:
                print(end='\r')

    def to_helps_as_sec(self, helps):
        """Add me to helps as a section"""
        helps.add(self.name + ': ' + self.desc, sectype='option', prefix='')

    def to_helps_with_procs(self, helps, nmods):
        """Add me to helps with procs information"""
        self.to_helps_as_sec(helps)
        with self.load_procs("%s/%s" % (len(helps), nmods)) as procs:
            for proc in procs.values():
                proc.to_helps(helps.select(self.name + ': '))

    def add_to_completions(self, comp):
        """Add me to completions"""
        comp.add_command(self.name, self.desc)
        with self.load_procs() as procs:
            for pname, proc in procs.items():
                comp.add_command(self.name + '.' + pname, proc.desc)
                proc.add_to_completions(comp.command(self.name + '.' + pname))


class Process:
    """A bioprocs process"""
    def __init__(self, proc, module, doc, aliasof=None):
        self.name = proc.id
        self.module = module
        self.proc = proc
        self.desc = self.proc.desc
        self.doc = doc
        self.aliasof = aliasof
        self._helps = Helps()
        self._parsed = {}

    def to_helps(self, helpsec):
        """Send me to a Helps section"""
        helpsec.prefix = ''
        if self.aliasof:
            helpsec.add((self.name, '', 'Alias of `%s`' % self.aliasof))
        elif not self.proc.origin or self.proc.origin == self.name:
            helpsec.add((self.name, '', self.desc))
        else:
            helpsec.add((self.name, '', 'Alias of `%s`' % self.proc.origin))

    def requirements(self):
        """Get requirements of process"""
        return self.doc.section('requires') or '[Not documented.]'

    def add_to_completions(self, comp):
        """Add me to completions"""
        for inname in self.doc.input or []:
            comp.add_option('-i.' + inname, self.doc.input[inname]['desc'])
        for outname in self.doc.output or []:
            comp.add_option('-o.' + outname, self.doc.output[outname]['desc'])
        for key in self.doc.args or []:
            comp.add_option('--args.' + key, self.doc.args[key]['desc'])
        # add -config.
        comp.add_option('--config.',
                        'Plugin configrations, such as --config.export_dir')

    def helps(self):
        """Construct help page using doc"""
        if self._helps:
            return self._helps

        self._helps.add(
            'Name', '%s.%s %s[lang = %s]' %
            (Path(self.module.__file__).stem, self.name,
             '(%s)' % self.proc.origin if self.proc.origin
             and self.name != self.proc.origin else '', self.proc.lang))
        self._helps.add('Description', self.doc.description or self.desc)

        # input
        self._helps.add('Input options (Use \':list\' for multi-jobs)',
                        sectype='option')
        for inname in (self.doc.input or []):
            self._helps.select('Input options').add(
                ('-i.' + inname, '<%s>' % self.doc.input[inname]['type'],
                 self.doc.input[inname]['desc']))

        # output
        self._helps.add('Output options (\'--config.export_dir\' '
                        'implied if path specified)',
                        sectype='option')

        for outname in (self.doc.output or []):
            self._helps.select('Output options').add((
                '-o.' + outname,
                '<%s>' % self.doc.output[outname]['type'],
                self.doc.output[outname]['desc'] + \
                'Default: ' + repr(self.doc.output[outname]['default'])))

        # args
        self._helps.add('Process arguments', sectype='option')
        for key in (self.doc.args or {}):
            self._helps.select('Process arguments').add((
                '--args.' + key,
                '<%s>' % self.doc.args[key]['type'],
                self.doc.args[key]['desc'] + \
                'Default: ' + str(self.doc.args[key]['default'])))

        self._helps.add('Other options', sectype='option')

        # pipeline configurations
        self._helps.select('Other options').add(
            ('--config.<subconf>[.subconf]', '',
             'Plugin configrations, such as --config.export_dir'))

        # help
        self._helps.select('Other options').add_param(
            bioprocs.params[bioprocs.params._hopts[0]],
            bioprocs.params._hopts,
            ishelp=True)

        self._helps.add('Requirements', self.requirements(), sectype='plain')
        return self._helps

    @staticmethod
    def _update_args(args):
        # replace ',' with '.' in key
        # a bug of python-box, copy lost metadata
        #ret = args.copy() # try to keep the type
        ret = Diot()
        for key, val in args.items():
            ret[key.replace(',', '.')] = Process._update_args(val) \
                if isinstance(val, dict) else val
        return ret

    def print_helps(self, error=None, halt=True):
        """Print helps"""
        assembler = HelpAssembler()
        error = error or []
        if isinstance(error, str):
            error = [error]
        for err in error:
            print(assembler.error(err))
        print('\n'.join(assembler.assemble(self.helps())), end='')
        if halt:
            sys.exit(1)

    def _logger(self, # pylint: disable=too-many-arguments
                msg,
                hiword='',
                hicolor=Back.GREEN,
                end='\n',
                prefix=True):
        if prefix:
            modname = Path(self.module.__file__).stem
            prefix_str = f'[{modname}.{self.name}] '
            prefix_str = highlight(prefix_str,
                                   prefix_str[1:-2],
                                   incase=False,
                                   hicolor=Back.MAGENTA)
        else:
            prefix_str = ''
        if not hiword:
            print(f'{prefix_str}{msg}', end=end)
        else:
            himsg = highlight(msg, hiword, hicolor=hicolor)
            print(f'{prefix_str}{himsg}', end=end)

    def run(self, opts): # pylint: disable=too-many-branches
        """Construct a pipeline with the process and run it"""
        if (any(opts.get(h) for h in bioprocs.params._hopts)
                or all(key in bioprocs.params._hopts for key in opts)):
            self.print_helps()

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

        if opts.get('o') is not None:
            if not isinstance(opts.o, dict):
                self.print_helps(error='Malformat output specification.')

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

        args = opts.pop('args', {})
        if not isinstance(args, dict):
            self.print_helps(error='Malformat args specification.')
        self.proc.args.update(Process._update_args(args))

        self.proc.config.update(opts.pop('config', {}))
        # config = {
        # 	'logger'  : {'file': None },
        # 	'ppldir': Path(gettempdir()) / 'bioprocs.workdir'
        # }
        # config.update(Process._update_args(opts.get('config', {})))

        for key, val in opts.items():
            if key in bioprocs.params._hopts + ['i', 'o']:
                continue
            setattr(self.proc, key, val)

        PyPPL(logger_file=False,
              ppldir=Path(gettempdir()) / 'bioprocs.workdir').start(
                  self.proc
              ).run()


class Pipeline:
    """Assembled pipeline"""
    @staticmethod
    def pipelines():
        """Get all available pipelines"""
        return [
            pplfile.stem[9:]
            for pplfile in (Path(bioprocs.__file__).parent / 'console' /
                            'pipeline').glob('bioprocs_*.py')
            if not pplfile.stem.startswith('_')
        ]

    def __init__(self, name):
        self.name = name
        self.module = __import__('bioprocs.console.pipeline',
                                 fromlist=['bioprocs_' + name])
        self.module = getattr(self.module, 'bioprocs_' + name)
        self.desc = (self.module.__doc__.strip()
                     if self.module.__doc__
                     else '[ Not documented ]')

    def run(self):
        """Run the pipeline"""
        prog = 'bioprocs ' + self.name
        sys.argv = [prog] + sys.argv[2:]
        bioprocs.params._prog = prog
        bioprocs.params._assembler.progname = prog

        self.module.main()

    def add_to_completions(self, comp):
        """Add me to completions"""
        comp.add_command(self.name, self.desc)
        self.module.params._add_to_completions(comp.command(self.name),
                                               withtype=False,
                                               alias=True)
