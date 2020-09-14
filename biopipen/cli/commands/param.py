"""List params or get the value of a param."""
import sys
from pyparam import Params, POSITIONAL
from rich import print as rich_print, box
from rich.table import Table
from rich.highlighter import RegexHighlighter
from ..utils import logger

params = Params(desc=__doc__, help_on_void=False)
params.add_param('g,get', default=False, desc='Get the value of a parameter.')
params.add_param(POSITIONAL, required=True, type=str,
                 desc=['List params with names containing given string.',
                       'When -g or --get is given, this has to be the '
                       'exact param name'])

class QueryHighlighter(RegexHighlighter):
    """Apply style for the query word."""

    def __init__(self, kw: str):
        super().__init__()
        self.highlights = [rf"(?P<code>{kw})"]

def main(opts):
    from ..._params import params, opts as pmopts
    if opts.get:
        if opts[POSITIONAL] not in pmopts:
            logger.error('No such parameter: %r', opts[POSITIONAL])
            sys.exit(1)

        rich_print(pmopts[opts[POSITIONAL]])
    else:
        for qkey in vars(pmopts):
            if not opts[POSITIONAL] in qkey:
                continue

            param = params.get_param(qkey)
            highlighter = QueryHighlighter(opts[POSITIONAL])
            table = Table(show_header=False, box=box.SQUARE, expand=True)
            table.add_column(style="green", width=10, ratio=1)
            table.add_column(style="italic", overflow='fold', ratio=9)
            table.add_row('Name', highlighter(qkey), style="bold")
            table.add_row('Type', param.typestr())
            table.add_row('Value', repr(pmopts[qkey]))

            rich_print(table)



if __name__ == '__main__':
    parsed = params.parse()
    main(parsed)
