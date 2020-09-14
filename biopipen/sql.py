"""
Some sql database utilities
"""
from diot import Diot
from . import opts, proc_factory
from .utils import sql

# pylint: disable=invalid-name

pCreateTable = proc_factory(
    desc='Create tables in the database.',
    config=Diot(annotate="""
    @name:
        pCreateTable
    @description:
        Create tables in the database
    @input:
        `dsn`: The dsn to connect to the database
            - currently support `sqlite:file=...`
        `schema:file`: The schema file
            - could be a pure schema file:
            ```
            Field	Type	Statement
            ID	INT	PRIMARY KEY
            ...
            ```
            - or a data file with header
    @output:
        `dsn`: The dsn
    @args:
        `intype`: The input file schema file or a data file. Default: `schema`
        `drop`:  Force creating the table (drop the pre-existing table)
        `delimit`:The delimit of input file. Default: `\\t`
    """),
    input="dsn, schema:file",
    output="dsn:var:{{i.dsn}}",
    lang=opts.python,
    args=Diot(
        intype='schema',
        drop=False,
        delimit="\t",
    ),
    envs=Diot(
        dsnparse=sql.dsnparse,
        schemaparse=sql.schemaparse,
    )
)

pImportData = proc_factory(
    desc='Create tables and import the data',
    config=Diot(annotate="""
    @name:
        pImportData
    @description:
        Create tables and import the data
    @input:
        `dsn`: The dsn to connect to the database
            - currently support `sqlite:file=...`
        `datafile:file`: The schema file
            - must have header
    @output:
        `dsn`: The dsn
    @args:
        `delimit`:The delimit of input file. Default: `\\t`
    """),
    input="dsn, datafile:file",
    output="dsn:var:{{i.dsn}}",
    lang=opts.python,
    args=Diot(
        delimit="\t",
        drop=False,
    ),
    envs=Diot(
        dsnparse=sql.dsnparse,
        schemaparse=sql.schemaparse,
        tablename=lambda fn: __import__('os').path.basename(
            fn).split('.')[0]
    )
)

pUpdateTable = proc_factory(
    desc='Update table using sql.',
    config=Diot(annotate="""
    @name:
        pUpdateTable
    @description:
        Update table using sql.
    @input:
        `dsn`: The dsn to connect to the database
            - currently support `sqlite:file=...`
    @output:
        `dsn`: The dsn
    @args:
        `sql`: The sql to update the table (list)
    """),
    input='dsn',
    output='dsn:var:{{i.dsn}}',
    lang=opts.python,
    args=Diot(sql=[]),
    envs=Diot(dsnparse=sql.dsnparse)
)

pSelectTable = proc_factory(
    desc='Select data from table and dump it.',
    config=Diot(annotate="""
    @name:
        pSelectTable
    @description:
        Select data from table and dump it.
    @input:
        `dsn`: The dsn to connect to the database
            - currently support `sqlite:file=...`
    @output:
        `outfile:file`: The dumped file
    @args:
        `sql`: The sql to select data from the table (list)
    """),
    input='dsn',
    output='outfile:file:sqlSelectTable-{{i.dsn, args.sql | :_1+_2 \
                                                          | md5 \
                                                          | .hexdigest: \
                                                          | [:16]}}.dumped.txt',
    lang=opts.python,
    args=Diot(sql=[]),
    envs=Diot(dsnparse=sql.dsnparse, md5=__import__('hashlib').md5)
)
