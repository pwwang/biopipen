"""Download/Get data from Websites instead of APIs"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name
# Web utils
pDownloadForm = proc_factory(
    desc="Download results by submitting to a form",
    config=Diot(annotate="""
    @name:
        pDownloadForm
    @description:
        Download results by submitting a form, supporting pagination.
    @input:
        `url`   : the URL contains the form
        `data`  : the data used to fill the form (JSON string or transformed from dict by json.dumps).
        `submit`: the submit button to submit the form (use Xpath).
        `next`  : the button for next page (use Xpath)
    @output:
        `outdir:file`: The directory saves the results
    @args:
        `interval`: seconds to wait between fetching each page. Default: 1
    @requires:
        [`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
        [`Phantomjs`](http://phantomjs.org/)
    """),
    input="url, data, submit, next",
    output="outdir:dir:{{i.url | bn | lambda x: x if x else 'outdir'}}",
    lang=opts.python,
    args=Diot(interval=1)
)

pDownloadGet = proc_factory(
    desc="Download from URLs",
    config=Diot(annotate="""
    @name:
        pDownloadGet
    @description:
        Download results by urls.
    @input:
        `url`: the URLs to download
    @output:
        `outfile:file`: The output file
    """),
    input="url",
    output="outfile:file:{{i.url | bn | ?!:'outfile' \
                                      | $.replace: '?', '__Q__' \
                                      | .replace: '&', '__N__'}}",
    lang=opts.python,
)

pDownloadPost = proc_factory(
    desc="Download from URLs",
    config=Diot(annotate="""
    @name:
        pDownloadPost
    @description:
        Download results by POST.
    @input:
        `url` : the URLs to download
        `data`: the POST data.
    @output:
        `outfile:file`: The output file
    """),
    input="url, data",
    output="outfile:file:{{i.url | bn | ?!:'outfile'}}",
    lang=opts.python,
)
