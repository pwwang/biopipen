"""Get data from the web"""
from ..core.proc import Proc
from ..core.config import config


class Download(Proc):
    """Download data from URLs

    Input:
        url: The URL to download data from

    Output:
        outfile: The file downloaded

    Envs:
        tool (choice): Which tool to use to download the data
            - wget: Use wget
            - aria2c: Use aria2c
            - urllib: Use python's urllib
            - aria: Alias for aria2c
        wget: Path to wget
        aria2c: Path to aria2c
        args: The arguments to pass to the tool
        ncores: The number of cores to use

    Requires:
        wget: Only required when envs.tool == "wget"
            - check: {{proc.envs.wget}} --version
        aria2c: Only required when envs.tool == "aria2c"
            - check: {{proc.envs.aria2c}} --version
    """
    input = "url"
    output = (
        # Need to replace http:// and https:// to avoid cloudpathlib.AnyPath to get
        # the basename for something like "https://example.com/data/?file=datafile.txt"
        # as data, but "?file=datafile.txt"
        "outfile:file:"
        """{{in.url
            | replace: 'http://', ''
            | replace: 'https://', ''
            | basename
            | url_decode
            | slugify: separator='.', lowercase=False, regex_pattern='[^-a-zA-Z0-9_]+'
        }}"""
    )
    lang = config.lang.python
    envs = {
        "tool": "wget",  # or aria2c, python
        "wget": config.exe.wget,
        "aria2c": config.exe.aria2c,
        "args": {},
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/web/Download.py"


class DownloadList(Proc):
    """Download data from URLs in a file.

    This does not work by iterating over the URLs in the file. The whole file is
    passed to `wget` or `aria2c` at once.

    Input:
        urlfile: The file containing the URLs to download data from

    Output:
        outdir: The directory containing the downloaded files

    Envs:
        tool (choice): Which tool to use to download the data
            - wget: Use wget
            - aria2c: Use aria2c
            - urllib: Use python's urllib
            - aria: Alias for aria2c
        wget: Path to wget
        aria2c: Path to aria2c
        args: The arguments to pass to the tool
        ncores: The number of cores to use

    Requires:
        wget: Only required when envs.tool == "wget"
            - check: {{proc.envs.wget}} --version
        aria2c: Only required when envs.tool == "aria2c"
            - check: {{proc.envs.aria2c}} --version
    """
    input = "urlfile:file"
    output = "outdir:dir:{{in.urlfile | stem}}.downloaded"
    lang = config.lang.python
    envs = {
        "tool": "wget",  # or aria2c
        "wget": config.exe.wget,
        "aria2c": config.exe.aria2c,
        "args": {},
        "ncores": config.misc.ncores,
    }
    script = "file://../scripts/web/DownloadList.py"


class GCloudStorageDownloadFile(Proc):
    """Download file from Google Cloud Storage

    Before using this, make sure you have the `gcloud` tool installed and
    logged in with the appropriate credentials using `gcloud auth login`.

    Also make sure you have [`google-crc32c`](https://pypi.org/project/google-crc32c/)
    installed to verify the integrity of the downloaded files.

    Input:
        url: The URL to download data from.
            It should be in the format gs://bucket/path/to/file

    Output:
        outfile: The file downloaded

    Envs:
        gcloud: Path to gcloud
        args (ns): Other arguments to pass to the `gcloud storage cp` command
            - do_not_decompress (flag): Do not decompress the file.
            - <more>: More arguments to pass to the `gcloud storage cp` command
                See `gcloud storage cp --help` for more information
    """
    input = "url:var"
    output = "outfile:file:{{in.url | replace: 'gs://', '/' | basename}}"
    lang = config.lang.python
    envs = {
        "gcloud": config.exe.gcloud,
        "args": {"do_not_decompress": True},
    }
    script = "file://../scripts/web/GCloudStorageDownloadFile.py"


class GCloudStorageDownloadBucket(Proc):
    """Download all files from a Google Cloud Storage bucket

    Before using this, make sure you have the `gcloud` tool installed and
    logged in with the appropriate credentials using `gcloud auth login`.

    Note that this will not use the `--recursive` flag of `gcloud storage cp`.
    The files will be listed and downloaded one by one so that they can be parallelized.

    Also make sure you have [`google-crc32c`](https://pypi.org/project/google-crc32c/)
    installed to verify the integrity of the downloaded files.

    Input:
        url: The URL to download data from.
            It should be in the format gs://bucket

    Output:
        outdir: The directory containing the downloaded files

    Envs:
        gcloud: Path to gcloud
        keep_structure (flag): Keep the directory structure of the bucket
        ncores (type=int): The number of cores to use to download the files in parallel
        args (ns): Other arguments to pass to the `gcloud storage cp` command
            - do_not_decompress (flag): Do not decompress the file.
            - <more>: More arguments to pass to the `gcloud storage cp` command
                See `gcloud storage cp --help` for more information
    """
    input = "url:var"
    output = "outdir:dir:{{in.url | replace: 'gs://', ''}}"
    lang = config.lang.python
    envs = {
        "gcloud": config.exe.gcloud,
        "keep_structure": True,
        "ncores": config.misc.ncores,
        "args": {"do_not_decompress": True},
    }
    script = "file://../scripts/web/GCloudStorageDownloadBucket.py"
