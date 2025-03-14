from biopipen.utils.misc import run_command, dict_to_cli_args
from biopipen.scripts.web.gcloud_common import is_logged_in, is_valid_gs_file_url

url: str = {{in.url | repr}}  # pyright: ignore  # noqa: E999
outfile = {{out.outfile | repr}}  # pyright: ignore
gcloud: str = {{envs.gcloud | repr}}  # pyright: ignore
args: dict = {{envs.args | repr}}  # pyright: ignore

if not is_valid_gs_file_url(url):
    raise Exception(
        f"Invalid Google Cloud Storage URL for a file: {url}. "
        "URL should be in the format gs://bucket/path/to/file"
    )

if not is_logged_in(gcloud):
    raise Exception(
        "You need to be logged in to gcloud to download files. "
        "Please run `gcloud auth login` first."
    )

args[""] = [gcloud, "storage", "cp", url, outfile]

run_command(dict_to_cli_args(args, dashify=True), fg=True)
