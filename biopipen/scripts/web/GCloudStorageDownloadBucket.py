import concurrent.futures
from pathlib import Path

from biopipen.utils.misc import run_command, dict_to_cli_args, logger
from biopipen.scripts.web.gcloud_common import (
    is_logged_in,
    is_valid_gs_bucket_url,
    get_file_path,
)

url: str = {{in.url | quote}}  # pyright: ignore  # noqa: E999
outdir = Path({{out.outdir | repr}})  # pyright: ignore
gcloud: str = {{envs.gcloud | quote}}  # pyright: ignore
keep_structure = {{envs.keep_structure | repr}}  # pyright: ignore
ncores: int = {{envs.ncores | repr}}  # pyright: ignore
args: dict = {{envs.args | repr}}  # pyright: ignore

if not is_valid_gs_bucket_url(url):
    raise Exception(
        f"Invalid Google Cloud Storage URL for a bucket: {url}. "
        "URL should be in the format gs://bucket"
    )

if not is_logged_in(gcloud):
    raise Exception(
        "You need to be logged in to gcloud to download files. "
        "Please run `gcloud auth login` first."
    )


def create_folders(folder_lines):
    for folder_line in folder_lines:
        folder_path = get_file_path(folder_line)
        folder = outdir / folder_path
        folder.mkdir(parents=True, exist_ok=True)


def download_file(i: int, line: str, total: int):
    path = get_file_path(line)

    if total <= 50:
        logger.info(f"Downloading {path}")
    elif 50 < total <= 500:
        if i % 10 == 0:
            logger.info(f"Downloading {i}/{total} ...")
    else:
        if i % 100 == 0:
            logger.info(f"Downloading {i}/{total} ...")

    if keep_structure:
        target = (outdir / path)
    else:
        name = Path(path).name
        target = outdir / name
        if target.exists():
            new_name = f"g{i}-{name}"
            logger.warning(f"{name} already exists. Renaming to {new_name}.")
            target = outdir / new_name

    gs_args = args.copy()
    gs_args[""] = [gcloud, "storage", "cp", line, target]
    run_command(dict_to_cli_args(gs_args, dashify=True), fg=True)


def download_bucket():
    out = run_command([gcloud, "storage", "ls", "--recursive", url], stdout="RETURN")
    # remove empty lines and skip the root
    out = list(filter(None, out.splitlines()[1:]))  # type: ignore
    if keep_structure:
        # create folders first
        logger.info(f"Creating folders to keep structure.")
        folder_lines = [line[:-2] for line in out if line.endswith("/:")]
        create_folders(folder_lines)

    out = [line for line in out if not line.endswith("/:")]
    length = len(out)
    with concurrent.futures.ProcessPoolExecutor(max_workers=ncores) as executor:
        executor.map(download_file, range(length), out, [length] * length)


if __name__ == "__main__":
    download_bucket()
