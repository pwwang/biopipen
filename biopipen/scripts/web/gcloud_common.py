"""Provides common functions for interacting with Google Cloud Storage."""
from biopipen.utils.misc import run_command


def is_logged_in(gcloud: str) -> bool:
    """Check if the user is logged in to Google Cloud Storage.

    Args:
        gcloud: Path to the `gcloud` executable.

    Returns:
        bool: True if the user is logged in, False otherwise.
    """
    out = run_command([gcloud, "auth", "list"], stdout="RETURN")
    return "ACTIVE" in out  # type: ignore


def is_valid_gs_bucket_url(url: str) -> bool:
    """Check if a URL is a valid Google Cloud Storage bucket URL.

    Such as `gs://bucket`.
    """
    if not url.startswith("gs://"):
        return False

    url = url.rstrip("/")
    return "/" not in url[5:]


def get_file_path(url: str) -> str:
    """Get the file path from a Google Cloud Storage file URL, without bucket.

    For example: gs://bucket/path/to/file -> path/to/file

    Args:
        url: The Google Cloud Storage file URL.

    Returns:
        str: The file path.
    """
    return url[5:].split("/", 1)[1]


def is_valid_gs_file_url(url: str) -> bool:
    """Check if a URL is a valid Google Cloud Storage file URL.

    Such as `gs://bucket/path/to/file`.
    """
    return url.startswith("gs://")
