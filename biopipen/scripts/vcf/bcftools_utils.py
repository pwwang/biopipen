"""Utilities for bcftools"""

from biopipen.utils.misc import run_command, dict_to_cli_args
from biopipen.utils.reference import tabix_index


def bcftools_version(bcftools: str) -> tuple[int, ...]:
    """Get the version of bcftools

    Args:
        bcftools (str): Path to bcftools

    Returns:
        tuple[int, ...]: The version of bcftools
    """
    bversion = (
        run_command([bcftools, "version"], stdout="return")
        .splitlines()[0]  # bcftools 1.20  # type: ignore
        .replace("bcftools", "")
        .strip()  # 1.20
        .split(".")
    )
    return tuple(map(int, bversion))


def run_bcftools(
    args: dict,
    bcftools: str,  # TODO: get from the first argument of args
    index: bool,
    tabix: str
) -> None:
    """Run bcftools with the given arguments

    Args:
        args: Arguments to pass to bcftools
        bcftools (str): Path to bcftools
        index (bool): Whether to index the output
        tabix (str): Path to tabix
    """
    if not index:
        run_command(dict_to_cli_args(args, dashify=True), fg=True)
    else:
        bversion = bcftools_version(bcftools)
        if bversion >= (1, 20):
            # requires bcftools 1.20+
            # '--write-index tbi' not working
            # it has to be '--write-index=tbi'
            args["write_index=tbi"] = True
            run_command(dict_to_cli_args(args, dashify=True), fg=True)
        else:
            run_command(dict_to_cli_args(args, dashify=True), fg=True)
            tabix_index(args["o"], "vcf", tmpdir=False, tabix=tabix)
