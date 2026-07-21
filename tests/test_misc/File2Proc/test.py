from pathlib import Path

from biopipen.ns.misc import File2Proc as File2Proc_
from biopipen.core.testing import get_pipeline


class File2ProcSymlink(File2Proc_):
    ...


class File2ProcCopy(File2Proc_):
    envs = {"copy": True}


def pipeline(**kwargs):
    return (
        get_pipeline(__file__, **kwargs)
        .set_start(File2ProcSymlink, File2ProcCopy)
        .set_data(__file__, __file__)
    )


def testing(pipen):
    # assert pipen._succeeded
    outfile_symlink = pipen.procs[-2].workdir.joinpath(
        "0", "output", Path(__file__).name
    )
    assert outfile_symlink.exists()
    assert outfile_symlink.is_symlink()

    # copied
    outfile_copied = pipen.procs[-1].workdir.joinpath(
        "0", "output", Path(__file__).name
    )
    assert outfile_copied.exists()
    assert not outfile_copied.is_symlink()


if __name__ == "__main__":
    pipen = pipeline(loglevel="debug")
    assert pipen.run()
    testing(pipen)
