"""Tests for the Conda test runner's environment selection."""

import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from unittest import TestCase, main


RUNNER = Path(__file__).with_name("run_test.sh")


class TestRunTest(TestCase):
    def test_environment_selection(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            proc = root / "process"
            proc.mkdir()
            conda = root / "conda"
            conda.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
            conda.chmod(0o755)
            env = dict(os.environ, PATH=f"{root}{os.pathsep}{os.environ['PATH']}")

            def run():
                return subprocess.run(
                    ["bash", str(RUNNER), str(proc), "VERBOSE=true"],
                    env=env,
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout

            self.assertIn("\n-n\nbiopipen\npython\n", run())
            (proc / "run.env").write_text('ENVNAME="custom"\n')
            self.assertIn("\n-n\ncustom\npython\n", run())


if __name__ == "__main__":
    main()
