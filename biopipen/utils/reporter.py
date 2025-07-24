from __future__ import annotations
from typing import Sequence
from os import PathLike
from pathlib import Path

"""An implementation of reporter in python
"https://pwwang.github.io/biopipen.utils.R/reference/Reporter.html

to generate a json file for pipen-report to build a report for a process.
"""

import json


class Reporter:

    def __init__(self):
        self.report = {}

    def add(
        self,
        *args,
        h1: str,
        h2: str = "#",
        h3: str = "#",
        ui: str = "flat",
    ) -> None:
        """Add a content to the report

        Args:
            *args: The content of the report
            h1 (str): The first level header
            h2 (str): The second level header
            h3 (str): The third level header
            ui (str): The user interface of the report
        """

        self.report.setdefault(h1, {})
        self.report[h1].setdefault(h2, {})
        self.report[h1][h2].setdefault(h3, {})
        self.report[h1][h2][h3][ui] = []

        for arg in args:
            self.report[h1][h2][h3][ui].append(arg)

    def add2(
        self,
        *args,
        hs: Sequence[str],
        hs2: Sequence[str] = (),
        ui: str = "flat",
        collapse: str = ": ",
    ) -> None:
        """Add a content to the report

        Args:
            *args: The content of the report
            hs: The headings of the case
            hs2: The headings that must be shown.
                When there are more items in `hs`, they will be concatenated.
                For example, if `hs = c("Section1", "Case1")`, and `hs2 = c("A", "B")`,
                then headings will be `h1 = "Section1: Case1"` and `h2 = "A"` and
                `h3 = "B"`.
            ui: The user interface of the report
            collapse: The separator to concatenate the headings
        """
        if len(hs2) > 2:
            raise ValueError("hs2 must have 2 or less items")

        if len(hs2) == 2:
            h1 = collapse.join(hs)
            h2 = hs2[0]
            h3 = hs2[1]
        elif len(hs2) == 1:
            h1 = hs[0]
            hs = hs[1:]
            if hs:
                h2 = collapse.join(hs)
                h3 = hs2[0]
            else:
                h2 = hs2[0]
                h3 = "#"
        else:
            h1 = hs[0]
            hs = hs[1:]
            if hs:
                h2 = hs[0]
                hs = hs[1:]
            else:
                h2 = "#"

            if hs:
                h3 = collapse.join(hs)
            else:
                h3 = "#"

        self.add(*args, h1=h1, h2=h2, h3=h3, ui=ui)

    def image(
        self,
        prefix: str,
        more_formats: str | Sequence[str],
        save_code: bool,
        kind: str = "image",
        **kwargs,
    ) -> dict:
        """Generate a report for an image to be added.

        Args:
            prefix: The prefix of the image.
            more_formats: More formats of the image available.
            save_code: Whether to save the code to reproduce the plot.
            kind: The kind of the report, default is "image".
            **kwargs: Other arguments to add to the report.

        Returns:
            dict: The structured report for the image

        Examples:
            >>> reporter = Reporter()
            >>> reporter.add(
            >>>   {
            >>>     "name": "Image 1",
            >>>     "contents": [
            >>>       reporter.image("/path/to/image1", "pdf", save_code=True)
            >>>     ]
            >>>   },
            >>>   h1="Images",
            >>>   h2="Image 1",
            >>> )
        """
        out = {
            "kind": kind,
            "src": f"{prefix}.png",
            **kwargs,
        }

        if more_formats or save_code:
            out["download"] = []

        if more_formats:
            for mf in more_formats:
                out["download"].append(f"{prefix}.{mf}")

        if save_code:
            out["download"].append(
                {
                    "src": f"{prefix}.code.zip",
                    "tip": "Download the code to reproduce the plot",
                    "icon": "Code",
                }
            )

        return out

    def clear(self):
        """Clear the report"""
        self.report = {}

    def save(self, path: str | PathLike, clear: bool = True) -> None:
        """Save the report to a file

        Args:
            path: The path to save the report
                If the path is a directory, the report will be saved as `report.json`
                in the directory. Otherwise, the report will be saved to the file.
            clear: Whether to clear the report after saving.
        """
        path = Path(path)
        if path.is_dir():
            path = path / "report.json"

        with open(path, "w") as f:
            json.dump(self.report, f, indent=2)

        if clear:
            self.clear()
