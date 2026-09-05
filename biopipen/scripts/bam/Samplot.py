from __future__ import annotations

import os
import warnings
import shutil
from filelock import FileLock
from pathlib import Path
from subprocess import check_output
from biopipen.utils.misc import run_command, dict_to_cli_args

bamfiles: list[str] = {{ in.bamfiles | each: str | repr }}  # pyright: ignore # noqa
outfile: str = {{ out.outfile | quote }}  # pyright: ignore

envs: dict = {{ envs | dict | repr }}  # pyright: ignore
samplot: str = envs.pop("samplot", "samplot")
titles: list[str] = envs.pop("titles", [])
chrom: str | None = envs.pop("chrom", None)
start: int | None = envs.pop("start", None)
end: int | None = envs.pop("end", None)

if not titles:
    titles = [Path(bamfile).stem for bamfile in bamfiles]

if len(titles) != len(bamfiles):
    raise ValueError("The number of titles must match the number of bamfiles.")

envs["b"] = bamfiles
envs["o"] = outfile
envs["titles"] = titles
envs["chrom"] = chrom
envs["start"] = start
envs["end"] = end

if not chrom or start is None or end is None:
    raise ValueError("chrom, start, and end must be specified.")

# patch samplot to fix:
# Traceback (most recent call last):
#   File "/user/someone//bin/samplot", line 10, in <module>
#     sys.exit(main())
#              ^^^^^^
#   File "/user/someone//lib/python3.12/site-packages/samplot/__main__.py", line 31, in main
#     args.func(parser, args, extra_args)
#   File "/user/someone//lib/python3.12/site-packages/samplot/samplot.py", line 3576, in plot
#     read_data, max_coverage = get_read_data(
#                               ^^^^^^^^^^^^^^
#   File "/user/someone//lib/python3.12/site-packages/samplot/samplot.py", line 2773, in get_read_data
#     read_data["all_pairs"] = downsample_pairs(
#                              ^^^^^^^^^^^^^^^^^
#   File "/user/someone//lib/python3.12/site-packages/samplot/samplot.py", line 2790, in downsample_pairs
#     all_pairs[bam_i][hp_i] = sample_normal(
#                              ^^^^^^^^^^^^^^
#   File "/user/someone//lib/python3.12/site-packages/samplot/samplot.py", line 539, in sample_normal
#     for read_name in random.sample(inside_norm.keys(), max_depth):
#                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/user/someone//lib/python3.12/random.py", line 413, in sample
#     raise TypeError("Population must be a sequence.  "
# TypeError: Population must be a sequence.  For dicts or sets, use sorted(d).

samplot_path = shutil.which(samplot)
if samplot_path is None:
    raise ValueError(f"samplot executable not found: {samplot}")

# what's in samplot:
# #!/bin/sh
# '''exec' /home/pwwang/miniconda3/bin/python "$0" "$@"
# ' '''
# # -*- coding: utf-8 -*-
# import re
# import sys
# from samplot.__main__ import main
# if __name__ == '__main__':
#     sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
#     sys.exit(main())

# Get python path from samplot executable file

python_path = None
with open(samplot_path, "r") as f:
    for line in f:
        if line.startswith("'''exec'"):
            python_path = line.split()[1]
            break

if python_path is None:
    warnings.warn(
        "Could not find python path in samplot executable. "
        "The samplot patch may not be applied."
    )
else:
    samplot_file = check_output(
        [
            python_path,
            "-c",
            "from samplot import samplot; print(samplot.__file__)",
        ]
    )
    # do we have write permission to the samplot file? If not, we cannot patch it.
    samplot_file = Path(samplot_file.decode().strip())
    if not samplot_file.exists():
        warnings.warn(
            f"samplot file does not exist: {samplot_file}. "
            "The samplot patch may not be applied."
        )

    elif not os.access(samplot_file, os.W_OK):
        warnings.warn(
            f"No write permission to the samplot file: {samplot_file}. "
            "The samplot patch may not be applied."
        )

    else:
        with FileLock(f"/tmp/{samplot_file.name}.lock"):
            old_code = code = samplot_file.read_text()
            changed = False

            # Fix Python 3.12 TypeError in random.sample
            old = "inside_norm.keys()"
            new = "sorted(inside_norm.keys())"
            if old in code and new not in code:
                code = code.replace(old, new)
                changed = True

            # Fix overlapping axis text with matplotlib >= 3.7
            # (samplot pins matplotlib <3.7 upstream, see ryanlayer/samplot#204)
            # 1. Hide the duplicate x tick labels on the coverage twin axis,
            #    which shares the x-axis formatter with the insert-size axis
            old = '    ax2.tick_params(axis="x", length=0)'
            new = '    ax2.tick_params(axis="x", length=0, labelbottom=False)'
            if old in code and new not in code:
                code = code.replace(old, new)
                changed = True

            # 2. Blank the unused scaffold axis that leaks default
            #    matplotlib tick labels over the real ones
            old = (
                "        #ax is never used, annotating this for readability\n"
                "        ax = plt.subplot(grid[ax_i])"
            )
            new = old + "\n        ax.set_xticklabels([])\n        ax.set_yticklabels([])"
            if old in code and new not in code:
                code = code.replace(old, new)
                changed = True

            # 3. Pin the x tick positions (FixedLocator) before setting the
            #    labels, so they always match the drawn ticks
            old = (
                "            curr_ax.set_xticklabels("
                "labels, fontsize=xaxis_label_fontsize)"
            )
            new = (
                "            curr_ax.set_xticks("
                "curr_ax.xaxis.get_majorticklocs())\n"
                "            curr_ax.set_xticklabels("
                "labels, fontsize=xaxis_label_fontsize)"
            )
            if old in code and new not in code:
                code = code.replace(old, new)
                changed = True

            if changed:
                samplot_file.with_suffix(".py.bak").write_text(old_code)
                samplot_file.write_text(code)

cmd = [
    samplot,
    "plot",
    *dict_to_cli_args(envs, dashify=False, dup_key=False),
]

run_command(cmd)
