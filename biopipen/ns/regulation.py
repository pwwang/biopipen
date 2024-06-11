"""Provides processes for the regulation related"""

from ..core.proc import Proc
from ..core.config import config


class MotifScan(Proc):
    """Scan the input sequences for binding sites using motifs.

    Currently only [fimo](https://meme-suite.org/meme/tools/fimo) from MEME suite
    is supported, based on the research/comparisons done by the following reference.

    Reference:
        - [Evaluating tools for transcription factor binding site prediction](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6889335/)

    Input:
        motiffile: File containing motif names.
            The file contains the motif and regulator names.
            The motif names should match the names in the motif database.
            This file must have a header.
            If multiple columns are present, it should be delimited by tab.
        seqfile: File containing sequences in FASTA format.

    Output:
        outdir: Directory containing the results.
            Especially `fimo_output.txt` extending from `fimo.tsv`, which contains:
            1. the results with the regulator information if `envs.regulator_col`
                is provided, otherwise, the `regulator` columns will be filled with
                the motif names.
            2. the original sequence from the fasta file (in.seqfile)
            3. corrected genomic coordinates if the genomic coordinates are included
                in the sequence names.

            See also the `Output` section of
            <https://meme-suite.org/meme/doc/fimo.html>.
            Note that `--no-pgc` is passed to fimo to not parse the genomic coordinates
            from the sequence names by fimo. When fimo parses the genomic coordinates,
            `DDX11L1` in `>DDX11L1::chr1:11869-14412` will be lost.
            The purpose of this is to keep the sequence names as they are in the output.
            If the sequence names are in the format of `>NAME::chr1:START-END`, we will
            correct the coordinates in the output.
            Also note that it requires meme/fimo v5.5.5+ to do this
            (where the --no-pgc option is available).

    Envs:
        tool (choice): The tool to use for scanning.
            Currently only fimo is supported.
            - fimo: Use fimo from MEME suite.
        fimo: The path to fimo binary.
        motif_col: The column name in the motif file containing the motif names.
        regulator_col: The column name in the motif file containing the regulator names.
            Both `motif_col` and `regulator_col` should be the direct column names or
            the index (1-based) of the columns.
            If no `regulator_col` is provided, no regulator information is written in
            the output.
        notfound (choice): What to do if a motif is not found in the database.
            - error: Report error and stop the process.
            - ignore: Ignore the motif and continue.
        motifdb: The path to the motif database. This is required
            It should be in the format of MEME motif database.
            Databases can be downloaded here: <https://meme-suite.org/meme/doc/download.html>.
            See also introduction to the databases: <https://meme-suite.org/meme/db/motifs>.
        cutoff (type=float): The cutoff for p-value to write the results.
            When `envs.q_cutoff` is set, this is applied to the q-value.
            This is passed to `--thresh` in fimo.
        q (flag): Calculate q-value.
            When `False`, `--no-qvalue` is passed to fimo.
            The q-value calculation is that of Benjamini and Hochberg (BH) (1995).
        q_cutoff (flag): Apply `envs.cutoff` to q-value.
        args (ns): Additional arguments to pass to the tool.
            - <more>: Additional arguments for fimo.
                See: <https://meme-suite.org/meme/doc/fimo.html>
    """  # noqa: E501
    input = "motiffile:file, seqfile:file"
    output = "outdir:dir:{{in.motiffile | stem}}.fimo"
    lang = config.lang.python
    envs = {
        "tool": "fimo",
        "fimo": config.exe.fimo,
        "motif_col": 1,
        "regulator_col": None,
        "notfound": "error",
        "motifdb": config.tf_motifdb,
        "cutoff": 1e-4,
        "q": False,
        "q_cutoff": False,
        "args": {},
    }
    script = "file://../scripts/regulation/MotifScan.py"
