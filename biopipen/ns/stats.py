"""Provides processes for statistics."""

from ..core.proc import Proc
from ..core.config import config


class ChowTest(Proc):
    """Massive Chow tests.

    See Also https://en.wikipedia.org/wiki/Chow_test

    Input:
        infile: The input data file. The rows are samples and the columns are
            features. It must be tab-delimited.
            ```
            Sample   F1   F2   F3   ...   Fn
            S1       1.2  3.4  5.6        7.8
            S2       2.3  4.5  6.7        8.9
            ...
            Sm       5.6  7.8  9.0        1.2
            ```
        groupfile: The group file. The rows are the samples and the columns
            are the groupings. It must be tab-delimited.
            ```
            Sample   G1   G2   G3   ...   Gk
            S1       0    1    0          0
            S2       2    1    0          NA  # exclude this sample
            ...
            Sm       1    0    0          0
            ```
        fmlfile: The formula file. The first column is grouping and the
            second column is the formula. It must be tab-delimited.
            ```
            Group   Formula  ...  # Other columns to be added to outfile
            G1      Fn ~ F1 + Fx + Fy   # Fx, Fy could be covariates
            G1      Fn ~ F2 + Fx + Fy
            ...
            Gk      Fn ~ F3 + Fx + Fy
            ```

    Output:
        outfile: The output file. It is a tab-delimited file with the first
            column as the grouping and the second column as the p-value.
            ```
            Group  Formula  ...  Pooled  Groups  SSR  SumSSR  Fstat  Pval  Padj
            G1     Fn ~ F1       0.123   2       1    0.123   0.123  0.123 0.123
            G1     Fn ~ F2       0.123   2       1    0.123   0.123  0.123 0.123
            ...
            Gk     Fn ~ F3       0.123   2       1    0.123   0.123  0.123 0.123
            ```

    Envs:
        padj (choice): The method for p-value adjustment.
            - none: No p-value adjustment (no Padj column in outfile).
            - holm: Holm-Bonferroni method.
            - hochberg: Hochberg method.
            - hommel: Hommel method.
            - bonferroni: Bonferroni method.
            - BH: Benjamini-Hochberg method.
            - BY: Benjamini-Yekutieli method.
            - fdr: FDR correction method.
        transpose_input (flag): Whether to transpose the input file.
        transpose_group (flag): Whether to transpose the group file.
    """
    input = "infile:file, groupfile:file, fmlfile:file"
    output = "outfile:file:{{in.infile | stem}}.chowtest.txt"
    lang = config.lang.rscript
    envs = {
        "padj": "none",
        "transpose_input": False,
        "transpose_group": False,
    }
    script = "file://../scripts/stats/ChowTest.R"


