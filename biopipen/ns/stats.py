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


class LiquidAssoc(Proc):
    """Liquid association tests.

    See Also https://github.com/gundt/fastLiquidAssociation
    Requieres https://github.com/pwwang/fastLiquidAssociation

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
            The features (columns) will be tested pairwise, which will be the X and
            Y columns in the result of `fastMLA`
        covfile: The covariate file. The rows are the samples and the columns
            are the covariates. It must be tab-delimited.
            If provided, the data in `in.infile` will be adjusted by covariates by
            regressing out the covariates and the residuals will be used for
            liquid association tests.
        groupfile: The group file. The rows are the samples and the columns
            are the groupings. It must be tab-delimited.
            ```
            Sample   G1   G2   G3   ...   Gk
            S1       0    1    0          0
            S2       2    1    0          NA  # exclude this sample
            ...
            Sm       1    0    0          0
            ```
            This will be served as the Z column in the result of `fastMLA`
            This can be omitted. If so, `envs.nvec` should be specified, which is
            to select column from `in.infile` as Z.
        fmlfile: The formula file. The 3 columns are X3, X12 and X21. The results
            will be filtered based on the formula. It must be tab-delimited without
            header.

    Output:
        outfile: The output file.
            ```
            X12  X21  X3  rhodiff  MLA value  estimates  san.se  wald  Pval  model
            C38  C46  C5  0.87  0.32  0.67  0.20  10.87  0  F
            C46  C38  C5  0.87  0.32  0.67  0.20  10.87  0  F
            C27  C39  C4  0.94  0.34  1.22  0.38  10.03  0  F
            ```

    Envs:
        nvec: The column index (1-based) of Z in `in.infile`, if `in.groupfile` is
            omitted. You can specify multiple columns by comma-seperated values, or
            a range of columns by `-`. For example, `1,3,5-7,9`. It also supports
            column names. For example, `F1,F3`. `-` is not supported for column
            names.
        x: Similar as `nvec`, but limit X group to given features.
            The rest of features (other than X and Z) in `in.infile` will
            be used as Y.
            The features in `in.infile` will still be tested pairwise, but only
            features in X and Y will be kept.
        topn (type=int): Number of results to return by `fastMLA`, ordered from
            highest `|MLA|` value descending.
            The default of the package is 2000, but here we set to 1e6 to return as
            many results as possible (also good to do pvalue adjustment).
        rvalue (type=float): Tolerance value for LA approximation. Lower values of
            rvalue will cause a more thorough search, but take longer.
        cut (type=int): Value passed to the GLA function to create buckets
            (equal to number of buckets+1). Values placing between 15-30 samples per
            bucket are optimal. Must be a positive integer>1. By default,
            `max(ceiling(nrow(data)/22), 4)` is used.
        ncores (type=int): Number of cores to use for parallelization.
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
        transpose_cov (flag): Whether to transpose the covariate file.
        xyz_names: The names of X12, X21 and X3 in the final output file. Separated
            by comma. For example, `X12,X21,X3`.
    """
    input = "infile:file, covfile:file, groupfile:file, fmlfile:file"
    output = "outfile:file:{{in.infile | stem}}.liquidassoc.txt"
    lang = config.lang.rscript
    envs = {
        "nvec": None,
        "x": None,
        "topn": 1e6,
        "rvalue": 0.5,
        "cut": 20,
        "ncores": config.misc.ncores,
        "padj": "none",
        "transpose_input": False,
        "transpose_group": False,
        "transpose_cov": False,
        "xyz_names": None,
    }
    script = "file://../scripts/stats/LiquidAssoc.R"


