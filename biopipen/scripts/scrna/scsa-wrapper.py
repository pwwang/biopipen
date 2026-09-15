"""Run SCSA on an AnnData and save the per-cell labels.

Used by CellTypeAnnotation-scsa.R (through
biopipen.utils::RunCellTypeAnnotation()), the same way as the other
python-based annotation tools.

Prerequisites
-------------
SCSA is not on CRAN/Bioconductor/PyPI, so the repo has to be cloned manually
(https://github.com/bioinfo-ibms-pumc/SCSA); `--scsa-dir` points to the clone,
which must hold `SCSA.py` and its reference database `whole.db`. The python of
`--python` (default: the python running this script) needs the SCSA
dependencies: `pandas`, `numpy`, `scipy` and `openpyxl`.

The input of SCSA.py
--------------------
`SCSA.py -i` is *not* an expression matrix (and it has no h5ad support): it is
a per-cluster marker (differential expression) table, with `Cluster`, `Gene`,
`avg_logFC` and `p_val_adj` columns for `-s seurat`, which SCSA filters with
`-f`/`-p` and tests against the marker sets. So this wrapper computes the
table from the h5ad with scanpy's `rank_genes_groups` (the h5ad is the query
matrix here, as for the other h5ad tools), writes it next to the other
intermediate files, and passes it to SCSA.py with `-s seurat -E` (the marker
genes, in the universal marker table and in the h5ad, are gene symbols).

The clusters come from `--ident` (a metadata column of the h5ad); when it is
not given, the only categorical column of the h5ad is used.

SCSA labels the clusters, so the per-cell labels are the labels of the
clusters the cells belong to. `SCSA.py` exits with status 0 even when it
cannot do anything (e.g. no marker survives the fold-change/p-value
thresholds), so the output file is checked here instead.
"""
from argparse import ArgumentParser
import os
import subprocess
import sys
import tempfile


def cluster_column(adata, ident):
    """The metadata column holding the clusters to annotate."""
    if ident:
        if ident not in adata.obs.columns:
            sys.exit(
                f"Cluster column '{ident}' is not in the h5ad: "
                f"{', '.join(map(str, adata.obs.columns))}"
            )
        return ident

    # No --ident: the h5ad of a Seurat object carries all of its metadata, so
    # the only way to pick the clusters is the only categorical column
    cats = [
        col for col in adata.obs.columns
        if str(adata.obs[col].dtype) == "category"
    ]
    if len(cats) != 1:
        sys.exit(
            "Cannot tell which metadata column holds the clusters "
            f"(candidate column(s): {', '.join(map(str, cats)) or 'none'}). "
            "Pass it with --ident."
        )
    return cats[0]


def write_deg_table(adata, ident, path):
    """SCSA's `-i` input: the per-cluster markers of the h5ad, in the `seurat`
    convention of SCSA.py (`cluster`, `gene`, `avg_logFC`, `p_val_adj`)."""
    import numpy as np
    import scanpy as sc

    values = adata.X.data if hasattr(adata.X, "data") else np.asarray(adata.X)
    if values.size and np.isfinite(values.max()) and np.all(values % 1 == 0):
        # raw counts (log1p-normalized data is not integral)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    sc.tl.rank_genes_groups(adata, groupby=ident, method="wilcoxon")
    deg = sc.get.rank_genes_groups_df(adata, group=None)
    deg = deg.rename(
        columns={
            "names": "gene",
            "logfoldchanges": "avg_logFC",
            "pvals_adj": "p_val_adj",
            # the cluster column is `group` in scanpy >= 1.9
            "group" if "group" in deg.columns else "cluster": "cluster",
        }
    )
    deg = deg[["cluster", "gene", "avg_logFC", "p_val_adj"]]
    # SCSA.py's parser strips the ".1" suffixes of the gene column
    # (`str.replace("\\.\d+", "")`), which needs strings
    deg["gene"] = deg["gene"].astype(str)
    deg.to_csv(path, index=False)


def main():
    parser = ArgumentParser(description="Run SCSA")
    parser.add_argument(
        "-i", "--input", required=True, help="Input H5AD file (AnnData)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file (barcode, cell type)"
    )
    parser.add_argument(
        "-m", "--marker", required=True,
        help="Marker file: a 2-column headerless TSV of cell type and gene"
    )
    parser.add_argument(
        "--scsa-dir", required=True,
        help="The path to the cloned SCSA repo (holding SCSA.py and whole.db)"
    )
    parser.add_argument(
        "--ident", default=None,
        help="The metadata column of the h5ad holding the clusters "
             "(default: the only categorical column)"
    )
    parser.add_argument(
        "--norefdb", action="store_true",
        help="Use only the marker file, not SCSA's reference database"
    )
    parser.add_argument("-g", "--species", default=None, help="Species")
    parser.add_argument("-k", "--tissue", default=None, help="Tissue")
    parser.add_argument(
        "--python", default=sys.executable,
        help="The python to run SCSA.py with (default: this python)"
    )
    args = parser.parse_args()

    scsa_py = os.path.join(args.scsa_dir, "SCSA.py")
    if not os.path.exists(scsa_py):
        sys.exit(
            f"No SCSA.py in '{args.scsa_dir}'. SCSA is not on "
            "CRAN/Bioconductor/PyPI: clone "
            "https://github.com/bioinfo-ibms-pumc/SCSA and pass the clone "
            "with `envs.scsa.scsa_dir`."
        )
    if not os.path.exists(args.marker):
        sys.exit(f"Marker file does not exist: {args.marker}")

    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(args.input)
    ident = cluster_column(adata, args.ident)
    print(
        f"SCSA: annotating {adata.n_obs} cells in "
        f"{adata.obs[ident].nunique()} cluster(s) of '{ident}'"
    )

    workdir = tempfile.mkdtemp(prefix="scsa-")
    deg_file = os.path.join(workdir, "deg.csv")
    out_file = os.path.join(workdir, "scsa.txt")
    write_deg_table(adata, ident, deg_file)

    command = [
        args.python, scsa_py,
        "-i", deg_file,
        "-o", out_file,
        "-s", "seurat",
        "-E",  # the genes are symbols, not ensembl IDs
        "-m", "txt",
        "-M", args.marker,
        # check_db() resolves `whole.db` relative to SCSA.py but opens it
        # relative to the cwd, so use the absolute path
        "-d", os.path.join(args.scsa_dir, "whole.db"),
    ]
    if args.norefdb:
        command.append("-N")
    if args.species:
        command += ["-g", args.species]
    if args.tissue:
        command += ["-k", args.tissue]

    print(f"SCSA: running {' '.join(command)}")
    proc = subprocess.run(
        command, cwd=args.scsa_dir, capture_output=True, text=True
    )
    print(proc.stdout, end="", file=sys.stderr)
    if proc.returncode != 0:
        sys.exit(f"SCSA.py failed with status {proc.returncode}:\n{proc.stderr}")
    if not os.path.exists(out_file):
        # SCSA.py exits with status 0 when it cannot annotate anything (e.g. no
        # marker passes -f/-p, or no marker matches the reference database)
        sys.exit(
            "SCSA.py produced no result. Check the markers (and, without "
            "`--norefdb`, that `whole.db` matches `--species`/`--tissue`):\n"
            f"{proc.stdout}"
        )

    # "Cell Type\tZ-score\tCluster", one row per cluster
    result = pd.read_csv(out_file, sep="\t")
    if not {"Cell Type", "Cluster"} <= set(result.columns):
        sys.exit(f"Unexpected SCSA output columns: {list(result.columns)}")
    cluster2celltype = dict(
        zip(result["Cluster"].astype(str), result["Cell Type"].astype(str))
    )
    print(f"SCSA: annotated {len(cluster2celltype)} cluster(s)")

    labels = adata.obs[ident].astype(str).map(cluster2celltype)
    pd.DataFrame(
        {
            "cell": adata.obs_names,
            "scsa_celltype": labels.fillna("unassigned").to_numpy(),
        }
    ).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
