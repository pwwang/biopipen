"""Run SCSA on an AnnData and save the per-cell labels.

Used by CellTypeAnnotation-scsa.R (through
biopipen.utils::RunCellTypeAnnotation()), the same way as the other
python-based annotation tools.

Prerequisites
-------------
SCSA (https://github.com/bioinfo-ibms-pumc/SCSA) is not on
CRAN/Bioconductor/PyPI, and its `SCSA.py` is 2019-era python: it uses
`numpy.asfarray`/`numpy.mat` and `DataFrame.append` (both removed by numpy 2 /
pandas 2), and its reference database `whole.db` is a pandas<1.0 pickle that
modern pandas cannot unpickle (it is read even with `--norefdb`). The scoring
of its `-s seurat` path is therefore ported here, as maca-wrapper.py ports MACA,
so biopipen runs SCSA without a clone of the repo -- and without the reference
database, which only the user marker table (`-N -M` of `SCSA.py`) needs.

The port keeps SCSA's scoring and its thresholds (`-f`/`--foldchange` and
`-p`/`--pvalue`, with the defaults of `SCSA.py`): every cell type of the marker
table is scored against the cluster's markers, weighted by their fold change,
and the best scored cell type labels the cluster. SCSA prints the two top cell
types when they are close (`A|B`), but that is only its stdout report: the cell
type it labels the cluster with is the first one, which is what is used here.

The input of SCSA
-----------------
`SCSA.py -i` is *not* an expression matrix (and it has no h5ad support): it is
a per-cluster marker (differential expression) table, with `Cluster`, `Gene`,
`avg_logFC` and `p_val_adj` columns for `-s seurat`. So this wrapper computes
the table from the h5ad with scanpy's `rank_genes_groups` (the h5ad is the
query matrix here, as for the other h5ad tools), and scores it with the marker
genes as gene symbols (`-E` of `SCSA.py`).

The clusters come from `--ident` (a metadata column of the h5ad); when it is
not given, the only categorical column of the h5ad is used.

SCSA labels the clusters, so the per-cell labels are the labels of the clusters
the cells belong to. A cluster with no marker passing the thresholds (or with
none of them in the marker table) is left unannotated.
"""
from argparse import ArgumentParser
import sys

# SCSA's own defaults for the input markers (`-f`/`-p` of `SCSA.py`)
FOLDCHANGE = 2.0
PVALUE = 0.05


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


def deg_table(adata, ident):
    """SCSA's input: the per-cluster markers of the h5ad, in the `seurat`
    convention of SCSA.py (`cluster`, `gene`, `avg_logFC`, `p_val_adj`)."""
    import numpy as np
    import scanpy as sc

    # `.data` holds only the stored non-zeros of a sparse matrix (the implicit
    # zeros are integral anyway); anything else has to be materialized: a dense
    # `ndarray.data` is a memoryview (no `.size`/`.max()`) and a backed X is an
    # h5py dataset
    values = np.asarray(adata.X.data if hasattr(adata.X, "tocsr") else adata.X)
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
    # SCSA matches the markers by gene name, so the comparison has to be on
    # strings, whatever the h5ad holds in `var_names`
    deg["gene"] = deg["gene"].astype(str)
    return deg


def scsa_labels(markers, deg, foldchange=FOLDCHANGE, pvalue=PVALUE):
    """The cell type of every cluster: the port of SCSA's `-s seurat -N -M`."""
    import numpy as np
    import pandas as pd

    markers = markers[["cellName", "gene"]]
    cell_types = sorted(markers["cellName"].unique())
    deg = deg[(deg["avg_logFC"] >= foldchange) & (deg["p_val_adj"] <= pvalue)]

    labels = {}
    for cluster, cluster_deg in deg.groupby("cluster", observed=True):
        genes = sorted(set(cluster_deg["gene"]) & set(markers["gene"]))
        if not genes:
            continue
        # one column per marker gene the cluster has, one row per cell type of
        # the marker table: how many of those genes the cell type shares with
        # the cluster, log2-scaled (SCSA's `log2(count + 0.05)`, which the
        # sparse matrix of SCSA applies to the shared genes only)
        shared = markers[markers["gene"].isin(genes)]
        cell_matrix = (
            pd.crosstab(shared["cellName"], shared["gene"])
            .reindex(index=cell_types, columns=genes, fill_value=0)
            .to_numpy(dtype=float)
        )
        cell_matrix[cell_matrix > 0] = np.log2(cell_matrix[cell_matrix > 0] + 0.05)
        # the genes are weighted by their fold change, scaled by the mean fold
        # change of the cluster's markers (SCSA's `gene_matrix * mean`)
        logfc = cluster_deg.set_index("gene").loc[genes, "avg_logFC"]
        logfc = logfc.to_numpy(dtype=float)
        scores = pd.Series(cell_matrix @ (logfc * logfc.mean()), index=cell_types)
        scores = scores.abs().sort_values(ascending=False)
        # SCSA z-scores the scores across the cell types (with more than one)
        if len(scores) > 1 and scores.std(ddof=1) > 0:
            scores = (scores - scores.mean()) / scores.std(ddof=1)
        labels[str(cluster)] = scores.index[0]
    return labels


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
        "--ident", default=None,
        help="The metadata column of the h5ad holding the clusters "
             "(default: the only categorical column)"
    )
    parser.add_argument(
        "-f", "--foldchange", type=float, default=FOLDCHANGE,
        help=f"Minimum fold change of the markers (default: {FOLDCHANGE})"
    )
    parser.add_argument(
        "-p", "--pvalue", type=float, default=PVALUE,
        help=f"Maximum adjusted p-value of the markers (default: {PVALUE})"
    )
    args = parser.parse_args()

    import pandas as pd
    import scanpy as sc

    markers = pd.read_csv(
        args.marker, sep="\t", header=None, names=["cellName", "gene"],
        dtype=str,
    )
    adata = sc.read_h5ad(args.input)
    ident = cluster_column(adata, args.ident)
    print(
        f"SCSA: annotating {adata.n_obs} cells in "
        f"{adata.obs[ident].nunique()} cluster(s) of '{ident}'"
    )

    clusters = adata.obs[ident].astype(str)
    labels = scsa_labels(
        markers, deg_table(adata, ident), args.foldchange, args.pvalue
    )
    if not labels:
        sys.exit(
            "SCSA annotated no cluster: no marker of "
            f"{args.marker} passes the fold change ({args.foldchange}) and "
            f"p-value ({args.pvalue}) thresholds of SCSA. Check the markers."
        )
    print(
        "SCSA: annotated "
        f"{len(labels)}/{clusters.nunique()} cluster(s): "
        + ", ".join(f"{clust}={label}" for clust, label in sorted(labels.items()))
    )

    pd.DataFrame({
        "cell": adata.obs_names,
        "scsa_celltype": clusters.map(labels).fillna("unassigned").to_numpy(),
    }).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
