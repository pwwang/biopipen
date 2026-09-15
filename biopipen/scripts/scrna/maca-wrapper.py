"""Run MACA on an AnnData and save the per-cell labels.

Used by CellTypeAnnotation-maca.R (through
biopipen.utils::RunCellTypeAnnotation()), the same way as the other
python-based annotation tools.

Prerequisites
-------------
MACA is on PyPI (`MACA-Python`), but it pins `scanpy==1.6.0` and
`anndata==0.7.5`, which cannot coexist with a modern scanpy, so it needs its
own environment, for example:

    conda create -n maca python=3.9
    conda activate maca
    pip install scanpy==1.6.0 anndata==0.7.5 scikit-learn MACA-Python

and this wrapper has to be run with that environment's python (the
`envs.maca.python` of the process). MACA is imported lazily and a clear error
is raised naming the environment when the import fails.

`maca.singleMACA()` scores the cells against the marker sets (`cell_markers`,
a dict of cell type to marker genes, built here from the 2-column marker file),
clusters the score matrix at every combination of `res` and `n_neis`, and maps
the clusters to cell types by the majority of their cells. The per-cell labels
are aligned with the rows of the AnnData.
"""
from argparse import ArgumentParser
import sys


def csv_list(value, cast):
    """Parse a comma/space separated list, e.g. `1,2,3`."""
    if value is None:
        return None
    return [cast(item) for item in value.replace(",", " ").split()]


def main():
    parser = ArgumentParser(description="Run MACA")
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
        "--n-pcs", type=int, default=None,
        help="Number of principal components of the marker scores"
    )
    parser.add_argument(
        "--res", default=None,
        help="Comma separated Louvain resolutions (default: MACA's [1, 2, 3])"
    )
    parser.add_argument(
        "--n-neis", default=None,
        help="Comma separated numbers of neighbors (default: MACA's [5, 10])"
    )
    parser.add_argument(
        "--freq", type=float, default=None,
        help="Frequency threshold of the cluster mapping (default: 0.5)"
    )
    parser.add_argument(
        "--use-weight", action="store_true",
        help="Weight the markers by their order in the marker file"
    )
    args = parser.parse_args()

    try:
        import maca
    except ImportError as exc:
        sys.exit(
            f"Cannot import MACA ({exc}). MACA pins scanpy==1.6.0 and "
            "anndata==0.7.5, so it needs its own environment:\n"
            "  pip install scanpy==1.6.0 anndata==0.7.5 scikit-learn MACA-Python\n"
            "and the process has to run with that environment's python "
            "(`envs.maca.python`)."
        )
    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(args.input)
    markers = pd.read_csv(
        args.marker, sep="\t", header=None, names=["cell_type", "gene"]
    )
    cell_markers = {
        str(cell_type): genes.astype(str).tolist()
        for cell_type, genes in markers.groupby("cell_type")["gene"]
    }
    print(
        f"MACA: annotating {adata.n_obs} cells with "
        f"{len(cell_markers)} cell type(s)"
    )

    # None is MACA's own default for these: leave them out so the defaults stay
    # in one place (and `res`/`n_neis` must be lists)
    kwargs = {
        "n_pcs": args.n_pcs,
        "res": csv_list(args.res, float),
        "n_neis": csv_list(args.n_neis, int),
        "freq": args.freq,
        "use_weight": args.use_weight,
    }
    kwargs = {key: value for key, value in kwargs.items() if value is not None}

    _, labels = maca.singleMACA(
        ad=adata, cell_markers=cell_markers, **kwargs
    )
    pd.DataFrame(
        {"cell": adata.obs_names, "maca_celltype": labels}
    ).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
