"""Run scMapNet on an AnnData and save the per-cell labels.

Used by CellTypeAnnotation-scmapnet.R (through
biopipen.utils::RunCellTypeAnnotation()), the same way as the other
python-based annotation tools.

Prerequisites
-------------
scMapNet is not a package; clone https://github.com/Yuz7/scMapNet and pass the
clone with `--scmapnet-dir`. Its pipeline is driven from there:

1. treemap image generation: `generate_image_script.sh` (needs `Rscript` with
   the `treemap`, `data.table`, `magick`, `SingleCellExperiment` packages) ->
   one PNG per cell under `<data_dir>/train/<label>/`.
2. prediction: `main_finetune.py --test` (needs `torch`, `timm`,
   `torchvision`, `scanpy`, `umap`) -> `<output_dir>/class_pred.csv`.

The pre-trained weights are a manual download (see the repo's README, the
Google Drive link) and are licensed **CC BY-NC 4.0 (non-commercial)**: they are
not part of the repo, so `--weights` is required and has to point to a
checkpoint the user downloaded and fine-tuned.

Notes on the pipeline
---------------------
`main_finetune.py --test` reads the images with `torchvision`'s ImageFolder,
so the predicted class of a cell is an index into the *checkpoint's* label
space — the cell types the model was fine-tuned on, in the order ImageFolder
saw them at training time, which is their sorted order. This wrapper assumes
that space is the sorted cell types of the marker file, and maps the indices
back to them.

The query cells have no labels, so all the images are generated under one
`label` (the directory scMapNet puts them in); the `target` column of
`class_pred.csv` is therefore meaningless and only `pred` is used.

The gene ids of scMapNet come from the `treemap/ensemble_ID_transfer_new.csv`
of the clone (ENSEMBL -> NCBI gene id -> symbol); the genes of the h5ad and the
marker genes are both mapped through it, to the NCBI gene ids scMapNet works
with.
"""
from argparse import ArgumentParser
import os
import subprocess
import sys
import tempfile

MAIN = "main_finetune.py"
IMAGE_SCRIPT = "generate_image_script.sh"
TRANSFER = os.path.join("treemap", "ensemble_ID_transfer_new.csv")
# scMapNet's own example; the model has to be the one the checkpoint was
# trained with
MODEL = "vit_large_patch16"


def run(command, cwd, description):
    print(f"scMapNet: {description}: {' '.join(command)}")
    proc = subprocess.run(command, cwd=cwd, text=True)
    if proc.returncode != 0:
        sys.exit(f"{description} failed with status {proc.returncode}")


def main():
    parser = ArgumentParser(description="Run scMapNet")
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
        "--scmapnet-dir", required=True,
        help="The path to the cloned scMapNet repo"
    )
    parser.add_argument(
        "--weights", required=True,
        help="The path to the fine-tuned checkpoint (.pth)"
    )
    parser.add_argument(
        "--organ", default=None,
        help="The organ of the cells, for the organ -> cell type -> gene "
             "treemap hierarchy"
    )
    parser.add_argument(
        "--python", default=sys.executable,
        help="The python to run main_finetune.py with (default: this python)"
    )
    parser.add_argument(
        "--ncores", type=int, default=1,
        help="The number of processes to generate the treemap images with "
             "(default: 1)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.scmapnet_dir):
        sys.exit(
            f"scMapNet directory does not exist: {args.scmapnet_dir}. Clone "
            "https://github.com/Yuz7/scMapNet and pass the clone with "
            "`envs.scmapnet.scmapnet_dir`."
        )
    for required in (MAIN, IMAGE_SCRIPT, TRANSFER):
        path = os.path.join(args.scmapnet_dir, required)
        if not os.path.exists(path):
            sys.exit(f"Not a scMapNet clone (no {required}): {args.scmapnet_dir}")
    if not os.path.exists(args.weights):
        sys.exit(
            f"scMapNet checkpoint does not exist: {args.weights}. The "
            "pre-trained weights are a manual download (they are not part of "
            "the repo, see https://github.com/Yuz7/scMapNet) and are licensed "
            "CC BY-NC 4.0 (non-commercial). Pass a fine-tuned checkpoint with "
            "`envs.scmapnet.weights`; it has to be fine-tuned on the cell "
            "types of the marker file."
        )
    if not os.path.exists(args.marker):
        sys.exit(f"Marker file does not exist: {args.marker}")

    import pandas as pd
    import scanpy as sc

    transfer = pd.read_csv(os.path.join(args.scmapnet_dir, TRANSFER))
    transfer.columns = ["ensgid", "geneid", "symbol"]
    symbol2geneid = dict(zip(transfer["symbol"], transfer["geneid"].astype(str)))

    adata = sc.read_h5ad(args.input)
    markers = pd.read_csv(
        args.marker, sep="\t", header=None, names=["cell_type", "gene"]
    )
    cell_types = sorted(markers["cell_type"].astype(str).unique())
    markers["id"] = markers["gene"].astype(str).map(symbol2geneid)
    n_dropped = int(markers["id"].isna().sum())
    markers = markers.dropna(subset=["id"])
    if markers.empty:
        sys.exit(
            "None of the marker genes is in scMapNet's gene id transfer table "
            f"({TRANSFER}); the marker genes have to be gene symbols."
        )
    if n_dropped:
        print(
            f"scMapNet: ignoring {n_dropped} marker gene(s) with no gene id "
            "in the transfer table"
        )

    workdir = tempfile.mkdtemp(prefix="scmapnet-")
    data_dir = os.path.join(workdir, "data")
    out_dir = os.path.join(workdir, "finetune")

    # A cell-by-gene CSV: `generate_image_script.sh -t df` reads the gene ids
    # from the column names and the labels from the metadata file. The images
    # are named `<label>_<column index>`, so the column order is the barcode
    # order used below.
    expression_file = os.path.join(workdir, "expression.csv")
    geneids = [symbol2geneid.get(str(gene)) for gene in adata.var_names]
    keep = [index for index, geneid in enumerate(geneids) if geneid is not None]
    expression = adata.X[:, keep]
    expression = expression.toarray() if hasattr(expression, "toarray") else expression
    pd.DataFrame(
        expression, columns=[geneids[index] for index in keep]
    ).to_csv(expression_file, index=False)

    meta_file = os.path.join(workdir, "meta.csv")
    pd.DataFrame({"cell": adata.obs_names, "label": "query"}).to_csv(
        meta_file, index=False
    )

    markers_file = os.path.join(workdir, "markers.csv")
    markers_out = markers.rename(columns={"cell_type": "cell type"})
    markers_out["organ"] = args.organ or "Unknown"
    markers_out[["organ", "cell type", "id"]].to_csv(markers_file, index=False)

    run(
        [
            "bash", IMAGE_SCRIPT,
            "-e", expression_file,
            "-m", markers_file,
            "-o", data_dir,
            "-t", "df",
            "-i", meta_file,
            "-f", os.path.join(args.scmapnet_dir, TRANSFER),
            "-d", "0",
            "-s", "test",
            "-n", str(args.ncores),
        ],
        cwd=args.scmapnet_dir,
        description="generating the treemap images",
    )
    test_dir = os.path.join(data_dir, "test")
    if not os.path.isdir(test_dir):
        sys.exit(f"scMapNet generated no test images in {test_dir}")

    run(
        [
            args.python, MAIN,
            "--test",
            "--resume", args.weights,
            "--model", MODEL,
            "--nb_classes", str(len(cell_types)),
            "--data_path", data_dir,
            "--output_dir", out_dir,
        ],
        cwd=args.scmapnet_dir,
        description="predicting the cell types",
    )
    predictions_file = os.path.join(out_dir, "class_pred.csv")
    if not os.path.exists(predictions_file):
        sys.exit(f"scMapNet wrote no predictions to {predictions_file}")

    # ImageFolder walks the classes and the files of a class in sorted order,
    # and the data loader keeps that order, so the rows of class_pred.csv line
    # up with the sorted images; their names carry the (1-based) column index
    # of the cell in the expression file.
    images = [
        filename
        for cls in sorted(os.listdir(test_dir))
        for filename in sorted(os.listdir(os.path.join(test_dir, cls)))
    ]
    indexes = [int(name.rsplit("_", 1)[1].split(".")[0]) for name in images]
    predictions = pd.read_csv(predictions_file)
    if len(predictions) != len(indexes):
        sys.exit(
            f"scMapNet predicted {len(predictions)} of the {len(indexes)} "
            "images; the checkpoint does not match the data"
        )

    barcodes = list(adata.obs_names)
    labels = pd.Series("unassigned", index=barcodes)
    for colidx, pred in zip(indexes, predictions["pred"]):
        cell_type = (
            cell_types[pred] if 0 <= pred < len(cell_types) else "unassigned"
        )
        labels.iloc[colidx - 1] = cell_type
    print(
        f"scMapNet: annotated {int((labels != 'unassigned').sum())} of the "
        f"{len(barcodes)} cells"
    )

    pd.DataFrame(
        {"cell": barcodes, "scmapnet_celltype": labels.to_numpy()}
    ).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
