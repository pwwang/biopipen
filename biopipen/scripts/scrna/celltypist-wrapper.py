from argparse import ArgumentParser
from typing import Union
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist.classifier import logger, AnnData, Model, Classifier

parser = ArgumentParser(description="Run CellTypist")
parser.add_argument(
    "-i", "--input", required=True, help="Input H5AD file with AnnData object"
)
parser.add_argument("-o", "--output", required=True, help="Output file")
parser.add_argument("-m", "--model", required=True, help="Model file")
parser.add_argument(
    "-v", "--majority_voting", action="store_true", help="Majority voting"
)
parser.add_argument(
    "-c",
    "--over_clustering",
    required=False,
    default=None,
    help="Over clustering. Error if the column does not exist.",
)


def classifier_init(
    self, filename="", model="", transpose=False, gene_file=None, cell_file=None
):
    """Celltypist check if adata is in the range of log1p normalized data to 10000
    counts per cell. Otherwise it will use the raw data if available. However, in
    some cases, the raw data has invalid feature names (var_names) which causes errors.
    Here we check if the feature names of raw data is valid with intersection with
    model features, if not, we will use the adata.X instead of adata.raw.X
    """
    if isinstance(model, str):
        model = Model.load(model)
    self.model = model
    if not filename:
        logger.warn("📭 No input file provided to the classifier")
        return
    if isinstance(filename, str):
        self.filename = filename
        logger.info(f"📁 Input file is '{self.filename}'")
        logger.info("⏳ Loading data")
    if isinstance(filename, str) and filename.endswith(
        (".csv", ".txt", ".tsv", ".tab", ".mtx", ".mtx.gz")
    ):
        self.adata = sc.read(self.filename)
        if transpose:
            self.adata = self.adata.transpose()
        if self.filename.endswith((".mtx", ".mtx.gz")):
            if (gene_file is None) or (cell_file is None):
                raise FileNotFoundError(
                    "🛑 Missing `gene_file` and/or `cell_file`. Please provide both "
                    "arguments together with the input mtx file"
                )
            genes_mtx = pd.read_csv(gene_file, header=None)[0].values
            cells_mtx = pd.read_csv(cell_file, header=None)[0].values
            if len(genes_mtx) != self.adata.n_vars:
                raise ValueError(
                    f"🛑 The number of genes in {gene_file} does not match the number "
                    f"of genes in {self.filename}"
                )
            if len(cells_mtx) != self.adata.n_obs:
                raise ValueError(
                    f"🛑 The number of cells in {cell_file} does not match the number "
                    f"of cells in {self.filename}"
                )
            self.adata.var_names = genes_mtx
            self.adata.obs_names = cells_mtx
        if not float(self.adata.X[:1000].max()).is_integer():
            logger.warn(
                "⚠️ Warning: the input file seems not a raw count matrix. The "
                "prediction result may not be accurate"
            )
        if (
            (self.adata.n_vars >= 100000)
            or (len(self.adata.var_names[0]) >= 30)
            or (
                len(
                    self.adata.obs_names.intersection(
                        ["GAPDH", "ACTB", "CALM1", "PTPRC", "MALAT1"]
                    )
                )
                >= 1
            )
        ):
            logger.warn(
                "⚠️ The input matrix is detected to be a gene-by-cell matrix, will "
                "transpose it"
            )
            self.adata = self.adata.transpose()
        self.adata.var_names_make_unique()
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        self.indata = self.adata.X
        self.indata_genes = self.adata.var_names
        self.indata_names = self.adata.obs_names
    elif isinstance(filename, AnnData) or (
        isinstance(filename, str) and filename.endswith(".h5ad")
    ):
        self.adata = sc.read(filename) if isinstance(filename, str) else filename
        self.adata.var_names_make_unique()
        # When to use raw.X?
        # 1. if adata.raw exists
        # 2. if adata.raw.var_names has intersection with model genes
        # 3. if adata.X is not in the expected range
        use_raw = self.adata.raw and (
            self.adata.X[:1000].min() < 0 or self.adata.X[:1000].max() > 9.22
        ) and np.isin(
            self.adata.raw.var_names, self.model.classifier.features
        ).sum() > 0

        if use_raw:
            if not self.adata.raw:
                raise ValueError(
                    "🛑 Invalid expression matrix in `.X`, expect log1p normalized "
                    "expression to 10000 counts per cell"
                )
            elif (self.adata.raw.X[:1000].min() < 0) or (
                self.adata.raw.X[:1000].max() > 9.22
            ):
                raise ValueError(
                    "🛑 Invalid expression matrix in both `.X` and `.raw.X`, expect "
                    "log1p normalized expression to 10000 counts per cell"
                )
            else:
                logger.info(
                    "👀 Invalid expression matrix in `.X`, expect log1p normalized "
                    "expression to 10000 counts per cell; will use `.raw.X` instead"
                )
                self.indata = self.adata.raw.X
                self.indata_genes = self.adata.raw.var_names
                self.indata_names = self.adata.raw.obs_names
        else:
            self.indata = self.adata.X
            self.indata_genes = self.adata.var_names
            self.indata_names = self.adata.obs_names
        if np.abs(np.expm1(self.indata[0]).sum() - 10000) > 1:
            logger.warn(
                "⚠️ Warning: invalid expression matrix, expect ALL genes and log1p "
                "normalized expression to 10000 counts per cell. The prediction result "
                "may not be accurate"
            )
    else:
        raise ValueError(
            "🛑 Invalid input. Supported types: .csv, .txt, .tsv, .tab, .mtx, .mtx.gz "
            "and .h5ad, or AnnData loaded in memory"
        )

    logger.info(
        f"🔬 Input data has {self.indata.shape[0]} cells and {len(self.indata_genes)} "
        "genes"
    )


if __name__ == "__main__":
    Classifier.__init__ = classifier_init  # type: ignore

    args = parser.parse_args()
    adata = sc.read_h5ad(args.input)
    over_clustering = args.over_clustering
    if over_clustering and over_clustering not in adata.obs.columns:
        raise ValueError(
            f"Over clustering column '{over_clustering}' not found in AnnData object."
        )
    if "neighbors" in adata.uns and "params" in adata.uns["neighbors"]:
        adata.uns["neighbors"]["params"].setdefault("n_neighbors", 15)

    annotated = celltypist.annotate(
        adata,
        model=args.model,
        majority_voting=args.majority_voting,
        over_clustering=over_clustering,
    )

    out_adata = annotated.to_adata()
    # leave as is
    # if over_clustering and args.majority_voting:
    #     # rename majority_voting column to over_clustering
    #     out_adata.obs[over_clustering] = out_adata.obs["majority_voting"]

    if args.output.endswith(".h5ad"):
        try:
            out_adata._raw._var.rename(  # type: ignore
                columns={"_index": "features"}, inplace=True
            )
            del out_adata.raw
        except (KeyError, AttributeError):
            pass

        out_adata.write(args.output)
    else:
        out_adata.obs.to_csv(args.output, sep="\t", index=True)
