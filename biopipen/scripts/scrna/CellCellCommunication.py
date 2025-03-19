from pathlib import Path
from biopipen.utils.misc import run_command, logger
import os
import numpy as np
import pandas as pd
import scanpy
import liana
import liana.method.sc._liana_pipe as _liana_pipe

# monkey-patch liana.method.sc._liana_pipe._trimean due to the updates by scipy 1.14
# https://github.com/scipy/scipy/commit/a660202652deead0f3b4b688eb9fdcdf9f74066c
def _trimean(a, axis=0):
    try:
        arr = a.A
    except AttributeError:
        arr = a.toarray()

    quantiles = np.quantile(arr, q=[0.25, 0.75], axis=axis)
    median = np.median(arr, axis=axis)
    return (quantiles[0] + 2 * median + quantiles[1]) / 4


_liana_pipe._trimean = _trimean


sobjfile = Path({{in.sobjfile | quote}})  # pyright: ignore  # noqa: E999
outfile = Path({{out.outfile | quote}})  # pyright: ignore
envs: dict = {{envs | dict}}  # pyright: ignore

# https://github.com/h5py/h5py/issues/1082#issuecomment-1311498466
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
method = envs.pop("method")
assay = envs.pop("assay")
ncores = envs.pop("ncores")
species = envs.pop("species")
rscript = envs.pop("rscript")
subset = envs.pop("subset")
subset_using = envs.pop("subset_using", "auto")
if subset_using == "auto":
    subset_using = "python" if subset and "[" in subset else "r"
split_by = envs.pop("split_by")

if sobjfile.suffix.lower() == ".rds" or sobjfile.suffix.lower() == ".h5seurat":
    logger.info("Converting the Seurat object to h5ad ...")

    annfile = outfile.parent / f"{sobjfile.stem}.h5ad"
    if subset and subset_using == "r":
        r_script_convert_to_anndata = (
            "biopipen.utils::ConvertSeuratToAnnData"
            f"({str(sobjfile)!r}, {str(annfile)!r}, "
            f"assay = {{envs['assay'] | r}}, subset = {{envs['subset'] | r}})"
        )
    else:
        r_script_convert_to_anndata = (
            "biopipen.utils::ConvertSeuratToAnnData"
            f"({str(sobjfile)!r}, {str(annfile)!r}, assay = {{envs['assay'] | r}})"
        )
    run_command([rscript, "-e", r_script_convert_to_anndata], fg=True)
    sobjfile = annfile
elif subset and subset == "r":
    raise ValueError(
        "h5ad file is provided as input, ",
        "'subset' can only be a 'python' expression (`envs.subset_using = 'python'`)."
    )

logger.info("Reading the h5ad file ...")
adata = scanpy.read_h5ad(sobjfile)

if subset and subset_using == "python":
    logger.info("Subsetting the data ...")
    adata = adata[{{envs['subset']}}]  # pyright: ignore

method = method.lower()
if method == "log2fc":
    method_fun = liana.mt.logfc
else:
    method_fun = getattr(liana.mt, method)

envs["resource_name"] = "consensus" if species == "human" else "mouseconsensus"
envs["n_jobs"] = ncores
envs["inplace"] = True
envs["verbose"] = True
envs["key_added"] = "liana_ccc"

if split_by:
    split_vals = adata.obs[split_by].unique()
    result: pd.DataFrame = None  # type: ignore
    for split_val in split_vals:
        logger.info(f"Running {method} for {split_by} = {split_val} ...")
        adata_split = adata[adata.obs[split_by] == split_val]
        envs["adata"] = adata_split

        method_fun(**envs)
        res = adata_split.uns['liana_ccc']
        res[split_by] = split_val

        if result is None:
            result = res
        else:
            result = pd.concat([result, res], ignore_index=True)
else:
    logger.info(f"Running {method} ...")
    envs["adata"] = adata
    method_fun(**envs)

    result = adata.uns['liana_ccc']

mag_score_names = {
    "cellphonedb": "lr_means",
    "connectome": "expr_prod",
    "log2fc": None,
    "natmi": "expr_prod",
    "singlecellsignalr": "lrscore",
    "rank_aggregation": "magnitude_rank",
    "geometric_mean": "lr_gmeans",
    "scseqcomm": "inter_score",
    "cellchat": "lr_probs",
}

spec_score_names = {
    "cellphonedb": "cellphone_pvals",
    "connectome": "scaled_weight",
    "log2fc": "lr_logfc",
    "natmi": "spec_weight",
    "singlecellsignalr": None,
    "rank_aggregation": "specificity_rank",
    "geometric_mean": "gmean_pvals",
    "scseqcomm": None,
    "cellchat": "cellchat_pvals",
}

if mag_score_names[method] is not None:
    result['mag_score'] = result[mag_score_names[method]]
if spec_score_names[method] is not None:
    result['spec_score'] = result[spec_score_names[method]]

logger.info("Saving the result ...")
result.to_csv(outfile, sep="\t", index=False)
