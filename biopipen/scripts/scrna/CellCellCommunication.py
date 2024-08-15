from pathlib import Path
from biopipen.utils.misc import run_command, logger
import numpy as np
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


sobjfile = Path({{in.sobjfile | repr}})  # pyright: ignore  # noqa: E999
outfile = Path({{out.outfile | repr}})  # pyright: ignore
envs = {{envs | repr}}  # pyright: ignore

method = envs.pop("method")
assay = envs.pop("assay")
ncores = envs.pop("ncores")
species = envs.pop("species")
rscript = envs.pop("rscript")

if sobjfile.suffix.lower() == ".rds" or sobjfile.suffix.lower() == ".h5seurat":
    annfile = outfile.parent / f"{sobjfile.stem}.h5ad"
    r_script_convert_to_anndata = f"""
    {{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
    {{ biopipen_dir | joinpaths: "utils", "single_cell.R" | source_r }}

    seurat_to_anndata(
        "{sobjfile}",
        "{annfile}",
        assay = {{ envs.assay | r }},
        log_info = log_info
    )
    """
    run_command([rscript, "-e", r_script_convert_to_anndata], fg=True)

    sobjfile = annfile

logger.info("Reading the h5ad file ...")
adata = scanpy.read_h5ad(sobjfile)

method = method.lower()
if method == "log2fc":
    method_fun = liana.mt.logfc
else:
    method_fun = getattr(liana.mt, method)

logger.info(f"Running {method} ...")
envs["adata"] = adata
envs["resource_name"] = "consensus" if species == "human" else "mouseconsensus"
envs["n_jobs"] = ncores
envs["inplace"] = True
envs["verbose"] = True
envs["key_added"] = "liana_ccc"
method_fun(**envs)

res = adata.uns['liana_ccc']

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
    res['mag_score'] = res[mag_score_names[method]]
if spec_score_names[method] is not None:
    res['spec_score'] = res[spec_score_names[method]]

logger.info("Saving the result ...")
res.to_csv(outfile, sep="\t", index=False)
