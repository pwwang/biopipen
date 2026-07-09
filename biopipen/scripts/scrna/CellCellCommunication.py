from diot import Diot  # noqa
from pathlib import Path
from biopipen.utils.misc import logger
from biopipen.scripts.scrna.seurat_anndata_conversion import convert_seurat_to_anndata
import os
import numpy as np
import pandas as pd
import scanpy
import liana
import liana.method.sc._liana_pipe as _liana_pipe

# AttributeError: module 'numpy' has no attribute 'product'
if not hasattr(np, "product"):
    np.product = np.prod

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
assay = envs.pop("assay")
species = envs.pop("species")
rscript = envs.pop("rscript")
subset_using = envs.pop("subset_using", None)
subset = envs.pop("subset")
if subset and subset_using not in ["python", "r", "R"]:
    raise ValueError(
        f"Invalid value for `subset_using`: {subset_using}. "
        "It must be either 'python' or 'r' (or 'R')."
    )

cases = envs.pop("cases", {})

if not cases:
    cases = {"DEFAULT": envs.copy()}
else:
    for case_name, case_envs in cases.items():
        if not isinstance(case_envs, dict):
            raise ValueError(f"Case '{case_name}' must be a dictionary of envs.")

        if any(key in ["method", "species"] for key in case_envs.keys()):
            raise ValueError(
                f"Case '{case_name}' cannot override 'method' or 'species'."
            )

        tmp_envs = envs.copy()
        tmp_envs.update(case_envs)

        cases[case_name] = tmp_envs

DEFAULT_ONLY = "DEFAULT" in cases and len(cases) == 1

seurat_ident_col = None
if sobjfile.suffix.lower() in (".rds", ".qs", "qs2"):
    annfile = outfile.parent / f"{sobjfile.stem}.h5ad"
    logger.info("Converting the Seurat object to AnnData ...")
    seurat_ident_col = convert_seurat_to_anndata(
        input_file=str(sobjfile),
        output_file=str(annfile),
        assay=assay,
        subset=subset if subset_using == "r" else None,
        rscript=rscript,
        return_ident_col=True,
    )
    # groupby = groupby or seurat_ident_col
    sobjfile = annfile
elif subset and subset_using in ["r", "R"]:
    raise ValueError(
        "h5ad file is provided as input, "
        "'subset' can only be a 'python' expression (`envs.subset_using = 'python'`)."
    )

logger.info("Reading the h5ad file ...")
adata = scanpy.read_h5ad(sobjfile)

if subset and subset_using == "python":
    logger.info("Subsetting the data ...")
    adata = adata[adata.obs.query(subset).index]


def do_case(name):
    logger.info(f"- Running case '{name}' ...")
    case = cases[name]

    group_by = case.pop("group_by", None)
    case["groupby"] = case.pop("groupby", None) or group_by
    case_subset_using = case.pop("subset_using", None)
    case_subset = case.pop("subset", None)
    if case_subset and case_subset_using != "python":
        raise ValueError(
            f"Case '{name}' has 'subset' specified, but 'subset_using' is not 'python'. "
            "Only 'subset_using = \"python\"' is supported for the cases."
        )

    if case_subset:
        logger.info(f"  Subsetting the data for case '{name}' ...")
        case_adata = adata[adata.obs.query(case_subset).index]
    else:
        case_adata = adata

    if not case["groupby"]:
        if seurat_ident_col:
            logger.warning(
                "  `groupby` is not provided. "
                f"Using '{seurat_ident_col}' as the default groupby column."
            )
            case["groupby"] = seurat_ident_col
        else:
            logger.warning(
                "  `groupby` is not provided and no identity column can be determined from the input file. "
                "Using 'seurat_clusters' as the default groupby column, but it is recommended to provide a valid `groupby` column."
            )
            case["groupby"] = "seurat_clusters"

    method = case.pop("method").lower()
    if method == "log2fc":
        method_fun = liana.mt.logfc
    else:
        method_fun = getattr(liana.mt, method)

    case["resource_name"] = "consensus" if species == "human" else "mouseconsensus"
    case["n_jobs"] = case.pop("ncores", 1)
    case["inplace"] = True
    case["verbose"] = True
    case["key_added"] = "liana_ccc"

    split_by = case.pop("split_by", None)
    if split_by:
        split_vals = case_adata.obs[split_by].dropna().unique()
        result: pd.DataFrame = None  # type: ignore
        for split_val in split_vals:
            logger.info(f"  Running {method} for {split_by} = {split_val} ...")
            adata_split = case_adata[case_adata.obs[split_by] == split_val]
            # case_adata = adata_split

            case["adata"] = adata_split
            method_fun(**case)
            res = adata_split.uns['liana_ccc']
            res[split_by] = split_val

            if result is None:
                result = res
            else:
                result = pd.concat([result, res], ignore_index=True)
    else:
        logger.info(f"  Running {method} ...")
        case["adata"] = case_adata
        method_fun(**case)

        result = case_adata.uns['liana_ccc']

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

    if not DEFAULT_ONLY:
        result.insert(0, "Case", name)

    return result


final_result = pd.concat([do_case(name) for name in cases], ignore_index=True)

logger.info("Saving the final result ...")
final_result.to_csv(outfile, sep="\t", index=False)
