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

from importlib.metadata import version as _version
from packaging.version import parse as _v

# Monkey-patch anndata.AnnData.__init__ to handle `dtype` keyword on anndata >= 0.13.
# anndata 0.13 removed the `dtype` parameter from AnnData.__init__, but liana
# still passes `dtype="float32"` as a keyword argument. The legacy_api_wrap
# decorator only catches `dtype` when passed positionally, not as a keyword,
# so it leaks through and causes TypeError.
if _v(_version("anndata")) >= _v("0.13.0"):
    def _patch_anndata_dtype():
        import anndata
        if not hasattr(anndata.AnnData, "_immunopipe_patched"):
            _orig_init = anndata.AnnData.__init__

            def _patched_init(self, *args, **kwargs):
                import numpy as np
                dtype = kwargs.pop("dtype", None)
                if dtype is not None:
                    X = kwargs.get("X", args[0] if args else None)
                    if X is not None:
                        try:
                            from scipy import sparse
                        except ImportError:
                            sparse = None
                        if sparse is not None and (
                            sparse.issparse(X) or isinstance(X, np.ndarray)
                        ):
                            if X.dtype != np.dtype(dtype):
                                kwargs["X"] = X.astype(dtype)
                        elif isinstance(X, np.ndarray):
                            if X.dtype != np.dtype(dtype):
                                kwargs["X"] = X.astype(dtype)
                        else:
                            kwargs["X"] = np.asarray(X, dtype)
                return _orig_init(self, *args, **kwargs)

            anndata.AnnData.__init__ = _patched_init
            anndata.AnnData._immunopipe_patched = True

    _patch_anndata_dtype()

# patch liana to avoid the error:

# AttributeError: module 'numpy' has no attribute 'product'
if not hasattr(np, "product"):
    np.product = np.prod

# monkey-patch liana.method.sc._liana_pipe._trimean due to the updates by scipy 1.14
# https://github.com/scipy/scipy/commit/a660202652deead0f3b4b688eb9fdcdf9f74066c
# Also patch _get_lr to force mat_max to float32 — in scipy >= 1.14, sparse matrix
# .max() returns np.float64 instead of np.float32, breaking the isinstance(mat_max,
# np.float32) gate that triggers cellchat-specific trimean calculation.
if _v(_version("scipy")) >= _v("1.14.0"):
    def _trimean(a, axis=0):
        try:
            arr = a.A
        except AttributeError:
            arr = a.toarray()

        quantiles = np.quantile(arr, q=[0.25, 0.75], axis=axis)
        median = np.median(arr, axis=axis)
        return (quantiles[0] + 2 * median + quantiles[1]) / 4

    _liana_pipe._trimean = _trimean

    _original_get_lr = _liana_pipe._get_lr

    def _patched_get_lr(
        adata, resource, groupby_pairs, relevant_cols,
        mat_mean, mat_max, de_method, base, verbose,
    ):
        if mat_max is not None:
            mat_max = np.float32(mat_max)
        return _original_get_lr(
            adata, resource, groupby_pairs, relevant_cols,
            mat_mean, mat_max, de_method, base, verbose,
        )

    _liana_pipe._get_lr = _patched_get_lr


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
ALL_SPLIT_BY_COLS = []
for case_envs in cases.values():
    if isinstance(case_envs, dict) and case_envs.get("split_by"):
        split_by = case_envs["split_by"]
        if isinstance(split_by, str):
            split_by = [split_by]
        ALL_SPLIT_BY_COLS.extend(split_by)
ALL_SPLIT_BY_COLS = list(dict.fromkeys(ALL_SPLIT_BY_COLS))

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
        if isinstance(split_by, str):
            split_by = [split_by]

        # Keep only cells with all split columns non-NA (matches previous behavior)
        split_obs = case_adata.obs[split_by].dropna()

        # Pick a separator that no value contains, so the composite keys can
        # be separated back into the original column values unambiguously.
        sep = "_"
        all_vals = set()
        for col in split_by:
            all_vals.update(split_obs[col].astype(str).unique())
        while any(sep in v for v in all_vals):
            sep += "_"

        split_key = split_obs[split_by[0]].astype(str)
        for col in split_by[1:]:
            split_key = split_key + sep + split_obs[col].astype(str)

        result: pd.DataFrame = None  # type: ignore
        for split_val in split_key.unique():
            vals = split_val.split(sep)
            logger.info(f"  Running {method} for {split_by} = {vals} ...")
            # split_key only covers the dropna'd rows; align the mask with the obs index
            mask = (split_key == split_val).reindex(case_adata.obs.index, fill_value=False)
            adata_split = case_adata[mask]

            case["adata"] = adata_split
            method_fun(**case)
            res = adata_split.uns['liana_ccc']

            # Separate the key back into the original values of each column
            orig_vals = split_obs.loc[split_key == split_val].iloc[0]
            for col in reversed(split_by):
                res.insert(0, col, orig_vals[col])  # insert at the first column

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

    for split_col in ALL_SPLIT_BY_COLS:
        if split_col not in result.columns:
            # result[split_col] = "NA"
            result.insert(0, split_col, "NA")  # insert at the first column

    if not DEFAULT_ONLY:
        result.insert(0, "Case", name)

    return result


final_result = pd.concat([do_case(name) for name in cases], ignore_index=True)

logger.info("Saving the final result ...")
final_result.to_csv(outfile, sep="\t", index=False)
