from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from multiprocessing import Pool

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from scipy.io import mmread
from vireoSNP import BinomMixtureVB
from vireoSNP.plot import heat_matrix, anno_heat
from vireoSNP.plot.base_plot import vireo_colors
from biopipen.utils.misc import require_package, logger

require_package("vireoSNP", ">=0.5.8")

crdir = Path({{in.cellsnpout | quote}})  # noqa: E999 # pyright: ignore
outdir: str = {{out.outdir | quote}}  # pyright: ignore
ncores : int = {{envs.ncores | int}}  # pyright: ignore
n_init: int = {{envs.n_init | int}}  # pyright: ignore
n_clones: int | None = {{envs.n_clones | repr}}  # pyright: ignore
min_iter: int = {{envs.min_iter | int}}  # pyright: ignore
max_iter: int = {{envs.max_iter | int}}  # pyright: ignore
seed: int = {{envs.seed | int}}  # pyright: ignore

ad_file = crdir / "passed_ad.mtx"
if not ad_file.exists():
    ad_file = crdir / "cellSNP.tag.AD.mtx"
    if not ad_file.exists():
        raise FileNotFoundError(
            f"AD matrix file not found in {crdir}, "
            "expected at passed_ad.mtx or cellSNP.tag.AD.mtx"
        )
    mquad = False
else:
    mquad = True

if mquad:
    dp_file = crdir / "passed_dp.mtx"
    mtSNP_ids = np.genfromtxt(crdir / "passed_variant_names.txt", dtype='str')
else:
    dp_file = crdir / "cellSNP.tag.DP.mtx"
    mtSNP_ids = None

logger.info("Reading and converting matrices to CSC format...")
ad_data = mmread(ad_file).tocsc()
dp_data = mmread(dp_file).tocsc()

if ad_data.shape[0] < 2:
    logger.warning(
        f"AD matrix has only {ad_data.shape[0]} variants, "
        "VireoSNP requires at least 2 variants to run. "
    )
    logger.warning("Generating empty output files and exiting...")
    Path(f"{outdir}/cell_clone_assignment.tsv").write_text("cell_id\tassigned_clone\tassignment_prob\n")
    Path(f"{outdir}/best_n_clones.txt").write_text("0\n")
    plt.figure()
    # Empty ELBO plot, with text: Only {ad_data.shape[0]} variants found
    plt.text(0.5, 0.5, f"Only {ad_data.shape[0]} variants found,\nVireoSNP requires at least 2 variants to run.",
             horizontalalignment='center', verticalalignment='center', fontsize=12)
    plt.savefig(f"{outdir}/ELBO_n_clones.png")

    plt.figure()
    # Empty model fitting plot
    plt.text(0.5, 0.5, "No model fitting performed.",
             horizontalalignment='center', verticalalignment='center', fontsize=12)
    plt.savefig(f"{outdir}/model_fitting.png")

    plt.figure()
    # Empty heatmap plot
    plt.text(0.5, 0.5, "No heatmap generated.",
             horizontalalignment='center', verticalalignment='center', fontsize=12)
    plt.savefig(f"{outdir}/clone_allele_heatmap.png")
    exit(0)

# Check for cell_number_info.txt to get sample information
cell_number_info_file = crdir / "cell_number_info.txt"
if cell_number_info_file.exists():
    logger.info("Found cell_number_info.txt, loading sample information...")
    cell_info_df = pd.read_csv(cell_number_info_file, sep="\t")
    # Create cell-to-sample mapping based on cumulative cell counts
    cell_to_sample = {}
    cell_ids = []
    col_offset = 0
    for _, row in cell_info_df.iterrows():
        sample = row["Sample"]
        n_cells = row["CellNumber"]
        # Map cells to samples
        for i in range(n_cells):
            cell_to_sample[col_offset + i] = sample
            cell_ids.append(f"{sample}_cell{i}")
        col_offset += n_cells
    logger.info(f"Loaded sample information for {len(cell_info_df)} samples, {col_offset} total cells")
else:
    logger.info("No cell_number_info.txt found, proceeding without sample information")
    cell_to_sample = None
    cell_ids = None


def fit_model_for_k(k):
    """Fit BinomMixtureVB model for a given number of clones."""
    logger.info(f"Fitting model with {k} clones...")
    _model = BinomMixtureVB(
        n_var=ad_data.shape[0],
        n_cell=ad_data.shape[1],
        n_donor=k,
    )
    _model.fit(
        ad_data,
        dp_data,
        n_init=n_init,
        min_iter=min_iter,
        max_iter=max_iter,
        random_seed=seed,
    )
    return _model.ELBO_inits


if isinstance(n_clones, (tuple, list)):
    if len(n_clones) != 2:
        raise ValueError(
            "n_clones should be an integer or a list/tuple of two integers for range"
        )

    n_clone_list = np.arange(n_clones[0], n_clones[1] + 1)
    logger.info(f"Testing n_clones from {n_clones[0]} to {n_clones[1]} using {ncores} cores...")

    with Pool(processes=ncores) as pool:
        _ELBO_mat = pool.map(fit_model_for_k, n_clone_list)

    # Select the best K based on the elbow
    mean_ELBOs = [np.mean(elbos) for elbos in _ELBO_mat]
    deltas = np.diff(mean_ELBOs)
    second_deltas = np.diff(deltas)
    best_k_index = np.argmax(second_deltas) + 1  # +1 due to diff reducing length by 1
    best_k = n_clones[0] + best_k_index
    logger.info(f"Selected best n_clones: {best_k}, based on ELBO elbow method.")
    Path(f"{outdir}/best_n_clones.txt").write_text(str(best_k))

    # Visualize the ELBOs
    plt.plot(np.arange(1, len(n_clone_list)+1), np.max(_ELBO_mat, axis=1))
    plt.boxplot(_ELBO_mat)
    plt.xticks(np.arange(1, len(n_clone_list)+1), n_clone_list)
    plt.ylabel("ELBO")
    plt.xlabel("n_clones")

    plt.savefig(f"{outdir}/ELBO_n_clones.png")
else:
    best_k = n_clones

logger.info(f"Fitting final model with {best_k} clones...")
_model = BinomMixtureVB(n_var=ad_data.shape[0], n_cell=ad_data.shape[1], n_donor=best_k)
_model.fit(
    ad_data,
    dp_data,
    n_init=n_init,
    min_iter=min_iter,
    max_iter=max_iter,
    random_seed=seed,
)

# model fitting
logger.info("Generating model fitting plots...")
fig = plt.figure(figsize=(11, 4))
plt.subplot(1, 2, 1)
plt.hist(_model.ELBO_inits)
plt.ylabel("Frequency")
plt.xlabel("ELBO in multiple initializations")

plt.subplot(1, 2, 2)
plt.plot(_model.ELBO_iters)
plt.xlabel("Iterations")
plt.ylabel("ELBO in a single initialization")

plt.tight_layout()
plt.savefig(f"{outdir}/model_fitting.png")

# Visualize assignment probability and allele frequency
raw_col = matplotlib.colormaps.get_cmap('pink_r').resampled(200)
new_col = np.vstack(
    (raw_col(np.linspace(0, 0.7, 10)), raw_col(np.linspace(0.7, 1, 90)))
)
segpink = ListedColormap(new_col.tolist(), name='segpink')

fig = plt.figure(figsize=(7, 4), dpi=100)
plt.subplot(1, 2, 1)
im = heat_matrix(_model.ID_prob, cmap="Blues", alpha=0.8,
                 display_value=False, row_sort=True)
plt.colorbar(im, fraction=0.046, pad=0.04)
plt.title("Assignment probability")
plt.xlabel("Clone")
plt.ylabel("%d cells" %(_model.n_cell))
plt.xticks(range(_model.n_donor))


plt.subplot(1, 2, 2)
im = heat_matrix(_model.beta_mu, cmap=segpink, alpha=0.8,
                 display_value=False, row_sort=True)
plt.colorbar(im, fraction=0.046, pad=0.04)
plt.title("Mean allelic ratio")
plt.xlabel("Clone")
plt.ylabel("%d SNPs" %(_model.n_var))
plt.xticks(range(_model.n_donor))

plt.tight_layout()

plt.savefig(f"{outdir}/assignment_allele_freq.png")

# Visualize clones
logger.info("Preparing clone visualization...")
if mtSNP_ids is None:
    mtSNP_ids = ['variant%d' %x for x in range(ad_data.shape[0])]

cell_label = np.array(['clone%d' %x for x in np.argmax(_model.ID_prob, axis=1)])
id_uniq = np.unique(cell_label)
clone_ids = np.argmax(_model.ID_prob, axis=1)
af_mat = ad_data.toarray()/(dp_data.toarray() + 0.0001)
var_idx = np.argsort(af_mat @ (5**(max(clone_ids) - clone_ids)))[::-1]

# Create sample labels if available
if cell_to_sample is not None:
    logger.info("Creating sample labels...")
    sample_label = np.array([cell_to_sample[i] for i in range(_model.n_cell)])
    sample_uniq = sorted(list(set(cell_to_sample.values())))
    logger.info("Sample labels created.")
else:
    sample_label = None
    sample_uniq = None


def anno_heat(X, row_anno=None, col_anno=None, col_anno2=None,
              row_order_ids=None, col_order_ids=None, col_order_ids2=None,
              xticklabels=False, yticklabels=False,
              row_cluster=False, col_cluster=False,
              **kwargs):
    """
    Heatmap with column or row annotations. Based on seaborn.clustermap()
    Row or column will be ordered by the annotation group.
    Supports dual column annotations (col_anno and col_anno2).

    Note, haven't tested if input both row_anno and col_anno.
    """

    import seaborn as sns

    # Extend vireo_colors if needed by cycling through the available colors
    def get_color_for_index(idx):
        """Get color for an index, cycling through vireo_colors if needed."""
        return vireo_colors[idx % len(vireo_colors)]

    # prepare row annotation
    if row_anno is not None:
        if row_order_ids is None:
            row_order_ids = list(np.unique(row_anno))
        else:
            row_order_ids = [x for x in row_order_ids]
        row_num = np.array([row_order_ids.index(x) for x in row_anno])

        # dot_row = np.array(np.nansum(X, axis=1)).reshape(-1)
        idx_row = np.argsort(row_num * 2**X.shape[1])# + dot_row / dot_row.max())

        # Create color array for each row
        row_colors = np.array([get_color_for_index(num) for num in row_num])[idx_row]
    else:
        row_colors = None
        row_order_ids = []
        idx_row = range(X.shape[0])

    # prepare col annotation
    if col_anno is not None:
        if col_order_ids is None:
            col_order_ids = list(np.unique(col_anno))
        else:
            col_order_ids = [x for x in col_order_ids]
        col_num = np.array([col_order_ids.index(x) for x in col_anno])

        # dot_col = np.array(np.nansum(X, axis=0)).reshape(-1)
        # Python int too large to convert to C long
        # idx_col = np.argsort(col_num * 2**X.shape[0])# + dot_row / dot_row.max())
        idx_col = np.argsort(col_num * (10**6))  # + dot_row / dot_row.max())
        # Create color array for each column
        col_colors = np.array([get_color_for_index(num) for num in col_num])[idx_col]

        # Add second column annotation if provided
        if col_anno2 is not None:
            if col_order_ids2 is None:
                col_order_ids2 = list(np.unique(col_anno2))
            else:
                col_order_ids2 = [x for x in col_order_ids2]
            col_num2 = np.array([col_order_ids2.index(x) for x in col_anno2])
            col_colors2 = np.array([get_color_for_index(num + len(col_order_ids)) for num in col_num2])[idx_col]
            col_colors = [col_colors, col_colors2]
    else:
        col_colors = None
        col_order_ids = []
        col_order_ids2 = []
        idx_col = range(X.shape[1])

    ## plot with seaborn clustermap
    g = sns.clustermap(X[idx_row, :][:, idx_col],
                       row_colors=row_colors, col_colors=col_colors,
                       col_cluster=col_cluster, row_cluster=row_cluster,
                       xticklabels=xticklabels, yticklabels=yticklabels,
                       **kwargs)

    if row_anno is not None:
        for i in range(len(row_order_ids)):
            g.ax_row_dendrogram.bar(0, 0, color=get_color_for_index(i),
                                    label=row_order_ids[i], linewidth=0)
        g.ax_row_dendrogram.legend(loc="center", ncol=1, title="")

    if col_anno is not None:
        # Create first legend for clones
        clone_handles = []
        for i in range(len(col_order_ids)):
            handle = g.ax_col_dendrogram.bar(0, 0, color=get_color_for_index(i),
                                    label=col_order_ids[i], linewidth=0)
            clone_handles.append(handle)

        # Create second legend for samples if provided
        if col_anno2 is not None:
            sample_handles = []
            for i in range(len(col_order_ids2)):
                handle = g.ax_col_dendrogram.bar(0, 0, color=get_color_for_index(i + len(col_order_ids)),
                                        label=col_order_ids2[i], linewidth=0)
                sample_handles.append(handle)

            # Create two separate legends
            legend1 = g.ax_col_dendrogram.legend(handles=clone_handles, loc="upper center",
                                                 ncol=min(len(col_order_ids), 6),
                                                 title="Clones", frameon=True)
            g.ax_col_dendrogram.add_artist(legend1)  # Add first legend back
            g.ax_col_dendrogram.legend(handles=sample_handles, loc="lower center",
                                      ncol=min(len(col_order_ids2), 6),
                                      title="Samples", frameon=True)
        else:
            g.ax_col_dendrogram.legend(loc="center", ncol=6, title="Clones")

    g.cax.set_position([1.01, .2, .03, .45])

    return g


logger.info("Generating clone allele heatmap with dual annotations...")
# Calculate figure size based on number of cells and variants
n_cells = ad_data.shape[1]
n_vars = len(var_idx)
# Width: scale with number of cells (min 10, max 50 inches)
fig_width = max(10, min(50, n_cells / 200))
# Height: scale with number of variants (min 8, max 40 inches)
fig_height = max(8, min(40, n_vars / 20))
logger.info(f"Heatmap dimensions: {fig_width:.1f}x{fig_height:.1f} inches ({n_cells} cells x {n_vars} variants)")

if sample_label is not None:
    im = anno_heat(
        (ad_data/dp_data)[var_idx, :],
        col_anno=cell_label,
        col_anno2=sample_label,
        col_order_ids=id_uniq,
        col_order_ids2=sample_uniq,
        cmap=segpink,
        yticklabels=[mtSNP_ids[idx] for idx in var_idx],  # type: ignore
        figsize=(fig_width, fig_height),
    )
else:
    im = anno_heat(
        (ad_data/dp_data)[var_idx, :],
        col_anno=cell_label,
        col_order_ids=id_uniq,
        cmap=segpink,
        yticklabels=[mtSNP_ids[idx] for idx in var_idx],  # type: ignore
        figsize=(fig_width, fig_height),
    )

plt.savefig(f"{outdir}/clone_allele_heatmap.png", dpi=100, bbox_inches='tight')
logger.info("Heatmap saved.")

# Save assignments
logger.info("Saving cell-clone assignments...")
if cell_ids is not None and cell_to_sample is not None:
    # With sample information
    df_assign = pd.DataFrame({
        "cell_id": cell_ids,
        "sample": [cell_to_sample[i] for i in range(_model.n_cell)],
        "assigned_clone": np.argmax(_model.ID_prob, axis=1),
        "assignment_prob": np.max(_model.ID_prob, axis=1),
    })
else:
    # Without sample information
    df_assign = pd.DataFrame({
        "cell_id": [f"cell{x}" for x in range(_model.n_cell)],
        "assigned_clone": np.argmax(_model.ID_prob, axis=1),
        "assignment_prob": np.max(_model.ID_prob, axis=1),
    })
df_assign.to_csv(f"{outdir}/cell_clone_assignment.tsv", sep="\t", index=False)

logger.info("VireoSNP analysis completed successfully.")
