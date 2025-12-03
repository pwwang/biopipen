from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.io import mmread
from vireoSNP import BinomMixtureVB
from vireoSNP.plot import heat_matrix, anno_heat
from biopipen.utils.misc import require_package

require_package("vireoSNP", ">=0.5.8")

crdir = Path({{in.cellsnpout | quote}})  # noqa: E999 # pyright: ignore
outdir = {{out.outdir | quote}}  # pyright: ignore
n_init = {{envs.n_init | int}}  # pyright: ignore
n_clones = {{envs.n_clones | repr}}  # pyright: ignore
min_iter = {{envs.min_iter | int}}  # pyright: ignore
max_iter = {{envs.max_iter | int}}  # pyright: ignore
seed = {{envs.seed | int}}  # pyright: ignore

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

ad_data = mmread(ad_file).tocsc()
dp_data = mmread(dp_file).tocsc()

if isinstance(n_clones, (tuple, list)):
    if len(n_clones) != 2:
        raise ValueError(
            "n_clones should be an integer or a list/tuple of two integers for range"
        )

    _ELBO_mat = []
    n_clone_list = np.arange(n_clones[0], n_clones[1] + 1)
    for k in n_clone_list:
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
        _ELBO_mat.append(_model.ELBO_inits)

    # Select the best K based on the elbow
    mean_ELBOs = [np.mean(elbos) for elbos in _ELBO_mat]
    deltas = np.diff(mean_ELBOs)
    second_deltas = np.diff(deltas)
    best_k_index = np.argmax(second_deltas) + 1  # +1 due to diff reducing length by 1
    best_k = n_clones[0] + best_k_index
    print(f"Selected best n_clones: {best_k}, based on ELBO elbow method.")
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

if mtSNP_ids is not None:
    mtSNP_ids = ['mt_variant%d' %x for x in range(ad_data.shape[0])]
else:
    mtSNP_ids = ['variant%d' %x for x in range(ad_data.shape[0])]

cell_label = np.array(['clone%d' %x for x in np.argmax(_model.ID_prob, axis=1)])
id_uniq = np.unique(cell_label)
clone_ids = np.argmax(_model.ID_prob, axis=1)
af_mat = ad_data.toarray()/(dp_data.toarray() + 0.0001)
var_idx = np.argsort(af_mat @ (5**(max(clone_ids) - clone_ids)))[::-1]

im = anno_heat(
    (ad_data/dp_data)[var_idx, :], col_anno=cell_label,
    col_order_ids=id_uniq,
    cmap=segpink,
    yticklabels=[mtSNP_ids[idx] for idx in var_idx],
)

plt.savefig(f"{outdir}/clone_allele_heatmap.png")

# Save assignments
df_assign = pd.DataFrame({
    "cell_id": [f"cell{x}" for x in range(_model.n_cell)],
    "assigned_clone": np.argmax(_model.ID_prob, axis=1),
    "assignment_prob": np.max(_model.ID_prob, axis=1),
})
df_assign.to_csv(f"{outdir}/cell_clone_assignment.tsv", sep="\t", index=False)
print("VireoSNP analysis completed successfully.")
