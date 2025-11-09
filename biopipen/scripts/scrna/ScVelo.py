from __future__ import annotations
import os
import warnings
from pathlib import Path

from diot import Diot  # type: ignore[import]
import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from biopipen.utils.misc import logger, require_package
from biopipen.scripts.scrna.seurat_anndata_conversion import (
    convert_seurat_to_anndata,
    convert_anndata_to_seurat,
)

require_package("scvelo", ">=0.3.3")
from biopipen.scripts.scrna import scvelo_paga  # noqa: F401

warnings.simplefilter("ignore", category=UserWarning)
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=DeprecationWarning)


def SCVELO(
    adata,
    group_by,
    dirpath,
    logger,
    palette=None,
    linear_reduction=None,
    nonlinear_reduction=None,
    basis=None,
    mode=["deterministic", "stochastic", "dynamical"],
    fitting_by="stochastic",
    min_shared_counts=30,
    n_pcs=30,
    n_neighbors=30,
    stream_smooth=None,
    stream_density=2,
    arrow_size=5,
    arrow_length=5,
    arrow_density=0.5,
    denoise=False,
    denoise_topn=3,
    kinetics=False,
    kinetics_topn=100,
    calculate_velocity_genes=False,
    top_n=6,
    ncores=1,
    dpi=100,
    fileprefix="",
):
    os.chdir(os.path.expanduser(dirpath))
    if linear_reduction is None:
        sc.pp.pca(adata, n_comps=n_pcs)
        linear_reduction = "X_pca"
    elif linear_reduction not in adata.obsm.keys():
        logger.warning(
            f"Linear reduction '{linear_reduction}' not found in adata.obsm. "
            "Running PCA to generate it."
        )
        sc.pp.pca(adata, n_comps=n_pcs)
        linear_reduction = "X_pca"

    if basis is None:
        if nonlinear_reduction is not None:
            basis = nonlinear_reduction
        else:
            basis = "basis"
            adata.obsm["X_basis"] = adata.obsm[linear_reduction][
                :, 0:2
            ]
    scv.pl.utils.check_basis(adata, basis)

    if "spliced" not in adata.layers.keys():
        raise ValueError("'spliced' data must be provided.")

    if "unspliced" not in adata.layers.keys():
        raise ValueError("'unspliced' data must be provided.")

    if type(mode) is str:
        mode = [mode]

    mode.append(fitting_by)
    if kinetics is True or denoise is True:
        mode.append("dynamical")

    mode = list(set(mode))
    if "dynamical" in mode:
        mode.sort(key="dynamical".__eq__)

    adata.obs[group_by] = adata.obs[group_by].astype(dtype="category")
    scv.pl.proportions(adata, groupby=group_by, save=False, show=False)

    plt.savefig(
        ".".join(filter(None, [fileprefix, "proportions.png"])), dpi=dpi
    )

    logger.info("- Filtering and normalizing data ...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts)

    logger.info("- Running moments ...")
    # adata.var['highly_variable_genes'].astype(bool)
    # adata.var['highly_variable_genes'].fillna(False, inplace=True)
    scv.pp.moments(
        adata, n_pcs=n_pcs, n_neighbors=n_neighbors, use_rep=linear_reduction
    )

    highly_variable_genes = adata.var["highly_variable_genes"].index.tolist()
    adata.uns["layer_features_RNA"] = highly_variable_genes
    adata.uns["layer_features_spliced"] = highly_variable_genes
    adata.uns["layer_features_unspliced"] = highly_variable_genes

    for m in mode:
        vkey_list = [m]
        dk_list = [False]
        gene_subset_list = [None]
        autoscale_list = [True]

        logger.info(f"- mode: {m}")
        adata.uns["layer_features_" + m] = highly_variable_genes
        adata.uns["layer_features_variance_" + m] = highly_variable_genes
        if m == "dynamical":
            adata2 = adata[:, adata.var[fitting_by + "_genes"]].copy()
            Ms = adata2.layers["Ms"]
            Mu = adata2.layers["Mu"]
            adata2.layers.clear()
            adata2.layers["Ms"] = Ms
            adata2.layers["Mu"] = Mu
            connectivities = adata2.obsp["connectivities"]
            adata2.obsp.clear()
            adata2.obsp["connectivities"] = connectivities
            adata.uns["layer_features_Ms"] = highly_variable_genes
            adata.uns["layer_features_Mu"] = highly_variable_genes

            scv.tl.recover_dynamics(
                adata2,
                var_names=fitting_by + "_genes",
                use_raw=False,
                n_jobs=ncores,
            )

            var_add = [
                i
                for i in list(adata2.var.columns)
                if not i in list(adata.var.columns)
            ]
            adata.var = adata.var.merge(
                adata2.var[var_add], how="left", left_index=True, right_index=True
            )
            adata.uns["recover_dynamics"] = adata2.uns["recover_dynamics"]

            adata.varm["loss"] = np.empty(
                (adata.shape[1], adata2.varm["loss"].shape[1])
            )
            adata.varm["loss"][:] = np.nan
            adata.varm["loss"][adata.var[fitting_by + "_genes"], :] = adata2.varm[
                "loss"
            ]

            empty_layer = np.empty((adata.layers["spliced"].shape))
            empty_layer[:] = np.nan
            adata.layers["fit_t"] = adata.layers["fit_tau"] = adata.layers[
                "fit_tau_"
            ] = empty_layer
            adata.layers["fit_t"][:, adata.var[fitting_by + "_genes"]] = (
                adata2.layers["fit_t"]
            )
            adata.layers["fit_tau"][:, adata.var[fitting_by + "_genes"]] = (
                adata2.layers["fit_tau"]
            )
            adata.layers["fit_tau_"][:, adata.var[fitting_by + "_genes"]] = (
                adata2.layers["fit_tau_"]
            )
            adata.uns["layer_features_fit_t"] = highly_variable_genes
            adata.uns["layer_features_fit_tau"] = highly_variable_genes
            adata.uns["layer_features_fit_tau_"] = highly_variable_genes

            if kinetics is True:
                vkey_list.append("dynamical_kinetics")
                dk_list.append(True)
                gene_subset_list.append(None)
                autoscale_list.append(True)
                top_genes = (
                    adata.var["fit_likelihood"]
                    .sort_values(ascending=False)
                    .index[:kinetics_topn]
                )
                scv.tl.differential_kinetic_test(
                    adata, var_names=top_genes, groupby=group_by
                )

            if denoise is True:
                vkey_list.append("dynamical_denoise")
                dk_list.append(False)
                gene_subset_list.append(
                    adata.var["fit_likelihood"]
                    .sort_values(ascending=False)
                    .index[:denoise_topn]
                )
                autoscale_list.append(False)
                adata.layers["dynamical_denoise"] = adata.layers[m] + np.random.normal(
                    adata.layers[m], scale=adata.layers["Ms"].std(0)
                )
                adata.uns["layer_features_dynamical_denoise"] = highly_variable_genes

        for i in range(len(vkey_list)):
            vkey = vkey_list[i]
            dk = dk_list[i]
            gene_subset = gene_subset_list[i]
            autoscale = autoscale_list[i]

            # Velocity graph
            scv.tl.velocity(adata, mode=m, vkey=vkey, diff_kinetics=dk)
            scv.tl.velocity_graph(
                adata,
                vkey=vkey,
                gene_subset=gene_subset,
                n_neighbors=n_neighbors,
                n_jobs=ncores,
            )
            if m == "dynamical":
                adata.var["velocity_genes"] = adata.var[m + "_genes"]
                adata.layers["velocity"] = adata.layers[m]
                adata.layers["variance_u"] = adata.layers[m + "_u"]
                adata.uns["layer_features_velocity"] = highly_variable_genes
                adata.uns["layer_features_variance_u"] = highly_variable_genes
                adata.uns["layer_features_dynamical_u"] = highly_variable_genes
            else:
                adata.var["velocity_gamma"] = adata.var[m + "_gamma"]
                adata.var["velocity_r2"] = adata.var[m + "_r2"]
                adata.var["velocity_genes"] = adata.var[m + "_genes"]
                adata.layers["velocity"] = adata.layers[m]
                # adata.layers["variance_velocity"] = adata.layers["variance_" + m]
                adata.uns["layer_features_velocity"] = highly_variable_genes

            # Velocity embedding
            scv.tl.velocity_embedding(
                adata, basis=basis, vkey=vkey, autoscale=autoscale
            )
            scv.pl.velocity_embedding_stream(
                adata,
                vkey=vkey,
                basis=basis,
                title=vkey,
                color=group_by,
                palette=palette,
                smooth=stream_smooth,
                density=stream_density,
                legend_loc="none",
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(filter(None, [fileprefix, vkey + "_stream.png"])),
                dpi=dpi,
            )

            scv.pl.velocity_embedding(
                adata,
                vkey=vkey,
                basis=basis,
                title=vkey,
                color=group_by,
                palette=palette,
                arrow_length=arrow_length,
                arrow_size=arrow_size,
                density=arrow_density,
                linewidth=0.3,
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(filter(None, [fileprefix, vkey + "_arrow.png"])),
                dpi=dpi,
            )

            scv.pl.velocity_embedding_grid(
                adata,
                vkey=vkey,
                basis=basis,
                title=vkey,
                color=group_by,
                palette=palette,
                arrow_length=arrow_length / 2,
                arrow_size=arrow_size / 2,
                density=arrow_density * 2,
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(
                    filter(None, [fileprefix, vkey + "_embedding_grid.png"])
                ),
                dpi=dpi,
            )

            # Velocity confidence
            scv.tl.velocity_confidence(adata, vkey=vkey)
            scv.pl.scatter(
                adata,
                basis=basis,
                title=vkey + " length",
                color=vkey + "_length",
                cmap="coolwarm",
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(filter(None, [fileprefix, vkey + "_length.png"])),
                dpi=dpi,
            )

            scv.pl.scatter(
                adata,
                basis=basis,
                title=vkey + " confidence",
                color=vkey + "_confidence",
                cmap="magma",
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(filter(None, [fileprefix, vkey + "_confidence.png"])),
                dpi=dpi,
            )

            # Terminal states
            for term in [
                "root_cells",
                "end_points",
                vkey + "_root_cells",
                vkey + "_end_points",
            ]:
                if term in adata.obs.columns:
                    adata.obs.drop(term, axis=1, inplace=True)

            scv.tl.terminal_states(
                adata,
                vkey=vkey,
            )
            for term in ["root_cells", "end_points"]:
                adata.obs[vkey + "_" + term] = adata.obs[term]
                adata.obs.drop(term, axis=1, inplace=True)

            # scv.pl.scatter(adata,basis=basis,title=vkey+" terminal_states",color_gradients=[vkey+'_root_cells', vkey+'_end_points'], legend_loc="best", save=False, show=False)
            # if show_plot is True:
            #   plt.show()
            # if save:
            #   plt.savefig('.'.join(filter(None, [fileprefix, vkey+"_terminal_states.png"])), dpi=dpi)

            # Pseudotime
            scv.tl.velocity_pseudotime(
                adata,
                vkey=vkey,
                root_key=vkey + "_root_cells",
                end_key=vkey + "_end_points",
            )
            scv.pl.scatter(
                adata,
                basis=basis,
                title=vkey + " pseudotime",
                color=vkey + "_pseudotime",
                cmap="cividis",
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(filter(None, [fileprefix, vkey + "_pseudotime.png"])),
                dpi=dpi,
            )

            # Latent time
            if m == "dynamical":
                scv.tl.latent_time(
                    adata,
                    vkey=vkey,
                    root_key=vkey + "_root_cells",
                    end_key=vkey + "_end_points",
                )
                scv.pl.scatter(
                    adata,
                    basis=basis,
                    title=vkey + " latent time",
                    color="latent_time",
                    color_map="cividis",
                    save=False,
                    show=False,
                )
                plt.savefig(
                    ".".join(
                        filter(None, [fileprefix, vkey + "_latent_time.png"])
                    ),
                    dpi=dpi,
                )

            # PAGA
            adata.uns["neighbors"]["distances"] = adata.obsp["distances"]
            adata.uns["neighbors"]["connectivities"] = adata.obsp["connectivities"]
            scv.tl.paga(
                adata,
                groups=group_by,
                vkey=vkey,
                root_key=vkey + "_root_cells",
                end_key=vkey + "_end_points",
            )
            scv.pl.paga(
                adata,
                title=vkey + " PAGA (" + group_by + ")",
                node_colors=palette,
                basis=basis,
                alpha=0.5,
                min_edge_width=2,
                node_size_scale=1.5,  # type: ignore
                legend_loc="none",
                save=False,
                show=False,
            )
            plt.savefig(
                ".".join(filter(None, [fileprefix, vkey + "_paga.png"])),
                dpi=dpi,
            )

            # Velocity genes
            if calculate_velocity_genes is True:
                if m != "dynamical":
                    scv.tl.rank_velocity_genes(adata, vkey=vkey, groupby=group_by)
                    adata.var[vkey + "_score"] = adata.var["spearmans_score"]
                    df1 = scv.get_df(adata.uns["rank_velocity_genes"]["names"])
                    adata.uns["rank_" + vkey + "_genenames"] = df1
                    df2 = scv.get_df(adata.uns["rank_velocity_genes"]["scores"])
                    adata.uns["rank_" + vkey + "_genescores"] = df2
                    del adata.uns["rank_velocity_genes"]
                else:
                    scv.tl.rank_dynamical_genes(adata, groupby=group_by)
                    df1 = scv.get_df(adata.uns["rank_dynamical_genes"]["names"])
                    adata.uns["rank_" + vkey + "_genenames"] = df1
                    df2 = scv.get_df(adata.uns["rank_dynamical_genes"]["scores"])
                    adata.uns["rank_" + vkey + "_genescores"] = df2
                    del adata.uns["rank_dynamical_genes"]

                for cluster in df1.columns:
                    # df1[0:1].values.ravel()[:12] ### by row

                    scv.pl.scatter(
                        adata,
                        color=group_by,
                        palette=palette,
                        basis=df1[cluster].values[:top_n],
                        vkey=vkey,
                        size=10,
                        linewidth=2,
                        alpha=1,
                        ylabel="cluster: " + cluster + "\nunspliced",
                        add_linfit=True,
                        add_rug=True,
                        add_outline=True,
                        ncols=3,
                        frameon=True,
                        save=False,
                        show=False,
                    )
                    plt.savefig(
                        ".".join(
                            filter(
                                None,
                                [fileprefix, cluster, vkey + "_genes1.png"],
                            )
                        ),
                        dpi=dpi,
                    )

                    scv.pl.velocity(
                        adata,
                        color=group_by,
                        var_names=df1[cluster].values[:top_n],
                        vkey=vkey,
                        size=10,
                        linewidth=2,
                        alpha=1,
                        ylabel="cluster: " + cluster + "\nunspliced",
                        add_outline=True,
                        basis=basis,
                        color_map=["Blues", "YlOrRd"],
                        ncols=2,
                        save=False,
                        show=False,
                    )
                    plt.savefig(
                        ".".join(
                            filter(
                                None,
                                [fileprefix, cluster, vkey + "_genes2.png"],
                            )
                        ),
                        dpi=dpi,
                    )

    try:
        adata.__dict__["_raw"].__dict__["_var"] = (
            adata.__dict__["_raw"]
            .__dict__["_var"]
            .rename(columns={"_index": "features"})
        )
    except:
        pass

    return adata


sobjfile: str = {{in.sobjfile | quote}}  # pyright: ignore  # noqa: E999
outfile: str = {{out.outfile | quote}}  # pyright: ignore  # noqa: E999
outdir: str = os.path.dirname(outfile)

ncores: int = {{envs.ncores | repr}}  # pyright: ignore  # noqa: E999
group_by: str | None = {{envs.group_by | repr}}  # pyright: ignore  # noqa: E999
mode: str | list[str] = {{envs.mode | repr}}  # pyright: ignore  # noqa: E999
fitting_by: str = {{envs.fitting_by | repr}}  # pyright: ignore  # noqa: E999
min_shared_counts: int = {{envs.min_shared_counts | repr}}  # pyright: ignore  # noqa: E999
n_pcs: int = {{envs.n_pcs | repr}}  # pyright: ignore  # noqa: E999
n_neighbors: int = {{envs.n_neighbors | repr}}  # pyright: ignore  # noqa: E999
denoise: bool = {{envs.denoise | repr}}  # pyright: ignore  # noqa: E999
denoise_topn: int = {{envs.denoise_topn | repr}}  # pyright: ignore  # noqa: E999
kinetics: bool = {{envs.kinetics | repr}}  # pyright: ignore  # noqa: E999
kinetics_topn: int = {{envs.kinetics_topn | repr}}  # pyright: ignore  # noqa: E999
calculate_velocity_genes: bool = {{envs.calculate_velocity_genes | repr}}  # pyright: ignore  # noqa: E999
top_n: int = {{envs.top_n | repr}}  # pyright: ignore  # noqa: E999
rscript: str = {{envs.rscript | repr}}  # pyright: ignore  # noqa: E999


if sobjfile.endswith(".h5ad"):
    h5ad_file = Path(sobjfile)
else:
    h5ad_file = Path(outfile).with_suffix(".input.h5ad")
    logger.info("Converting Seurat object to AnnData (h5ad) format...")
    seurat_ident_col = convert_seurat_to_anndata(
        input_file=sobjfile,
        output_file=h5ad_file,
        rscript=rscript,
        return_ident_col=not group_by,
    )
    group_by = group_by or seurat_ident_col

if group_by is None:
    group_by = "seurat_clusters"
    logger.warning(
        "`envs.group_by` is not provided. "
        "Using 'seurat_clusters' as the default groupby column. "
        "It is recommended to provide the `envs.group_by` parameter."
    )

logger.info(f"Reading AnnData (h5ad) file ...")
adata = sc.read_h5ad(h5ad_file)

if group_by not in adata.obs.columns:
    raise ValueError(
        f"The group_by column envs.group_by = '{group_by}' is not found in the AnnData object."
    )

logger.info(f"Running scVelo analysis ...")

if isinstance(mode, str):
    mode = [mode]

if not all([m in ["deterministic","stochastic","dynamical"] for m in mode]):
    raise ValueError(
        "The 'envs.mode' parameter must be one or more of 'deterministic', 'stochastic', or 'dynamical'."
    )

if not fitting_by in ["deterministic","stochastic"]:
    raise ValueError(
        "The 'envs.fitting_by' parameter must be either 'deterministic' or 'stochastic'."
    )

adata = SCVELO(
    adata=adata,
    group_by=group_by,
    dirpath=outdir,
    linear_reduction="X_pca",
    mode=mode,
    fitting_by=fitting_by,
    min_shared_counts=min_shared_counts,
    n_pcs=n_pcs,
    n_neighbors=n_neighbors,
    stream_smooth=None,
    stream_density=2,
    arrow_size=5,
    arrow_length=5,
    arrow_density=0.5,
    denoise=denoise,
    denoise_topn=denoise_topn,
    kinetics=kinetics,
    kinetics_topn=kinetics_topn,
    calculate_velocity_genes=calculate_velocity_genes,
    top_n=top_n,
    ncores=ncores,
    logger=logger,
)

if outfile.endswith(".h5ad"):
    h5ad_file = Path(outfile)
else:
    h5ad_file = Path(outfile).with_suffix(".output.h5ad")

logger.info(f"Writing object to AnnData (h5ad) file ...")
adata.write_h5ad(h5ad_file)

if not outfile.endswith(".h5ad"):
    logger.info(f"Converting AnnData (h5ad) file to Seurat format ...")
    convert_anndata_to_seurat(
        input_file=h5ad_file,
        output_file=outfile,
        rscript=rscript,
    )
