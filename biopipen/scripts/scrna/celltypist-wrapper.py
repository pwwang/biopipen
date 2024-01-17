from argparse import ArgumentParser

parser = ArgumentParser(description="Run CellTypist")
parser.add_argument(
    "-i", "--input", required=True, help="Input H5AD file with AnnData object"
)
parser.add_argument("-o", "--output", required=True, help="Output file")
parser.add_argument("-m", "--model", required=True, help="Model file")
parser.add_argument(
    "-v", "--majority_voting",
    action="store_true",
    help="Majority voting"
)
parser.add_argument(
    "-c", "--over_clustering",
    default="seurat_clusters",
    help="Over clustering. Ignored if the column does not exist."
)


if __name__ == "__main__":
    import scanpy as sc
    import celltypist

    args = parser.parse_args()
    adata = sc.read_h5ad(args.input)
    over_clustering = args.over_clustering
    if over_clustering and over_clustering not in adata.obs.columns:
        print("WARNING: Over clustering column not found. Ignoring over clustering.")
        over_clustering = None

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
            out_adata._raw._var.rename(columns={"_index": "features"}, inplace=True)
            del out_adata.raw
        except (KeyError, AttributeError):
            pass

        out_adata.write(args.output)
    else:
        out_adata.obs.to_csv(args.output, sep="\t", index=True)
