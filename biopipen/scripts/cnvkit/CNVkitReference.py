import cmdy

covfiles = {{in.covfiles | repr}}  # pyright: ignore
target_file = {{in.target_file | repr}}  # pyright: ignore
antitarget_file = {{in.antitarget_file | repr}}  # pyright: ignore
sample_sex = {{in.sample_sex | repr}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
reffile = {{envs.ref | repr}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
cluster = {{envs.cluster | repr}}  # pyright: ignore
min_cluster_size = {{envs.min_cluster_size | repr}}  # pyright: ignore
male_reference = {{envs.male_reference | repr}} # pyright: ignore
no_gc = {{envs.no_gc | repr}}  # pyright: ignore
no_edge = {{envs.no_edge | repr}}  # pyright: ignore
no_rmask = {{envs.no_rmask | repr}}  # pyright: ignore


def main():

    cmdy.cnvkit.reference(
        f=reffile,
        o=outfile,
        c=cluster,
        min_cluster_size=min_cluster_size,
        x=sample_sex or False,
        y=male_reference,
        t=False if not target_file or covfiles else target_file,
        a=False if not antitarget_file or covfiles else antitarget_file,
        no_gc=False,
        no_edge=False,
        no_rmask=False,
        _=covfiles or [],
        _exe=cnvkit,
    ).fg()


if __name__ == "__main__":
    main()
