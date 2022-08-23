import cmdy

bamfile = {{in.bamfile | quote}}  # pyright: ignore
target_file = {{in.target_file | quote}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
reffile = {{envs.reffile | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
count = {{envs.count | repr}}  # pyright: ignore
min_mapq = {{envs.min_mapq | repr}}  # pyright: ignore
ncores = {{envs.ncores | repr}}  # pyright: ignore


def main():

    cmdy.cnvkit.coverage(
        f=reffile,
        c=count,
        q=min_mapq,
        p=ncores,
        o=outfile,
        _=[bamfile, target_file],
        _exe=cnvkit,
    ).fg()


if __name__ == "__main__":
    main()
