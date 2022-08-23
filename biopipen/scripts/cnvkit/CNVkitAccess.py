import cmdy

excfiles = {{in.excfiles | repr}}  # pyright: ignore
outfile = {{out.outfile | quote}}  # pyright: ignore
reffile = {{envs.ref | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
min_gap_size = {{envs.min_gap_size | quote}}  # pyright: ignore


def main():
    other_args = {}
    if excfiles:
        other_args["exclude"] = excfiles

    cmdy.cnvkit.access(
        **other_args,
        s=min_gap_size,
        o=outfile,
        _=reffile,
        _exe=cnvkit,
    ).fg()


if __name__ == "__main__":
    main()
