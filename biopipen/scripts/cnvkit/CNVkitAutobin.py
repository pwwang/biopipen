import cmdy

bamfiles = {{in.bamfiles | repr}}  # pyright: ignore
accfile = {{in.accfile | quote}}  # pyright: ignore
baitfile = {{in.baitfile | repr}}  # pyright: ignore
target_file = {{out.target_file | quote}}  # pyright: ignore
antitarget_file = {{out.antitarget_file | quote}}  # pyright: ignore
reffile = {{envs.ref | quote}}  # pyright: ignore
cnvkit = {{envs.cnvkit | quote}}  # pyright: ignore
method = {{envs.method | quote}}  # pyright: ignore
bp_per_bin = {{envs.bp_per_bin | repr}}  # pyright: ignore
target_max_size = {{envs.target_max_size | repr}}  # pyright: ignore
target_min_size = {{envs.target_min_size | repr}}  # pyright: ignore
antitarget_max_size = {{envs.antitarget_max_size | repr}}  # pyright: ignore
antitarget_min_size = {{envs.antitarget_min_size | repr}}  # pyright: ignore
annotate = {{envs.annotate | repr}}  # pyright: ignore
short_names = {{envs.short_names | repr}}  # pyright: ignore


def main():

    cmdy.cnvkit.autobin(
        f=reffile,
        m=method,
        g=accfile,
        t=baitfile or False,
        b=bp_per_bin,
        target_max_size=target_max_size,
        target_min_size=target_min_size,
        antitarget_max_size=antitarget_max_size,
        antitarget_min_size=antitarget_min_size,
        annotate=annotate or False,
        short_names=short_names,
        target_output_bed=target_file,
        antitarget_output_bed=antitarget_file,
        _=bamfiles,
        _exe=cnvkit,
    ).fg()


if __name__ == "__main__":
    main()
