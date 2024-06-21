# shellcheck disable=all
export infile={{in.infile | quote}}
export outfile={{out.outfile | quote}}
is_outdir={{envs.outdir | int}}
cmd_given={{envs.cmd | bool | int}}
{% set _ = out.outfile | dirname | joinpath: "cmd.sh" | as_path | attr: 'write_text' | call: envs.cmd %}
cmd="{{proc.lang}} {{out.outfile | dirname | joinpath: 'cmd.sh'}}"
if [[ "$cmd_given" -eq 0 ]]; then
    echo "No command given." 1>&2
    exit 1
fi
if [[ $is_outdir -eq 1 ]]; then
    mkdir -p "$outfile"
fi
eval "$cmd"
