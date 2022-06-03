vcfs={{in.vcfs | join: " " | quote}}
outfile={{out.outfile | quote}}
truvari={{envs.truvari | quote}}

$truvari consistency \
    $vcfs > "$outfile"
