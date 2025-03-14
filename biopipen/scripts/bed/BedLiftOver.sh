# shellcheck disable=SC2148,SC1083
inbed={{ in.inbed | quote }}
outbed={{ out.outbed | quote }}
rejfile={{ job.outdir | joinpaths: "rejected.bed" | quote }}
liftover={{ envs.liftover | quote }}
chain={{ envs.chain | quote }}

# shellcheck disable=SC2154
$liftover \
    "$inbed" \
    "$chain" \
    "$outbed" \
    "$rejfile"
