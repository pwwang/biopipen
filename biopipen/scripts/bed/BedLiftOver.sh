inbed={{ in.inbed | quote }}
outbed={{ out.outbed | quote }}
rejfile={{ job.outdir | joinpaths: "rejected.bed" | quote }}
liftover={{ envs.liftover | quote }}
chain={{ envs.chain | quote }}

$liftover \
    $inbed \
    $chain \
    $outbed \
    $rejfile
