invcf={{ in.invcf | quote }}
outvcf={{ out.outvcf | quote }}
rejfile={{ job.outdir | joinpaths: "rejected.vcf" | quote }}
gatk={{ envs.gatk | quote }}
chain={{ envs.chain | quote }}
reffa={{ envs.reffa | quote }}
args={{ envs.args | dict_to_cli_args | quote }}

refdict="${reffa%.fa}.dict"
if [[ ! -e "$refdict" ]]; then
    echo "Sequence dictionary does not exist: $refdict" 1>&2
    exit 1
fi

$gatk LiftoverVcf \
    $args \
    --INPUT "$invcf" \
    --OUTPUT "$outvcf" \
    --REJECT "$rejfile" \
    --REFERENCE_SEQUENCE "$reffa" \
    --CHAIN "$chain" \
