# shellcheck disable=SC2148
# shellcheck disable=SC1083
invcf={{ in.invcf | quote }}
outvcf={{ out.outvcf | quote }}
rejfile={{ job.outdir | joinpaths: "rejected.vcf" | quote }}
gatk={{ envs.gatk | quote }}
chain={{ envs.chain | quote }}
reffa={{ envs.reffa | quote }}
args={{ envs.args | dict_to_cli_args: join=True }}

# shellcheck disable=SC2154
refdict="${reffa%.fa}.dict"
if [[ ! -e "$refdict" ]]; then
    echo "Sequence dictionary does not exist: $refdict" 1>&2
    exit 1
fi

# shellcheck disable=SC2154
# shellcheck disable=SC2086
$gatk LiftoverVcf \
    $args \
    --INPUT "$invcf" \
    --OUTPUT "$outvcf" \
    --REJECT "$rejfile" \
    --REFERENCE_SEQUENCE "$reffa" \
    --CHAIN "$chain" \
