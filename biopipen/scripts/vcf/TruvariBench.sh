# shellcheck disable=SC1083
compvcf={{in.compvcf | quote}}
basevcf={{in.basevcf | quote}}
outdir={{out.outdir | quote}}
truvari={{envs.truvari | quote}}
ref={{envs.ref | quote}}
refdist={{envs.refdist | quote}}
pctseq={{envs.pctseq | quote}}
pctsize={{envs.pctsize | quote}}
pctovl={{envs.pctovl | quote}}
sizemax={{envs.sizemax | default: 50000 | quote}}
# shellcheck disable=SC1054
{% if envs.typeignore %}
typeignore="--typeignore"
{% else %}
typeignore=""
{% endif %}
{% if envs.multimatch %}
multimatch="--multimatch"
# shellcheck disable=SC1009
{% else %}
multimatch=""
# shellcheck disable=SC1073
{% endif %}

rm -rf $outdir
cmd="$truvari bench \
    -c '$compvcf' \
    -b '$basevcf' \
    -f '$ref' \
    --refdist $refdist \
    --pctseq $pctseq \
    --pctsize $pctsize \
    --pctovl $pctovl \
    --sizemax $sizemax \
    $typeignore \
    $multimatch \
    -o $outdir"

echo "$cmd"
eval "$cmd"
