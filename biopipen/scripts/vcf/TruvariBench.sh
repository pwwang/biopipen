compvcf={{in.compvcf | quote}}
basevcf={{in.basevcf | quote}}
outdir={{out.outdir | quote}}
truvari={{envs.truvari | quote}}
ref={{envs.ref | quote}}
refdist={{envs.refdist | quote}}
pctsim={{envs.pctsim | quote}}
pctsize={{envs.pctsize | quote}}
pctovl={{envs.pctovl | quote}}
sizemax={{envs.sizemax | default: 50000 | quote}}
{% if envs.typeignore %}
typeignore="--typeignore"
{% else %}
typeignore=""
{% endif %}

rm -rf $outdir
$truvari bench \
    -c "$compvcf" \
    -b "$basevcf" \
    -f "$ref" \
    --refdist $refdist \
    --pctsim $pctsim \
    --pctsize $pctsize \
    --pctovl $pctovl \
    --sizemax $sizemax \
    $typeignore \
    -o $outdir
