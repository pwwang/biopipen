ref="{{args.ref}}"
fai="$ref.fai"
if [[ ! -e "$fai" ]]; then
	echo "Generating fai index for $ref ..." 1>&2
	{{args.samtools}} faidx "$ref"
fi