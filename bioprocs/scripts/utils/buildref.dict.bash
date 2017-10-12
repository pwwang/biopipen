ref="{{args.ref}}"
dict="${ref%.*}.dict"
if [[ ! -e "$dict" ]]; then
	echo "Generating dict index for $ref ..." 1>&2
	{{args.picard}} CreateSequenceDictionary R="$ref" O="$dict"
fi