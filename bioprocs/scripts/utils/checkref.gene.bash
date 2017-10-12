if [[ ! -e "{{args.refgene}}" ]]; then
	echo "Reference gene file is not specfied or not exists." 1>&2
	exit 1
fi