if [[ ! -e "{{args.ref}}" ]]; then
	echo "Reference file is not specfied or not exists." 1>&2
	exit 1
fi
