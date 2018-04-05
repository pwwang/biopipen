function reffai_check() {
	if [[ -e $1 ]]; then
		return 0
	else
		return 1
	fi
}

# $1: The ref file
# $2: The index file
# $3: samtools
function reffai_do() {
	samtools="samtools"
	if [[ $# -gt 2 ]]; then
		samtools=$3
	fi
	$samtools faidx $1
	if [[ "$1.fai" != "$2" ]]; then
		mv "$1.fai" "$2"
	fi
}