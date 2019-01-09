function reffai_check() {
	if [[ -e "$1" ]]; then
		echo 0
	else
		echo 1
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
	if [[ ! -e "$1" ]]; then
		echo "Reference file not exists: $1"
		exit 1
	fi
	$samtools faidx "$1"
	if [[ "$1.fai" != "$2" ]]; then
		mv "$1.fai" "$2"
	fi
}

function refdict_check() {
	if [[ -e "$1" ]]; then
		echo 0
	else
		echo 1
	fi
}

# $1: The ref file
# $2: The index file
# $3: picard
function refdict_do() {
	picard="picard"
	if [[ $# -gt 2 ]]; then
		picard=$3
	fi
	if [[ ! -e "$1" ]]; then
		echo "Reference file not exists: $1"
		exit 1
	fi
	$picard CreateSequenceDictionary R="$1" O="$2"
}

function refkallisto_check() {
	if [[ -e "$1" ]]; then
		echo 0
	else
		echo 1
	fi
}

function refkallisto_do() {
	kallisto="kallisto"
	if [[ $# -gt 2 ]]; then
		kallisto=$3
	fi
	if [[ ! -e "$1" ]]; then
		echo "Reference file not exists: $1"
		exit 1
	fi
	$kallisto index -i "$1" "$2"
}
