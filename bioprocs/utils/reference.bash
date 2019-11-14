
function logger() {
	flag="INFO"
	if [[ $# -gt 1 ]]; then
		flag="$2"
	fi
	printf "%s: %s\n" $flag "$1" 1>&2
}

function ref_check() {
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

function reference_bwa() {
	if [[ -z "$bwa" ]]; then bwa="bwa"; fi
	ref=$1
	ref2check="$1.pac"
	if [[ ! -e "$ref2check" ]]; then
		logger "BWA index not exists, creating ..." "WARNING"
		$bwa index -a bwtsw -b 10000000000 "$ref"
	fi
}

function reference_bowtie2() {
	if [[ -z "$bowtie2" ]]; then bowtie2="bowtie2"; fi
	if [[ -z "$nthread" ]]; then nthread=1; fi
	bowtie2_build="${bowtie2}-build"
	ref=$1
	if [[ ! -e "$ref.bt2" && ! -e "$ref.1.bt2" ]]; then
		logger "Bowtie2 index not exists, creating ..." "WARNING"
		$bowtie2_build --threads $nthread $ref $ref
	fi
}

function reference_ngm() {
	if [[ -z "$ngm" ]]; then ngm="ngm"; fi
	if [[ -z "$nthread" ]]; then nthread=1; fi
	ref=$1
	if [[ ! -e "$ref-enc.2.ngm" ]]; then
		logger "NGM index not exists, creating ..." "WARNING"
		$ngm -r $ref -t $nthread
	fi
}

function reference_star() {
	if [[ -z "$star" ]]; then star="star"; fi
	if [[ -z "$nthread" ]]; then nthread=1; fi
	ref=$1
	refdir="${ref%.fa*}.star"
	if [[ ! -d "$refdir" ]]; then
		logger "Star index not exists, creating ..." "WARNING"
		mkdir -p "$refdir"
		$star --runMode genomeGenerate --genomeDir "$refdir" --genomeFastaFiles $ref --sjdbGTFfile $refexon --runThreadN $nthread
	fi
}

function reference_fasta() {
	if [[ -z "$samtools" ]]; then samtools="samtools"; fi
	ref=$1
	if [[ ! -e "$ref.fai" ]]; then
		logger "Fasta index not exists, creating ..." "WARNING"
		$samtools faidx $ref "$ref.fai"
	fi
}

function reference_picard() {
	if [[ -z "$picard" ]]; then picard="picard"; fi
	ref=$1
	if [[ ! -e "${ref%.fa*}.dict" ]]; then
		logger "Picard index not exists, creating ..." "WARNING"
		$picard CreateSequenceDictionary REFERENCE=$ref OUTPUT="${ref%.fa*}.dict"
	fi
}

function reference_elprep() {
	if [[ -z "$elprep" ]]; then elprep="elprep"; fi
	ref=$1
	if [[ ! -e "$ref.elprep" ]]; then
		logger "Elprep index not exists, creating ..." "WARNING"
		$elprep fasta-to-elfasta $ref "$ref.elprep"
	fi
}

function reference_kallisto() {
	ref=$1
	if [[ ! -e "$ref" ]]; then
		logger "Kallisto index does not exist, please create it before run this process."
		logger "Check: kallisto index ..."
		logger "  Please use a reference with gene sequences only!"
		logger "  Whole genome sequence reference will take forever!"
		exit 1
	fi
}

function reference() {
	tool="reference_$1"
	shift;
	$tool $@
}
