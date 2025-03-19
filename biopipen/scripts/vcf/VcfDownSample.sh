# shellcheck disable=SC2148
# shellcheck disable=SC2036
# shellcheck disable=SC2030
# shellcheck disable=SC1083
# shellcheck disable=SC2288
infile={{in.infile | quote}}
outfile={{out.outfile | quote}}
n={{envs.n}}

# shellcheck disable=SC2031
if [[ $infile == *.gz ]]; then
    outfile=$(echo "$outfile" | sed -r "s/\.gz$//")
    # shellcheck disable=SC2126
    nheader=$(zcat "$infile" | head -n 9999 | grep "^#" | wc -l | cut -d' ' -f1)
    if [[ ! $n -gt 1 ]]; then
        nrows=$(zcat "$infile" | wc -l | cut -d' ' -f1)
        # shellcheck disable=SC2004
        nvars=$(($nrows - $nheader))
        n=$(echo "$nvars * $n" | bc)
    fi
    zcat "$infile" | head -n "$nheader" > "$outfile"
    # shellcheck disable=SC2004
    zcat "$infile" | tail -n +$(($nheader + 1)) | shuf -n "$n" | LC_ALL=C sort -k1,1V -k2,2n >> "$outfile"
    bgzip "$outfile"
else
    # shellcheck disable=SC2126
    nheader=$(head -n 9999 "$infile" | grep "^#" | wc -l | cut -d' ' -f1)
    if [[ ! $n -gt 1 ]]; then
        nrows=$(wc -l "$infile" | cut -d' ' -f1)
        # shellcheck disable=SC2004
        nvars=$(($nrows - $nheader))
        n=$(echo "$nvars * $n" | bc)
    fi
    head -n "$nheader" "$infile" > "$outfile"
    # shellcheck disable=SC2004
    tail -n +$(($nheader + 1)) "$infile" | shuf -n "$n" | LC_ALL=C sort -k1,1V -k2,2n >> "$outfile"
fi
