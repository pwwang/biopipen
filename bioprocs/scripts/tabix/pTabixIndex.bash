infile={{i.infile | quote}}
outfile={{o.outfile | quote}}
outidx={{o.outidx | quote}}

cat >"{{job.outdir}}/params.py" <<EOF
from sys import stdout
{{params2CmdArgs}}
stdout.write(params2CmdArgs({{args.params}}))
EOF

tabix={{args.tabix}} $({{args.python}} "{{job.outdir}}/params.py")
if [[ "$infile" == *".gz" ]]; then
	out=$($tabix "$infile" 2>&1)
	if [[ "$out" == "Not a BGZF file"* ]]; then
		tmpfile="{{job.outdir}}/{{i.infile | fn}}.tmp"
		gunzip "$infile" -c > "$tmpfile"
		bgzip "$tmpfile" -c > "$outfile"
	else
		ln -s "$outfile" "$infile" 
	fi
else
	bgzip "$infile" -c > "$outfile"
fi
exec $tabix "$outfile"