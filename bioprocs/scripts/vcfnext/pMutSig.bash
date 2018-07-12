mutsig={{args.mutsig | quote}}
mcr={{args.mcr | quote}}
infile={{in.infile | quote}}
cvrg={{args.cvrg | quote}}
cvrt={{args.cvrt | quote}}
out="{{out.outdir}}/{{in.infile | fn}}"
mutdict={{args.mutdict | quote}}
chrdir={{args.chrdir | quote}}
mafile="$infile"
if echo $infile | grep '.gz$' >/dev/null; then
	mafile="{{job.outdir}}/{{in.infile | fn}}" # no .gz
	gunzip -c "$infile" > "$mafile"
	infile="$mafile"
fi
cmd="$mutsig '$mcr' '$infile' '$cvrg' '$cvrt' '$out' '$mutdict' '$chrdir'"
echo $cmd
exec $cmd
