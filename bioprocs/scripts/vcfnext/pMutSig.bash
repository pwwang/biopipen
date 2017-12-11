cmd='{{args.mutsig}} {{args.mcr | quote}} {{in.infile | quote}} {{args.cvrg | quote}} {{args.cvrt | quote}} "{{out.outdir}}/{{in.infile | fn}}" {{args.mutdict | quote}} {{args.chrdir | quote}}'
echo $cmd
exec $cmd