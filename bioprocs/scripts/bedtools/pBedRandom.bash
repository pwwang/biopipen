{{args.bedtools}} random -l {{i.l}} -n {{i.n}} -seed {{args.seed}} -g {{args.gsize | quote}} > {{o.outfile | quote}}
