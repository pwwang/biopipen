{{args.bedtools}} random -l {{in.l}} -n {{in.n}} -seed {{args.seed}} -g {{args.gsize | quote}} > {{out.outfile | quote}}
