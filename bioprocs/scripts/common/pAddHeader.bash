head -n {{args.n}} {{in.infile1 | quote}} > {{out.outfile | quote}}
cat {{in.infile2 | quote}} >> {{out.outfile | quote}}