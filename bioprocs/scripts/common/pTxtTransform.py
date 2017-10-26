{{txtTransform}}

#txtTransform(infile, outfile, transformer = None, header = True, delimit = "\t")
txtTransform({{in.txtfile | quote}}, {{out.outfile | quote}}, transformer = {{args.transformer}}, header = {{args.header}}, delimit = {{args.delimit | quote}}, comment = {{args.comment | quote}})