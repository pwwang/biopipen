{{ "plot.R" | rimport }}

ggs <- {{ args.ggs | R }}
infile <- {{ i.infile | R }}
outfile <- {{ o.outfile | R }}
inopts <- {{ args.inopts | R }}
params <- {{ args.params | R }}
intype <- {{ args.intype | R }}
devpars <- {{ args.devpars | R }}
data <- read.table.inopts(infile, inopts)

plot.pie(data, outfile,
         params = list(intype = intype),
         ggs = ggs,
         devpars = devpars)
