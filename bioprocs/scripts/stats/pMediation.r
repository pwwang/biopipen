library(methods)
library(mediation)

set.seed(8525)

inopts.cnames = {{args.inopts | lambda x: x.get('cnames', True) | R}}
inopts.rnames = {{args.inopts | lambda x: x.get('rnames', True) | R}}

infile  = {{in.infile | quote}}
outfile = {{out.outfile | quote}}

indata  = read.table(infile, header = inopts.cnames, row.names = if(inopts.rnames) 1 else NULL, check.names = F, sep = "\t")

attach(indata)

modelm   = {{args.medopts.modelm | lambda x: 'r:' + x | R}}
modely   = {{args.medopts.modely | lambda x: 'r:' + x | R}}
mediator = {{args.medopts.mediator | quote}}
treat    = {{args.medopts.treat | quote}}
boot     = {{args.medopts.boot | R}}
sims     = {{args.medopts.sims | R}}

output   = matrix(NA, ncol = 2, nrow = 16)
colnames(output) = c("Treatment", "Control")
rownames(output) = c(
	"ACME", "ACME_CI1", "ACME_CI2", "ACME_PVAL", 
	"ADE" , "ADE_CI1" , "ADE_CI2" , "ADE_PVAL", 
	"TE"  , "TE_CI1"  , "TE_CI2"  , "TE_PVAL", 
	"PM"  , "PM_CI1"  , "PM_CI2"  , "PM_PVAL" 
)
tryCatch({
	results  = mediate(modelm, modely, treat = treat, mediator = mediator, boot = boot, sims = sims)
	output["ACME",      "Treatment"] = results$d1
	output["ACME_CI1",  "Treatment"] = results$d1.ci[1]
	output["ACME_CI2",  "Treatment"] = results$d1.ci[2]
	output["ACME_PVAL", "Treatment"] = results$d1.p
	output["ACME",      "Control"]   = results$d0
	output["ACME_CI1",  "Control"]   = results$d0.ci[1]
	output["ACME_CI2",  "Control"]   = results$d0.ci[2]
	output["ACME_PVAL", "Control"]   = results$d0.p
	output["ADE",       "Treatment"] = results$z1
	output["ADE_CI1",   "Treatment"] = results$z1.ci[1]
	output["ADE_CI2",   "Treatment"] = results$z1.ci[2]
	output["ADE_PVAL",  "Treatment"] = results$z1.p
	output["ADE",       "Control"]   = results$z0
	output["ADE_CI1",   "Control"]   = results$z0.ci[1]
	output["ADE_CI2",   "Control"]   = results$z0.ci[2]
	output["ADE_PVAL",  "Control"]   = results$z0.p
	output["TE",        "Treatment"] = results$tau.coef
	output["TE_CI1",    "Treatment"] = results$tau.ci[1]
	output["TE_CI2",    "Treatment"] = results$tau.ci[2]
	output["TE_PVAL",   "Treatment"] = results$tau.p
	output["TE",        "Control"]   = results$tau.coef
	output["TE_CI1",    "Control"]   = results$tau.ci[1]
	output["TE_CI2",    "Control"]   = results$tau.ci[2]
	output["TE_PVAL",   "Control"]   = results$tau.p
	output["PM",        "Treatment"] = results$n1
	output["PM_CI1",    "Treatment"] = results$n1.ci[1]
	output["PM_CI2",    "Treatment"] = results$n1.ci[2]
	output["PM_PVAL",   "Treatment"] = results$n1.p
	output["PM",        "Control"]   = results$n0
	output["PM_CI1",    "Control"]   = results$n0.ci[1]
	output["PM_CI2",    "Control"]   = results$n0.ci[2]
	output["PM_PVAL",   "Control"]   = results$n0.p
}, error = function(e){
	output["ACME",      "Treatment"] = 0
	output["ACME_CI1",  "Treatment"] = 0
	output["ACME_CI2",  "Treatment"] = 0
	output["ACME_PVAL", "Treatment"] = 1
	output["ACME",      "Control"]   = 0
	output["ACME_CI1",  "Control"]   = 0
	output["ACME_CI2",  "Control"]   = 0
	output["ACME_PVAL", "Control"]   = 1
	output["ADE",       "Treatment"] = 0
	output["ADE_CI1",   "Treatment"] = 0
	output["ADE_CI2",   "Treatment"] = 0
	output["ADE_PVAL",  "Treatment"] = 1
	output["ADE",       "Control"]   = 0
	output["ADE_CI1",   "Control"]   = 0
	output["ADE_CI2",   "Control"]   = 0
	output["ADE_PVAL",  "Control"]   = 1
	output["TE",        "Treatment"] = 0
	output["TE_CI1",    "Treatment"] = 0
	output["TE_CI2",    "Treatment"] = 0
	output["TE_PVAL",   "Treatment"] = 1
	output["TE",        "Control"]   = 0
	output["TE_CI1",    "Control"]   = 0
	output["TE_CI2",    "Control"]   = 0
	output["TE_PVAL",   "Control"]   = 1
	output["PM",        "Treatment"] = 0
	output["PM_CI1",    "Treatment"] = 0
	output["PM_CI2",    "Treatment"] = 0
	output["PM_PVAL",   "Treatment"] = 1
	output["PM",        "Control"]   = 0
	output["PM_CI1",    "Control"]   = 0
	output["PM_CI2",    "Control"]   = 0
	output["PM_PVAL",   "Control"]   = 1

})


write.table(round(output, 3), outfile, col.names = T, row.names = T, sep = "\t", quote = F)