from bioprocs.utils import alwaysList
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
suffix  = {{args.suffix | repr}}
if not suffix:
	suffix = ['']
suffix = alwaysList(suffix)

reader = TsvReader(infile, skip = 1, comment = 'CDE_ID:')
writer = TsvWriter(outfile)
writer.cnames = ['Sample', 'FID', 'IID', 'PID', 'MID', 'Sex', 'Pheno']
writer.writeHead()

for r in reader:
	for suf in suffix:
		rec        = TsvRecord()
		rec.Sample = r.bcr_patient_barcode + suf
		rec.FID    = rec.Sample
		rec.IID    = rec.Sample
		rec.PID    = 0
		rec.MID    = 0
		rec.Pheno  = -9
		rec.Sex    = 1 if r.gender == 'MALE' else 2 if r.gender == 'FEMALE' else 0
		writer.write(rec)
writer.close()

