from pyppl import proc
"""
A set of TFBS procs
"""

"""
@name:
	pMotifScanByMEME
@input:
	`mfile:file`: The motif file
	`sfile:file`: The sequence file
@output:
	`outdir:file`: The output dir
@args:
	`bin-fimo`: The path of `fimo` executable, default: "fimo"
	`params`:   Other parameters for `fimo`
@requires:
	[`fimo` from MEME Suite](http://meme-suite.org/tools/fimo)
"""
pMotifScanByMEME = proc ()
pMotifScanByMEME.input  = "mfile:file, sfile:file"
pMotifScanByMEME.output = "outdir:file:{{mfile | fn}}-{{sfile | fn}}.fimo"
pMotifScanByMEME.args   = {"params": "", "bin-fimo": "fimo"}
pMotifScanByMEME.script = """
if [[ -e "{{outdir}}" ]]; then rm -rf "{{outdir}}"; fi
{{proc.args.bin-fimo}} --o "{{outdir}}" {{proc.args.params}} "{{mfile}}" "{{sfile}}"
"""

pMS2Bed = proc ()
pMS2Bed.input  = "msdir:file"
pMS2Bed.output = "outfile:file:{{msdir | fn}}.bed"
pMS2Bed.lang   = "python"
pMS2Bed.script = """
#pattern name   sequence name   start   stop    strand  score   p-value q-value matched sequence
#TFDP1   TERT::chr5:1293184-1297184  2266    2287    -   22.1515 5.06e-09    1.95e-05    GCGAGCGGCGCGCGGGCGGGGA
#TFDP1   TERT::chr5:1293184-1297184  1986    2007    +   20.6212 2.97e-08    5.71e-05    CGCCGCGAGGAGAgggcggggc
#TFDP1   TERT::chr5:1293184-1297184  1706    1727    +   20.0758 5.26e-08    6.74e-05    CGGAAGGAgggggcggcggggg
#TBX15   TERT::chr5:1293184-1297184  3888    3906    -   19.422  1.24e-07    0.000529    ACGGGGGTGGGGGTGGGGT
# sometimes "TERT::" part is missing

with open ("{{msdir}}/fimo.txt") as f, open("{{outfile}}", "w") as fout:
	for line in f:
		line   = line.strip()
		if not line or line.startswith("#"): continue
		parts  = line.split()
		out    = [''] * 9
		range  = parts[1].split("::")[1] if "::" in parts[1] else parts[1]
		suffix = "::" + parts[1].split("::")[0] if "::" in parts[1] else ""
		ranges = range.split(':')
		out[0] = ranges[0]
		start  = ranges[1].split('-')[0]
		out[1] = str(int(start) + int(parts[2]))
		out[2] = str(int(start) + int(parts[3]))
		out[3] = parts[0] + suffix
		out[4] = parts[5]
		out[5] = parts[4]
		out[6] = parts[6]
		out[7] = parts[7]
		out[8] = parts[8]
		fout.write ("\\t".join(out) + "\\n")
"""

"""
@name:
	pMEMEMDB2Gene
@input:
	`memefile:file`: The meme motif downloaded from MEME website
	`species`:       The species, can be tax id, or common name
	- only can be one of fruitfly, mouse, human, frog, zebrafish, thale-cress, pig, rat, nematode
	- check them from http://http://mygene.info/v3/metadata
@output:
	`outfile:file`:  The output file containing the name pairs
@requires:
	[python-mygene](https://pypi.python.org/pypi/mygene)
"""
pMEMEmDB2Gene = proc ()
pMEMEmDB2Gene.input  = "memefile:file, species"
pMEMEmDB2Gene.output = "outfile:file:{{ memefile | readlink | dirname | bn }}-{{memefile | fn}}.m2gene"
pMEMEmDB2Gene.lang   = "python"
pMEMEmDB2Gene.script = """
memedir = "{{ memefile | readlink | dirname | bn }}"
from mygene import MyGeneInfo
mg = MyGeneInfo()
import re

def queryGene (q):
	g = mg.query(q, scopes="symbol,alias", fields="symbol", species="{{species}}", size=1)
	if not g.has_key("hits"): return q
	if not g["hits"]: return q
	return g['hits'][0]['symbol']

def getGene (t, md, extra=""):
	if md == "ARABD": return queryGene (t[2])
	if md == "CIS-BP" or md == "CISBP-RNA":
		if extra.startswith("Drosophila_"):            return queryGene (t[2])
		if extra.startswith("Mus_musculus"):           return queryGene (t[2])
		if extra.startswith("Homo_sapiens"):           return queryGene (t[2]) 
		if extra.startswith("Xenopus_"):               return queryGene (t[2]) 
		if extra.startswith("Danio_rerio"):            return queryGene (t[2]) 
		if extra.startswith("Arabidopsis_"):           return queryGene (t[2]) 
		if extra.startswith("Sus_scrofa"):             return queryGene (t[2]) 
		if extra.startswith("Rattus_norvegicus"):      return queryGene (t[2]) 
		if extra.startswith("Caenorhabditis_elegans"): return queryGene (t[2])
		return t[2]
	if md == "ECOLI": return t[1]
	if md == "EUKARYOTE": 
		if extra.startswith("hallikas2006"):           return t[1]
		if extra.startswith("homeodomain"):            return queryGene (t[2].split("_")[0])
		if extra.startswith("jolma2010"):              return queryGene (t[2].split("_")[0])
		if extra.startswith("jolma2013"):              return queryGene (t[1].split("_")[0])
		if extra.startswith("macisaac_theme"):	       return queryGene (t[1].split("_")[0])
		if extra.startswith("prodoric"):               return t[2]
		if extra.startswith("regtransbase"):           return t[1].split("_")[0]
		if extra.startswith("SwissRegulon_human_and_mouse"):
			# ADNP_IRX_SIX_ZHX / DMAP1_NCOR{1,2}_SMARC / CEBPA,B_DDIT3 / ELK1,4_GABP{A,B1} / GLI1..3 / POU5F1_SOX2{dimer} / TCF4_dimer
			ret = []
			q   = t[1][:-3].replace("{dimer}", "").replace("{mouse}", "").replace("_dimer", "")
			qs  = q.split("_")
			for q in qs:
				if "{" in q:
					base, rest = q.split('{')
					rest = rest[:-1]
					for r in rest.split(","):
						q   = base + r
						ret.append (queryGene(q))
				elif ".." in q and "," in q:
					base,med,end = re.split(r'(?:\\.\\.|,)', q)
					start = int (base[-1])
					base  = base[:-1]
					med   = int (med)
					end   = int (end)
					rs    = range(start, med+1) + [end]
					for r in rs:
						q = base + str(r)
						ret.append (queryGene(q))
				elif ".." in q:
					base, end = q.split("..")
					start = int (base[-1])
					base  = base[:-1]
					end   = int (end)
					for r in range (start, end+1):
						q = base + str(r)
						ret.append (queryGene(q))
				elif "," in q:
					parts = q.split(",")
					start = parts[0][-1]
					base  = parts[0][:-1]
					parts[0] = start
					for r in parts:
						q = base + r
						ret.append (queryGene(q))
				else:
					ret.append (queryGene(q))
			return filter(None, ret)
		if extra.startswith("wei2010_human_mws"):      return queryGene (t[1].split("-")[1])   
		if extra.startswith("wei2010_mouse_mws"):      return queryGene (t[1].split("-")[1])   
		if extra.startswith("wei2010_mouse_pbm"):      return queryGene (t[1].split(",")[0])
		if extra.startswith("zhao2011"):               return queryGene (t[1].split("_")[1])
	if md == 'FLY':
		if extra.startswith("dmmpmm2009"):             return queryGene (t[1])
		if extra.startswith("flyreg"):                 return queryGene (t[1])
		if extra.startswith("idmmpmm2009"):            return queryGene (t[1])
		if extra.startswith("fly_factor_survey"):      return queryGene (t[2].split("_")[0])
		if extra.startswith("OnTheFly_2014_Drosophi"): return queryGene (t[2].split("_")[0])
	if md == 'HUMAN':
		if extra.startswith("HOCOMOCOv9"):             return queryGene (t[2])
		if extra.startswith("HOCOMOCOv10_HUMAN_mono"): return queryGene (t[1].split("_")[0])
	if md == 'JASPAR':
		q   = t[2] if len(t) > 2 else t[1]
		q   = q.split('(')[0]
		ret = []
		qs  = q.split('::')
		for q in qs:
			ret.append (queryGene(q))
		return filter(None, ret)
	if md == 'MALARIA':             return queryGene (t[2])
	if md == 'MIRBASE':             return t[2]
	if md == 'MOUSE':
		if extra.startswith("chen2008"):               return queryGene (t[1])
		if extra.startswith("HOCOMOCOv10_MOUSE_mono"): return queryGene (t[1].split("_")[0])
		if extra.startswith("uniprobe_mouse"):         return queryGene (t[2].split("_")[0])
	if md == 'PROKARYOTE':
		if extra.startswith("collectf"):              return queryGene (t[2].split("_")[0])
		if extra.startswith("prodoric"):              return queryGene (t[2])
		if extra.startswith("regtransbase"):          return queryGene (t[1].split("_")[0])
	if md == 'RNA':                return t[2]
	if md == 'TFBSshape':       
		if extra.startswith("TFBSshape_JASPAR"):      return queryGene (t[2])
		if extra.startswith("TFBSshape_UniPROBE"):    return queryGene (t[2].split("_")[0])
	if md == 'WORM':
		if extra.startswith("TFBSshape_JASPAR"):      return queryGene (t[2])
	if md == 'YEAST':
		if extra.startswith("macisaac_yeast"):        return queryGene (t[1])
		if extra.startswith("scpd_matrix"):           return queryGene (t[1])
		if extra.startswith("SwissRegulon_s_cer"):    return queryGene (t[1])
		if extra.startswith("yeast_uniprobe_GR09"):   return queryGene (t[2])
		if extra.startswith("YEASTRACT_20130918"):    return queryGene (t[1].split("&")[0])
	return t[1]
	
with open ("{{memefile}}", "r") as f, open ("{{outfile}}", "w") as fout:
	for line in f:
		if not line.startswith ("MOTIF"): continue
		t    = line.strip().split(" ")
		gene = getGene (t, memedir, "{{memefile | fn}}")
		if not gene: continue
		if not isinstance (gene, list): gene = [gene]
		fout.write ("%s\\t%s\\n" % (t[1], "|".join(gene)))
"""