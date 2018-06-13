from pyppl import Proc

"""
Downstream analysis after mutations called
"""

"""
@name:
	pMutSig
@description:
	MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.

	For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
	
	See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)
@input:
	`maffile:file`: mutation table
	`cvgfile:file`: coverage table
	`cvrfile:file`: covariates table
	`mutdict:file`: mutation_type_dictionary_file
	`chrdir:file`:  chr_files_hg18 or chr_files_hg19 
@output:
	`outdir:dir`: The output directory
@args:
	`mutsig`: The path to `run_MutSigCV.sh`, default: 'mutsig'
	`mcr`: The Matlab MCR path
@requires:
	[MutSing](http://archive.broadinstitute.org/cancer/cga/mutsig_download)
"""
pMutSig = Proc()
pMutSig.input     = "maffile:file, cvgfile:file, cvrfile:file, mutdict:file, chrdir:file"
pMutSig.output    = "outdir:dir:mutsig.{{#}}"
pMutSig.args      = {"mutsig": "mutsig", "mcr": ""}
pMutSig.script    = """
{{args.mutsig}} "{{args.mcr}}" "{{maffile}}" "{{cvgfile}}" "{{cvrfile}}" "{{outdir}}/{{maffile | fn}}" "{{mutdict}}" "{{chrdir}}"
"""

"""
@name:
	pVcf2Maf
@description:
	Convert a snpEff-annotated somatic mutation vcf file (with normal and tumor samples) to [maf](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) file
@input:
	`infile:file`: vcf file
@output:
	`outfile:file`: The maf file
@args:
    `vepdata`: The path of vep data. Default: "" (default data dir of vep)
    `vep`: The path of vep excutable. Default: "vep"
    `vcf2maf`: The path of vcf2maf excutable. Default: "vcf2maf.pl"
    `reffile`: The reference fasta file.
    `nthread`: The number of threads used by vep. Default: 1
    `filtervcf`: The filter vcf
	`params`: Other parameters for `vcf2maf.pl`, default: ""
@requires:
	[vcf2maf.py](https://github.com/mskcc/vcf2maf)
"""
pVcf2Maf = Proc()
pVcf2Maf.input     = "infile:file"
pVcf2Maf.output    = "outfile:file:{{infile | fn}}.maf"
pVcf2Maf.args      = {"vepdata":"", "vep": "vep", "vcf2maf": "vcf2maf.pl", "reffile": "", "filtervcf": "", "nthread": 1, "params": ""}
pVcf2Maf.script    = """
if [[ -z "{{args.reffile}}" ]]; then
    echo "Reference file is required." 1>&2
    exit 1
fi
hline=$(grep -m1 "^#CHROM\\sPOS" "{{infile}}" | tr "\\t" " ")
IFS=" "
read -ra terms <<< $hline
gnormal=${terms[${#terms[@]}-1]}
gtumor=${terms[${#terms[@]}-2]}
if [[ "$gnormal" == "NORMAL" ]]; then
    read -ra parts <<< $(basename "{{infile}}" | tr "." " ")
    sample=${parts[0]}
    if [[ "$sample" == *"-"* ]]; then
        IFS="-"
        read -ra tn <<< "$sample"
        tumor=${tn[0]}
        normal=${tn[1]}
        if [[ $tumor == $normal ]]; then
            tumor="${tumor}_TUMOR"
            normal="${normal}_NORMAL"
        fi
    else
        tumor="${sample}_TUMOR"
        normal="${sample}_NORMAL"
    fi
else
    tumor=$gtumor
    normal=$gnormal
fi

veppath=$(dirname $(which "{{args.vep}}"))
{{args.vcf2maf}} --input-vcf "{{infile}}" --output-maf "{{outfile}}" --tumor-id $tumor --normal-id $normal --vcf-tumor-id $gtumor --vcf-normal-id $gnormal --vep-path $veppath --vep-data "{{args.vepdata}}" --ref-fasta "{{args.reffile}}" --filter-vcf "{{args.filtervcf}}" --vep-forks {{args.nthread}}
"""

"""
@name:
	pMergeMafs
@description:
	Merge MAF files
@input:
	`indir:file`: The directory containing MAF files to be merged
@output:
	`outfile:file`: The merged MAF file
"""
pMergeMafs = Proc ()
pMergeMafs.input  = "indir:file"
pMergeMafs.output = "outfile:file:{{indir | fn}}.{{#}}.maf"
pMergeMafs.script = """
mafs=({{indir}}/*.maf)
head -2 ${mafs[0]} > "{{outfile}}"
for maf in {{indir}}/*.maf; do
	tail -n +3 $maf >> "{{outfile}}"
done
"""

"""
@name:
	pMutsig4Plot
@description:
	Prepare somatic mutations for  plotting
@input:
	`msdir:file`:   The mutsig output directory
@output:
	`outfile:file`:  The file for plotting
	```
	#PANEL: Somatic mutations
	#INFO: MT|PI
	#DESC: Mutation type|Putative impact
	# could also be bordercolor, you can have up to 4 shape features
	#TYPE: shape|bgcolor
	# could also be continuous
	# expressions for set: a,b,c
	#                 norminal: no
	#                 continuous: [0,1]
	#DATA: set|norminal
	#NCOL: 2|2
	#NAME_MT: Frameshift|Missense|Nonsense|Silent|Splice_site|TSS|Nonstop
	#NAME_PI: HIGH|MODERATE|LOW|MODIFIER
	#VALUE_MT: 0|1|20|13|4|17|14
	#EXP_MT: frameshift_variant,inframe_deletion,inframe_insertion|missense_variant,initiator_codon_variant,stop_retained_variant,rare_amino_acid_variant|stop_gained|synonymous_variant|splice_acceptor_variant,splice_donor_variant|start_lost,start_retained|stop_lost
	#
	Sample1	Sample2	Sample3	Sample4	Sample5
	ABC	missense_variant|HIGH	missense_variant|HIGH	...
	...
	```
@args:
	`topn`:     the cutoff to select genes. If it is >= 1, top N genes will be selected, otherwise, it will be used as pvalue cutoff. Default: .05
@requires:
	[`pyvcf`](https://github.com/jamescasbon/PyVCF)
"""
pMutsig4Plot = Proc ()
pMutsig4Plot.input   = "msdir:file"
pMutsig4Plot.output  = "outfile:file:{{msdir | fn}}.4mutplot"
pMutsig4Plot.lang    = "python"
pMutsig4Plot.args    = {"topn": 5E-2}
pMutsig4Plot.script  = """
from glob import glob
fout = open("{{outfile}}", "w")
fout.write ("#PANEL: Somatic mutations\\n")
fout.write ("#INFO: MC|MT|PI\\n")
fout.write ("#DESC: Mutation class|Mutation type|Putative impact\\n")
fout.write ("#TYPE: shape|bordercolor|bgcolor\\n")
fout.write ("#VTYPE: numeric|character|character\\n")
fout.write ("#NCOL: 1|2|2\\n")
fout.write ("#DATA: norminal|norminal|norminal\\n")
fout.write ("#NAME_MC: 3'Flank|3'UTR|5'Flank|5'UTR|Frame_Shift_Del|Frame_Shift_Ins|IGR|In_Frame_Del|In_Frame_Ins|Intron|Missense_Mutation|Nonsense_Mutation|Nonstop_Mutation|RNA|Silent|Splice_Region|Splice_Site|Translation_Start_Site\\n")
fout.write ("#NAME_MT: INS|DEL|SNP\\n")
fout.write ("#NAME_PI: HIGH|MODERATE|LOW|MODIFIER\\n")
fout.write ("#VALUE_MC: 0|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17\\n")
fout.write ("#VALUE_MT: red|yellow|green\\n")
fout.write ("#VALUE_PI: red|yellow|gray|green\\n")

genes      = []
i          = 0
with open (glob("{{msdir}}/*.sig_genes.txt")[0]) as f:
	for line in f:
		line  = line.strip()
		if not line or line.startswith("gene"): continue
		parts = line.split("\\t")
		i    += 1
		if {{args.topn}} < 1 and float(parts[13]) >= {{args.topn}}: continue
		if {{args.topn}} >= 1 and i > {{args.topn}}: continue
		gene  = parts[0]
		genes.append(gene)

import sys, os
f = open (glob("{{msdir}}/*.mutations.txt")[0])
data = {}
impacts = {'': -1, 'MODIFIER': 0, 'LOW': 1, 'MODERATE': 2, 'HIGH': 3}
for line in f:
	line   = line.strip()
	if not line: continue
	parts  = line.split("\\t")
	gene   = parts[0]
	if gene not in genes: continue
	mt     = parts[9]
	mc     = parts[8]
	pi     = parts[92]
	sample = parts[15] + '-' + parts[16]
	if not data.has_key(gene): data[gene] = {}
	if not data[gene].has_key(sample): data[gene][sample] = ['', '', '']
	if impacts[pi] < impacts[data[gene][sample][2]]: continue # not significant enough
	data[gene][sample] = [mc, mt, pi]
	
samples    = []
samples    = [x.keys() for x in data.values()]
samples    = reduce (lambda x,y: x+y, samples)
samples    = sorted(list(set(samples)))

fout.write ("Gene\\t" + "\\t".join(samples) +"\\n")

for gene in genes:
	if not data.has_key(gene): continue
	fout.write (gene)
	for sample in samples:
		if not data[gene].has_key(sample):
			fout.write ("\\t|||")
		else:
			fout.write ("\\t%s" % '|'.join(data[gene][sample]))
	fout.write ("\\n")
fout.close()
"""

"""
@name:
	pMutPlot
@description:
	Plot mutations
	```
	|           |             |           |           |---
	|- ftWidth -|  s   s   s  |- pnWidth -|- lgWidth -| snHeight
	|           |             |           |           |---
	    feature1
		feature2
	```
@input:
	`indir:file`:    The input directory containing plot files
@output:
	`outfile:file`:  The plot png file
"""
pMutPlot = Proc ()
pMutPlot.input  = "indir:file"
pMutPlot.output = "outfile:file:mutplot.{{#}}.png"
pMutPlot.args   = {"snHeight": 40, "ftWidth": 40, "pnWidth": 40, "lgWidth": 160, "cex": 1}
pMutPlot.lang   = "Rscript"
pMutPlot.script = """
data   = list()
npanel = 0
nrows  = 0
gap    = 2
cwidth = 16
smcex  = 0.8
samples = c()
lgheight = {{args.snHeight}}
setwd ("{{indir}}")
for (pfile in list.files()) {

	con = file(pfile, open = "r")
	while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
		oneLine = sub("\\\\s+$", "", oneLine)
		if (oneLine == "#") break
		
		if (substr(oneLine, 1, 7) == "#PANEL:") {
			panel = sub("^\\\\s+", "", substring(oneLine, 8))
			data[[panel]] = list()
			next
		}
		
		if (substr(oneLine, 1, 6) == "#INFO:") {
			info  = sub("^\\\\s+", "", substring(oneLine, 7))
			data[[panel]][["info"]] = unlist(strsplit(info, "\\\\|"))
			next
		}
		
		if (substr(oneLine, 1, 6) == "#DESC:") {
			desc = sub("^\\\\s+", "", substring(oneLine, 7))
			data[[panel]][["desc"]] = unlist(strsplit(desc, "\\\\|"))
			lgheight = lgheight + cwidth + 11*gap + 1
			next
		}
		
		if (substr(oneLine, 1, 6) == "#TYPE:") {
			type = sub("^\\\\s+", "", substring(oneLine, 7))
			data[[panel]][["type"]] = unlist(strsplit(type, "\\\\|"))
			next
		}
		
		if (substr(oneLine, 1, 6) == "#NCOL:") {
			type = sub("^\\\\s+", "", substring(oneLine, 7))
			data[[panel]][["ncol"]] = as.numeric(unlist(strsplit(type, "\\\\|")))
			next
		}
		
		if (substr(oneLine, 1, 7) == "#VTYPE:") {
			vtype = sub("^\\\\s+", "", substring(oneLine, 8))
			data[[panel]][["vtype"]] = unlist(strsplit(vtype, "\\\\|"))
			next
		}
		
		if (substr(oneLine, 1, 6) == "#DATA:") {
			datastr = sub("^\\\\s+", "", substring(oneLine, 7))
			data[[panel]][["data"]] = unlist(strsplit(datastr, "\\\\|"))
			next
		}
		
		if (substr(oneLine, 1, 6) == "#NAME_") {
			ps = unlist(regmatches(oneLine, regexpr(":", oneLine), invert = TRUE))
			t  = substring(ps[1], 7)
			ns = unlist(strsplit(sub("^\\\\s+", "", ps[2]), "\\\\|"))
			idx   = which (data[[panel]][["info"]] == t)
			if (!idx %in% seq(data[[panel]]$names)) data[[panel]]$names[[idx]] = list()
			data[[panel]]$names[[idx]] = ns
			lgheight = lgheight + 2*cwidth + 2*gap + ceiling(length(ns) / data[[panel]][["ncol"]][idx]) * (cwidth + 2*gap)
			next
		}
		
		if (substr(oneLine, 1, 7) == "#VALUE_") {
			ps = unlist(regmatches(oneLine, regexpr(":", oneLine), invert = TRUE))
			t  = substring(ps[1], 8)
			ns = unlist(strsplit(sub("^\\\\s+", "", ps[2]), "\\\\|"))
			idx   = which (data[[panel]][["info"]] == t)
			if (!idx %in% seq(data[[panel]]$vals)) data[[panel]]$vals[[idx]] = list()
			data[[panel]]$vals[[idx]] = ns
			if (data[[panel]]$vtype[idx] == "numeric") data[[panel]]$vals[[idx]] = as.numeric(ns)
			next
		}
		
		if (substr(oneLine, 1, 5) == "#EXP_") {
			ps = unlist(regmatches(oneLine, regexpr(":", oneLine), invert = TRUE))
			t  = substring(ps[1], 6)
			ns = unlist(strsplit(sub("^\\\\s+", "", ps[2]), "\\\\|"))
			idx   = which (data[[panel]][["info"]] == t)
			if (!idx %in% seq(data[[panel]]$exps)) data[[panel]]$exps[[idx]] = list()
			if (data[[panel]][["data"]][idx] == "set") {
				exps = list()
				for (i in 1:length(ns)) {
					exps[[i]] = unlist(strsplit(ns[i], ","))
				}
			} else if (data[[panel]][["data"]][idx] == "norminal") {
				exps = ""
			} else {
				exps = list()
				for (x in 1:length(ns)) {
					exp = ns[x]
					e3  = ifelse(substr(exp, 1, 1) == '[', 1, 0)
					e4  = ifelse(substr(exp, nchar(exp)-1, nchar(exp)-1) == ']', 1, 0)
					e12 = as.numeric(unlist(strsplit(substr(exp, 2, nchar(exp)-1), ",")))
					e1  = ifelse(is.na(e12[1]), -Inf, e12[1])
					e2  = ifelse(is.na(e12[2]), Inf, e12[2])
					exps[[x]] = c(e1, e2, e3, e4)
				}
			}
			
			data[[panel]]$exps[[idx]] = exps
			next
		}
	} 
	
	close(con)
	
	content = read.table (pfile, header=T, row.names=1, check.names=F, stringsAsFactors=F)
	
	data[[panel]][["content"]] = content
	nrows   = nrows + nrow(content)
	npanel  = npanel + 1
	samples = union (samples, colnames(content))
}
if (!"ncol" %in% names(data[[panel]])) data[[panel]][["ncol"]] = c (1,1)

ncols = length(samples)
width  = {{args.ftWidth}}  + ncols * (cwidth + gap) + gap + {{args.lgWidth}} + {{args.pnWidth}}
height = {{args.snHeight}} + nrows * (cwidth + gap) + npanel * cwidth / 2
height = max (height, lgheight)

m = max(width, height)
res       = 300
basewidth = 2000
rwidth    = 400  # pixels in R
basewidth = basewidth / rwidth * width
if (width == m) {
	pngwidth  = basewidth
	pngheight = basewidth/width * height
} else {
	pngheight = basewidth
	pngwidth  = basewidth/height * width
}

png (file="{{outfile}}", res=res, width=pngwidth, height=pngheight)
# ploting
par(mar=c(0,0,0,0))
plot (c(0, width), c(0, height), type="n", xlab="", ylab="", axes=F)
# plot sample names
snx = {{args.ftWidth}} + gap + (1:ncols) * (cwidth + gap) - gap - cwidth / 2
sny = rep (height - {{args.snHeight}}, ncols)
text (snx, sny, samples, cex={{args.cex}}, adj = c(0, 0.5), srt=60)

lastheight = height - {{args.snHeight}}
for (panel in names(data)) {
	content = data[[panel]][["content"]] # matrix
	info    = data[[panel]][["info"]]    # vector
	desc    = data[[panel]][["desc"]]    # vector
	type    = data[[panel]][["type"]]    # vector
	dtype   = data[[panel]][["data"]]    # vector
	lncol   = data[[panel]][["ncol"]]    # vector # cols of legends
	names   = data[[panel]][["names"]]   # list
	vals    = data[[panel]][["vals"]]    # list
	exps    = data[[panel]][["exps"]]    # list
	len     = length(info)
	
	# Fill missing data with NA
	msamps  = samples[!samples %in% colnames(content)]
	fillna  = matrix(NA, ncol=length(msamps), nrow=nrow(content))
	colnames(fillna) = msamps
	content = cbind (content, fillna)
	content = content[, match(colnames(content), samples)]
	crows   = row(content)
	ccols   = col(content)
	cnrow   = nrow(content)

	# draw hline
	lines (c({{args.ftWidth}} + gap, {{args.ftWidth}} + ncols * (cwidth + gap)), rep(lastheight - cwidth/4, 2))
	lastheight = lastheight - cwidth / 2
	
	# plot feature names
	rnames  = rownames (content)
	fnx     = rep({{args.ftWidth}}, cnrow) - 2*gap
	fny     = lastheight - crows[, 1] * (cwidth + gap) + gap + cwidth/2
	text (fnx, fny, rnames, cex={{args.cex}}, adj = c(1, 0.5))
	
	# plot cells
	rx1     = {{args.ftWidth}} + gap + (ccols - 1) * (cwidth + gap)
	rx1     = as.vector(rx1)
	ry1     = lastheight - (crows - 1) * (cwidth + gap)
	ry1     = as.vector(ry1)
	       
	rx2     = rx1 + cwidth
	ry2     = ry1 - cwidth
	
	tmpvs   = as.matrix(as.data.frame(strsplit(as.matrix(content), "\\\\|")))
	colnames (tmpvs) = NULL

	values  = list()
	for (i in 1:nrow(tmpvs)) {
		tvs = tmpvs[i, ]
		if (dtype[i] == "continuous") tvs = as.numeric (tmpvs[i, ])
		
		if (!i %in% seq(exps) || is.null(exps[[i]])) {
			idx = match (tvs, names[[i]])
		} else {
			vs  = matrix (tvs, ncol=ncols)
			exp = exps[[i]]
			idx = matrix(ncol=ncols, nrow=cnrow)
			if (dtype[i] == 'set') {
				for (j in seq(exp)) {
					idxtmp = matrix(vs %in% exp[[j]], ncol=ncols)
					idx[idxtmp] = j
				}
				values[[i]] = matrix(vals[[i]][idx], ncol=ncols)
			} else if (dtype[i] == 'norminal') {
				for (j in seq(exp)) {
					idxtmp = which(vs == exp[[j]])
					idx[idxtmp] = j
				}
			} else {
				for (j in seq(exp)) {
					e = exp[[j]]
					
					if (is.infinite(e[1])) {
						idxtmp = which (vs < e[2])
						if (e[4] == 1) idxtmp = which (vs <= e[2])
					} else if (is.infinite(e[2])) {
						idxtmp =  which (vs > e[1])
						if (e[3] == 1) idxtmp = which (vs >= e[1])
					} else {
						if (e[3] == 1) {
							idxtmp = which (vs >= e[1] & vs < e[2])
							if (e[4] == 1) idxtmp = which (vs >= e[1] & vs <= e[2])
						} else {
							idxtmp = which (vs > e[1] & vs < e[2])
							if (e[4] == 1) idxtmp = which (vs > e[1] & vs <= e[2])
						}
					}
					idx[idxtmp] = j
				}
			}
		}
		values[[i]] = vals[[i]][idx]
	}
	
	bgcolors = c()
	bdcolors = c()
	bgidx = which(type == "bgcolor")
	if (length(bgidx) > 0) 
		bgcolors = values[[bgidx[1]]] # only one allowed

	bdidx = which (type == 'bordercolor')
	if (length(bdidx) > 0) 
		bdcolors = values[[bdidx[1]]] # only one allowed
	
	rect (rx1, ry1, rx2, ry2, col = bgcolors, border=bdcolors)
	
	# get number of shape features
	nsidx  = which (type == 'shape')
	nshape = length(nsidx)
	if (nshape > 2) {
		dx = rx1 + cwidth / 4 + 1
		dy = ry1 - cwidth / 4 - 1
		points (dx, dy, pch = values[[nsidx[1]]], cex = smcex)
		
		dx = rx1 + 3*cwidth / 4 - 1
		dy = ry1 - cwidth / 4 - 1
		points (dx, dy, pch = values[[nsidx[2]]], cex = smcex)
		
		dx = rx1 + cwidth / 4 + 1
		dy = ry1 - 3*cwidth / 4 + 1
		points (dx, dy, pch = values[[nsidx[3]]], cex = smcex)
		
		if (nshape > 3) {
			dx = rx1 + 3*cwidth / 4 - 1
			dy = ry1 - 3*cwidth / 4 + 1
			points (dx, dy, pch = values[[nsidx[4]]], cex = smcex)
		}
	} else if (nshape == 1) {
		dx = rx1 + cwidth / 2
		dy = ry1 - cwidth / 2
		points (dx, dy, pch = values[[nsidx[1]]])
	} else {
		dx = rx1 + cwidth / 4 + 1
		dy = ry1 - cwidth / 4 - 1
		points (dx, dy, pch = values[[nsidx[1]]], cex = smcex)
		
		dx = rx1 + 3*cwidth / 4 - 1
		dy = ry1 - 3*cwidth / 4 + 1
		points (dx, dy, pch = values[[nsidx[2]]], cex = smcex)
	}
	
	# plot panel name
	# draw vline
	lines (rep(max(rx2) + 2*gap, 2), c(lastheight, min(ry2)), lwd=2)
	text (max(rx2) + 4*gap, (lastheight + min(ry2))/2, panel, cex={{args.cex}}, adj = c(0.5, 0), srt=270)
	lastheight2 = lastheight
	
	# plot legends
	lgLeft  = max(rx2) + 2*gap  + {{args.pnWidth}}
	for (i in 1:len) {
		des = desc[i]
		t   = type[i]
		ns  = unlist(names[i])
		vs  = unlist(vals[i])
		nns = length(ns)
		nc  = lncol[i]
		
		text (lgLeft, lastheight, des, cex={{args.cex}}, adj = c(0, 0))
		lines (c(lgLeft, lgLeft + {{args.lgWidth}}), rep(lastheight - 3*gap, 2))
		lastheight = lastheight - 8*gap 
		
		nwidth = {{args.lgWidth}} / nc
		nsidx  = match (nsidx, nsidx)
		ncx    = lgLeft + ( (1:nns+1)%%nc ) * nwidth
		ncy    = lastheight - floor((1:nns - 1) / nc) * (cwidth + 2*gap)
		ncx2   = ncx + cwidth
		ncy2   = ncy - cwidth
		
		if (t == 'shape') {
			rect (ncx, ncy, ncx2, ncy2)
			# points
			pidx = length(which (type[1:i] == 'shape'))
			if (nshape > 2) {
				if (pidx == 1 || pidx == 3) {
					px = ncx + cwidth / 4 + 1
				} else {
					px = ncx + 3*cwidth / 4 - 1
				}
				if (pidx == 1 || pidx == 2) {
					py = ncy - cwidth / 4 - 1
				} else {
					py = ncy - 3*cwidth / 4 + 1
				}
				cex = smcex
			} else if (nshape == 1) {
				px = ncx + cwidth / 2
				py = ncy - cwidth / 2
				cex = 1
			} else {
				if (pidx == 1) {
					px = ncx + cwidth / 4 + 1
					py = ncy - cwidth / 4 - 1
				} else {
					px = ncx + 3*cwidth / 4 - 1
					py = ncy - 3*cwidth / 4 + 1
				}
				cex = smcex
			}
			points (px, py, pch = vs, cex = cex)
		} else if (t == 'bgcolor') {
			rect (ncx, ncy, ncx2, ncy2, col = vs)
		} else if (t == 'bordercolor') {
			rect (ncx, ncy, ncx2, ncy2, border = vs)
		}
		
		ntx   = ncx + cwidth + 2*gap
		text (ntx, ncy - cwidth/2, ns, cex = {{args.cex}}, adj=c(0, 0.5))		
		
		lastheight = min (ncy) - 2*cwidth - 2*gap
	}
	lastheight = min (lastheight, lastheight2 - cnrow * (cwidth + gap) - cwidth / 2)
}
dev.off()
"""

"""
@name:
	pCepip
@description:
	run CEPIP.
@input:
	`avinput:file`: The avinput file
	`cell`:         The cell
@output:
	`outfile:file`: The cepip result file
@args:
	`bin-cepip`:    The jar file path of cepip, default: /data2/junwenwang/shared/tools/cepip/cepip.jar
@requires:
	[`cepip`](http://jjwanglab.org/cepip/)
"""
pCepip = Proc()
pCepip.input  = "avinput:file, cell"
pCepip.output = "outfile:file:{{avinput | fn}}.cepip.flt.txt"
pCepip.args   = {"bin-cepip": "/data2/junwenwang/shared/tools/cepip/cepip.jar"}
pCepip.script = """
cd "{{args.bin-cepip | dirname}}"
java -jar "{{args.bin-cepip | bn}}" --annovar-file "{{avinput}}" --regulatory-causing-predict all --cell {{cell}} --db-score dbncfp --out "{{outfile | prefix | [:-4]}}"
"""

