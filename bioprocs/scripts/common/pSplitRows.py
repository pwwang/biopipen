from os import path

infile = {{i.infile | quote}}
outdir = {{o.outdir | quote}}
n      = {{args.n}}
header = {{args.cnames}}
nskip  = {{args.skip}}

# get # lines and header
nlines   = 0
headerow = ''
with open(infile) as f:
	for line in f: 
		if header and nlines == nskip: headerow = line
		nlines += 1
		
nlines -= nskip
if header: nlines -= 1

nrow, remain = divmod(nlines, n)
nrows = [0] * n
for i in range(n):
	nrows[i] = nrow + 1 if i < remain else nrow

with open(infile) as f:
	for _ in range(nskip): f.readline()
	if header: f.readline()
	for i in range(n):
		subfile = path.join(outdir, "{{i.infile | fn}}-%s{{i.infile | ext}}" % (i+1))
		with open(subfile, 'w') as fout:
			if header: fout.write(headerow)
			for _ in range(nrows[i]):
				fout.write(f.readline())

