from glob import glob
from os import path, symlink
from sys import stderr

for fname in {{i.infiles}}:
	bn  = path.basename (fname)
	dst = path.join ({{o.outdir | quote}}, bn)
	if path.exists (dst):
		if path.samefile(fname, dst):
			continue

		fn, ext = path.splitext(bn)
		dsts    = glob (path.join("{{o.outdir}}", fn + "[[]*[]]" + ext))

		if not dsts:
			dst = path.join("{{o.outdir}}", fn + "[1]" + ext)
		else:
			if any([path.samefile(fname, d) for d in dsts]):
				continue

			maxidx = max([int(path.basename(d)[len(fn)+1 : -len(ext)-1]) for d in dsts])
			dst = path.join("{{o.outdir}}", fn + "[" + str(maxidx+1) + "]" + ext)
		stderr.write ("pyppl.log.warning: rename %s to %s\n" % (bn, path.basename(dst)))
	symlink (fname, dst)
