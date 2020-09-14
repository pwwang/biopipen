from pathlib import Path
from diot import Diot
from bioprocs.utils import shell2 as shell, mem2

infiles = {{i.infiles | repr}}
outfile = {{o.outfile | quote}}
gatk3   = {{args.gatk | quote}}
ref     = {{args.ref | quote}}
tmpdir  = {{args.tmpdir | quote}}
mem     = {{args.mem | quote}}
params  = {{args.params | repr}}

shell.load_config(gatk3 = gatk3)

# "[" and "]" in filename will cause gatk error:
# MESSAGE: Invalid command line: Failed to parse value org.broadinstitute.gatk.utils.commandline.ArgumentMatchStringValue@d190639 for argument variantCollections. Message: Illegal character in path
# so symbolically link and rename them
for i, infile in enumerate(infiles):
	infile = Path(infile)
	if '[' in infile.name or ']' in infile.name:
		newfile = infile.parent.joinpath(infile.name.replace('[', '_').replace(']', '_'))
		newfile.symlink_to(infile)
		infiles[i] = newfile

params.T = 'CombineVariants'
params.R = ref
params.o = outfile
params.variant = infiles

shell.gatk3(f'-Djava.io.tmpdir={tmpdir}',
			*mem2(mem, 'java').split(), **params).fg
