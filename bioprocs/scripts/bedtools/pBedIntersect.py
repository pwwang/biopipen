from diot import Diot
from bioprocs.utils import shell2 as shell

params = {{args.params | repr}}

params['a']           = {{i.afile | quote}}
params['b']           = {{i.bfile | quote}}
params['wa']          = params.get('wa', True)
params['wb']          = params.get('wb', True)
params['nonamecheck'] = params.get('nonamecheck', True)
# params['_out']        = {{o.outfile | quote}}
params['_debug']      = True

shell.load_config(bedtools = {{args.bedtools | quote}})
shell.bedtools.intersect(**params) > {{o.outfile | quote}}

# attach header to the output file if input is a VCF file
# do it only when wa is True and wb is False and a is VCF file, or
# wb is True and wa is False and b is VCF file

def attach_header(ofile, vfile):
	"""Attach the header of vfile to ofile"""
	tmpfile = ofile + '.tmp'
	if vfile.endswith('.gz'):
		shell.zcat(vfile).p | shell.grep('^#') ^ tmpfile
	else:
		shell.grep('^#', vfile).r > tmpfile
	shell.cat(ofile, tmpfile) >> tmpfile
	shell.mv(tmpfile, ofile)

if params.wa and params.wb:
	pass
elif params.wb and (params.b.endswith('.vcf') or params.b.endswith('.vcf.gz')):
	attach_header(params._out, params.b)
elif (params.a.endswith('.vcf') or params.a.endswith('.vcf.gz')):
	attach_header(params._out, params.a)
