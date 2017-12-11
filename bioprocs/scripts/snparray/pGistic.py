from os import path
from subprocess import check_output

{{runcmd}}
{{params2CmdArgs}}

gisticbin = check_output(['which', {{args.gistic | quote}}]).strip()
gisticdir = path.dirname(path.realpath(gisticbin))
gisticbin = path.join(gisticdir, 'gp_gistic2_from_seg')

ldpath  = 'LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{{args.mcr}}/runtime/glnxa64:'
ldpath += '{{args.mcr}}/bin/glnxa64:'
ldpath += '{{args.mcr}}/sys/os/glnxa64'

gisticbin = ldpath + ' ' + gisticbin
refgene   = path.join(gisticdir, 'refgenefiles', {{args.genome | quote}} + ".mat")

params = {{args.params}}

params['b']       = {{out.outdir | quote}}
params['seg']     = {{in.segfile | quote}}
params['refgene'] = refgene
{% if in.mkfile %}
params['mk']      = {{in.mkfile | quote}}
{% endif %}
{% if in.alfile %}
params['alf']      = {{in.alfile | quote}}
{% endif %}
{% if in.cnvfile %}
params['cnv']      = {{in.cnvfile | quote}}
{% endif %}

cmd = "%s %s" % (gisticbin, params2CmdArgs(params, dash = '-', equal = ' '))
runcmd(cmd)
