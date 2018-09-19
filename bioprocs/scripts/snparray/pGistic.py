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

params['b']       = {{o.outdir | quote}}
params['seg']     = {{i.segfile | quote}}
params['refgene'] = refgene
{% if i.mkfile %}
params['mk']      = {{i.mkfile | quote}}
{% endif %}
{% if i.alfile %}
params['alf']      = {{i.alfile | quote}}
{% endif %}
{% if i.cnvfile %}
params['cnv']      = {{i.cnvfile | quote}}
{% endif %}

cmd = "%s %s" % (gisticbin, params2CmdArgs(params, dash = '-', equal = ' '))
runcmd(cmd)
