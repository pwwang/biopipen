{{runcmd}}
{{params2CmdArgs}}

params = {{args.params}}

####### htseq
{% if args.tool | lambda x: x == 'htseq' %}

{% if in.infile.endswith('.bam') %}
params['f'] = 'bam'
{% endif %}

cmd = '{{args.htseq}} %s "{{in.infile}}" "{{args.refgene}}" > "{{out.outfile}}"' % (params2CmdArgs(params))
runcmd (cmd)

{% endif %}