from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

params = {{args.params}}

####### htseq
{% if args.tool == 'htseq' %}

{% if i.infile.endswith('.bam') %}
params['f'] = 'bam'
{% endif %}

cmd = '{{args.htseq}} %s "{{i.infile}}" "{{args.refgene}}" > "{{o.outfile}}"' % (cmdargs(params))
runcmd (cmd)

{% endif %}
