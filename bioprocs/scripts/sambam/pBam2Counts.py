from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

params = {{args.params}}

####### htseq
{% if args.tool | lambda x: x == 'htseq' %}

{% if in.infile.endswith('.bam') %}
params['f'] = 'bam'
{% endif %}

cmd = '{{args.htseq}} %s "{{in.infile}}" "{{args.refgene}}" > "{{out.outfile}}"' % (cmdargs(params))
runcmd (cmd)

{% endif %}
