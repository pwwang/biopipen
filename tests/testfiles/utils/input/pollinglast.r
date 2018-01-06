{{pollingLast}}

cmd = 'sleep $((1-{{job.index}})); echo {{in.in}} > {{out.outfile}};'

pollingLast ({{proc.workdir | quote}}, {{proc.size}}, {{job.index}}, cmd, "flag.done")

{% if job.index, proc.size | lambda x, y: x == y - 1 %}
files = Sys.glob("{{proc.workdir}}/*/output/outfile")
ret = NULL
for (f in files) {
	ret = c(ret, readLines(f))
}
print (ret)
writeLines(ret, files[1])
{% endif %}