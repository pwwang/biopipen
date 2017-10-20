string = {{in.in | quote}}

{% for b in args.breakOn %}
string = string.replace({{b | quote}}, '\n')
{% endfor %}

{% if args.trimLine %}
string = '\n'.join([x.strip() for x in string.splitlines()])
{% endif %}

with open('{{out.outfile}}', 'w') as fout:
	fout.write(string)