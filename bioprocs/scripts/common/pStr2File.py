string = {{i.instr | quote}}

{% for b in args.breakOn %}
string = string.replace({{b | quote}}, '\n')
{% endfor %}

{% if args.trimLine %}
string = '\n'.join([x.strip() for x in string.splitlines()]) + '\n'
{% endif %}

with open('{{o.outfile}}', 'w') as fout:
	fout.write(string)