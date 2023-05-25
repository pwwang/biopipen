<script>
    import { Image } from "$libs";
</script>

{% for gwpng in job.out.outdir | glob: "*/*.png" %}
{%  set sample = gwpng | dirname | basename %}
<h1>{{sample}}</h1>
<Image src="{{gwpng}}" />
{% endfor %}
