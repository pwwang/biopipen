<script>
    import { Image, DataTable } from "@@";
</script>

<h1>Total number of arm-level gains/losses</h1>

<h2>Plots</h2>
<Image src="{{job.out.outdir}}/CAAs.png" />

{% if job.in.metafile %}
<Image src="{{job.out.outdir}}/CAAs_group.png" />
{% endif %}

<h2>Tables</h2>
<h3>arm</h3>
<DataTable src="{{ job.out.outdir | joinpaths: 'AS_arm.txt' }}"
    data={ {{ job.out.outdir | joinpaths: 'AS_arm.txt' | datatable: sep="\t" }} } />

<h3>seg</h3>
<DataTable src="{{ job.out.outdir | joinpaths: 'AS_seg.txt' }}"
    data={ {{ job.out.outdir | joinpaths: 'AS_seg.txt' | datatable: sep="\t" }} } />

<h2>Arm-level gain/loss breakdown</h2>

<h3>arm</h3>

<DataTable src="{{ job.out.outdir | joinpaths: 'CAA_arm.txt' }}"
    data={ {{ job.out.outdir | joinpaths: 'CAA_arm.txt' | datatable: sep="\t" }} } />

<h3>seg</h3>

<DataTable src="{{ job.out.outdir | joinpaths: 'CAA_seg.txt' }}"
    data={ {{ job.out.outdir | joinpaths: 'CAA_seg.txt' | datatable: sep="\t" }} } />
