{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "$libs";
    import { Tile } from "$ccs";
</script>


{%- macro report_job(job, h) -%}

    {%- set method_fullnames = {
        "chao1": "Chao1 index",
        "hill": "Hill numbers",
        "div": "True diversity",
        "gini.simp": "Gini-Simpson index",
        "inv.simp": "Inverse Simpson index",
        "gini": "Gini coefficient",
        "raref": "Rarefaction",
    } -%}

    {%- set method_desc = {
        "chao1": """
            <p>Chao1 index is a nonparameteric asymptotic estimator of species richness (number of species in a population).</p>
            <p><code>chao1 = S_obs + N_1(N_1-1)/(2*(N_2+1))</code></p>
            <p>where <code>N_1</code> and <code>N_2</code> are count of singletons and doubletons respectively, and <code>S_obs</code> is the abserved groups.</p>
        """,
        "hill": """
            <p>Hill numbers are a mathematically unified family of diversity indices (differing only by an exponent q).</p>
            <p>Also known as the Hill diversity index, is a measure of species diversity that takes into account both the richness and the evenness of a community.</p>
            <p>The Hill number is calculated using the following equation:</p>
            <p><code>H' = N^(1/q)</code></p>
            <p>where N is the total number of individuals in the sample and q is a parameter that controls the weight given to species richness versus species evenness.</p>
            <p>For q = 0, the Hill number reduces to the number of species in the community (i.e., species richness), while for q = 1, it reduces to the Shannon diversity index. For values of q between 0 and 1, the Hill number gives more weight to species evenness, while for values of q greater than 1, it gives more weight to species richness.</p>
        """,
        "div": """
            <p>True diversity, or the effective number of types, refers to the number of equally abundant types needed for the average proportional abundance of the types to equal that observed in the dataset of interest where all types may not be equally abundant.</p>
            <p>Also known as entropy or Shannon index. <code>H' = -Σ(p_i * log(p_i))</code></p>
            <p>It is zero if there is only one species present in the community.</p>
            <p>On the other hand, the Shannon diversity index will be maximized when all species in the community are equally abundant.</p>
        """,
        "gini.simp": """
            <p>The Gini-Simpson index is the probability of interspecific encounter, i.e., probability that two entities represent different types.</p>
            <p><code>Gini-Simpson = 1 - Σ(p_i^2)</code></p>
            <p>where <code>p_i</code> is the proportion of individuals in the sample belonging to the i-th species.</p>
            <p>The Gini-Simpson index ranges from 0 to 1, with higher values indicating higher diversity. A value of 0 indicates that all individuals belong to the same species, while a value of 1 indicates that all species are equally abundant.</p>
        """,
        "inv.simp": """
            <p>Inverse Simpson index is the effective number of types that is obtained when the weighted arithmetic mean is used to quantify average proportional abundance of types in the dataset of interest.</p>
            <p><code>1/D = Σ(p_i^2)</code></p>
            <p>where p_i is the proportion of individuals in the sample belonging to the i-th species.</p>
            <p>The Inverse Simpson index ranges from 1 to infinity, with higher values indicating lower diversity. A value of 1 indicates that all species are equally abundant, while a value of infinity indicates that all individuals belong to the same species.</p>
        """,
        "gini": """
            <p>The Gini coefficient measures the inequality among values of a frequency distribution (for example levels of income). A Gini coefficient of zero expresses perfect equality, where all values are the same (for example, where everyone has the same income). A Gini coefficient of one (or 100 percents ) expresses maximal inequality among values (for example where only one person has all the income).</p>
            <p>The Gini coefficient ranges from 0 to 1, with higher values indicating higher inequality. A value of 0 indicates perfect equality, while a value of 1 indicates perfect inequality.</p>
        """,
        "raref": """
            <p>Rarefaction is a technique to assess species richness from the results of sampling through extrapolation.</p>
        """,
    } -%}

    {% for divdir in job.out.outdir | glob: "*"  %}
        {% set method = divdir | basename %}
        <h{{h}}>Method: {{method}} ({{method_fullnames[method]}})</h{{h}}>

        {% if job.index == 0 %}
            <Tile>
                {{method_desc[method]}}
            </Tile>
        {% endif %}

        {% set imgs = [divdir + "/_nogrouping.png"] %}
        {% set caps = ["For individual samples"] %}
        {% for img in divdir | glob: "By *.png" %}
            {% set _ = imgs.append(img) %}
            {% set _ = caps.append(stem(img)) %}
        {% endfor %}

        {% if len(imgs) > 1 %}
            {{ table_of_images(imgs, caps) }}
        {% else %}
            <Image src={{imgs[0] | quote}} />
        {% endif %}

    {% endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.immdata | stem }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
