# {{report.title}}

MutSig[[1]] stands for "Mutation Significance".
MutSig analyzes lists of mutations discovered in DNA sequencing,
to identify genes that were mutated more often than expected by chance
given background mutation processes.

The idea could be simply elaborated as:

![MutSig Illustration](http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/mutsig_fig1.PNG)

{% for job in jobs %}
:::: {.panel}
## {{job.i.infile | fn}}

```table
caption: Top 100 significantly-mutated genes
file: {{job.o.outdir, '*.sig_genes.txt' | *glob1}}
download: true
```

{%  if job.index == 0 %}
{#  only show the tip once #}
::: {.admon .tip}
You can download the full list of genes using the `download` button on the right corner of the table.
:::

::: {.admon .info head="Description for each column"}
See more details in the paper's [supplementary documentation](https://static-content.springer.com/esm/art%3A10.1038%2Fnature12213/MediaObjects/41586_2013_BFnature12213_MOESM255_ESM.pdf)

- `gene` = name of the gene that this line reports coverage for
- `expr` = expression level of this gene, averaged across many cell lines in
the Cancer Cell Line Encylcopedia
- `reptime` = DNA replication time of this gene, ranging approximately
from 100 (very early) to 1000 (very late)
- `hic` = chromatin compartment of this gene, measured from HiC experment, ranging approximately from -50 (very closed) to +50 (very
open)
- `N_nonsilent` = Coverage of non-silent mutations
- `N_silent` = Coverage of silent mutations
- `N_noncoding` = Coverage of noncoding region mutations
- `n_nonsilent` = Number of non-silent mutations
- `n_silent` = Number of silent mutations
- `n_noncoding` = Number of noncoding region mutations
- `nnei` = Number of neighboring genes that are pooled together to compute the background mutation rate for that gene
- `x` = Number of mutated bases in these neigboring genes that are either silent or non-coding
- `X` = Total number of bases related to these neighboring genes
- `p` = P value
- `q` = Q value
:::
{%  endif %}


::::
{% endfor %}


[[1]]: Lawrence, Michael S., et al. "Mutational heterogeneity in cancer and the search for new cancer-associated genes." Nature 499.7457 (2013): 214-218.
