from pathlib import Path
import pytest
from pyppl import PyPPL
from bioprocs.tumhet import pClonEvol, pPyClone
from . import remotedata

def test_clonevol():
	pClonEvol1 = pClonEvol.copy()
	pClonEvol1.input = (
		remotedata.get('tumhet/aml1.txt'), remotedata.get('tumhet/aml1.sample.txt'))
	PyPPL().start(pClonEvol1).run()

def test_pyclone():
	vcfdir = remotedata.get('tumhet/SRR385940-individuals/')
	pPyClone1 = pPyClone.copy()
	pPyClone1.input = (
		list(vcfdir.glob('*.vcf.gz')),
		list(vcfdir.glob('*.vcf.gz')))
	pPyClone1.report = """
{% python from os import path %}
## {{proc.desc}}

PyClone[1] is a tool using Probabilistic model for inferring clonal population structure from deep NGS sequencing.

![Similarity matrix]({{path.join(jobs[0].o.outdir, "plots/loci/similarity_matrix.svg")}})

```table
caption: Clusters
file: "{{path.join(jobs[0].o.outdir, "tables/cluster.tsv")}}"
rows: 10
```

[1]: Roth, Andrew, et al. "PyClone: statistical inference of clonal population structure in cancer." Nature methods 11.4 (2014): 396.
"""
	PyPPL().start(pPyClone1).run().report(Path(pPyClone1.workdir) / 'report.html', title = 'Clonality analysis using PyClone')
