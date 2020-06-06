import random
import math
import cmdy
from pathlib import Path
from diot import Diot

infile    = Path({{i.infile | quote}})
trainfile = Path({{o.trainfile | quote}})
testfile  = Path({{o.testfile | quote}})
outdir    = Path({{job.outdir | quote}})
seed      = {{args.seed | ? | ! :None | =:_ | repr }}
inopts    = {{args.inopts | repr}}
method    = {{args.method | quote}}

random.seed(seed)

if inopts.cnames:
	cmdy.head(infile, n=1).r > trainfile
	cmdy.head(infile, n=1).r > testfile
	cmdy.tail(infile, n='+2').r > outdir.joinpath(infile.with_suffix('.main'))
	infile = outdir.joinpath(infile.with_suffix('.main'))

# get total lines
nlines = int(cmdy.wc(l=infile).split()[0])

# get the number of tests
if method.endswith('-fold'):
	fold = int(method[:-5])
	ntest = math.ceil(float(nlines) / float(fold))
elif method == 'leave-one':
	ntest = 1
else:
	ntrain, ntest = method.split(':', 1)
	ntest = math.ceil(int(ntest) * nlines / (int(ntrain) + int(ntest)))

tests = random.sample(list(range(nlines)), ntest)

with open(infile) as fin, open(trainfile, "a") as ftrain, open(testfile, "a") as ftest:
	for i, line in enumerate(fin):
		if i in tests:
			ftest.write(line)
		else:
			ftrain.write(line)

if infile.suffix == '.main':
	infile.unlink()
