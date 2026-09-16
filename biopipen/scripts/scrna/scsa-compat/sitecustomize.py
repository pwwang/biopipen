# Compatibility shims for SCSA.py, shipped with biopipen. `scsa-wrapper.py`
# prepends this directory to the PYTHONPATH of the SCSA process, so python
# imports this `sitecustomize` at startup and patches only that interpreter.
# SCSA.py is a 2019-era script: it uses numpy APIs removed in numpy 2.0
# (`asfarray`, `mat`) and unpickles whole.db, which holds pandas<1.0 objects
# whose classes (`pandas.core.indexes.numeric.{Int,Float,UInt}64Index`) modern
# pandas removed.
import sys, types
import numpy as np
import pandas as pd

if not hasattr(np, "asfarray"):
    np.asfarray = lambda a, dtype=float: np.asarray(a, dtype=dtype)
if not hasattr(np, "mat"):
    np.mat = np.matrix

_numeric = types.ModuleType("pandas.core.indexes.numeric")
for _name in ("Int64Index", "Float64Index", "UInt64Index", "NumericIndex"):
    setattr(_numeric, _name, pd.Index)
sys.modules.setdefault("pandas.core.indexes.numeric", _numeric)

# pandas 2.0 removed DataFrame.append; SCSA.py's GO annotation uses it
if not hasattr(pd.DataFrame, "append"):
    pd.DataFrame.append = (
        lambda self, other, ignore_index=False, **kw: pd.concat(
            [self, other], ignore_index=ignore_index, sort=False
        )
    )
