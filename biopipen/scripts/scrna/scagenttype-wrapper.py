"""Run scAgentType agent annotation on an AnnData and save per-cell labels.

Used by CellTypeAnnotation-scagenttype.R. The cluster memberships and the
AnnotationAgent configuration come from files; API credentials are expected
in the environment (injected by the R wrapper from `envs.scagenttype.api_key`
and `envs.scagenttype.base_url`, or exported in the job environment).
"""
from argparse import ArgumentParser
import json
import os
import shutil
import sys
import tempfile

# anndata.Settings (via scverse_misc) auto-loads a `.env` found from the
# current working directory at import time, with `extra="forbid"`. When the
# job's cwd chain contains an `.env` with e.g. OPENAI_* keys (like this
# test's), the import crashes with "Extra inputs are not permitted". Import
# from a clean cwd (no `.env` chain above /tmp), then restore.
_orig_cwd = os.getcwd()
os.chdir(tempfile.gettempdir())

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

os.chdir(_orig_cwd)

API_KEY_ENVS = {
    "openai": "OPENAI_API_KEY",
    "anthropic": "ANTHROPIC_API_KEY",
    "google": "GOOGLE_API_KEY",
}
DEFAULT_MODELS = {
    "openai": "gpt-4o-mini",
    "anthropic": "claude-sonnet-4-20250514",
}

parser = ArgumentParser(description="Run scAgentType")
parser.add_argument(
    "-i", "--input", required=True, help="Input H5AD file (AnnData)"
)
parser.add_argument(
    "-c", "--clusters", required=True,
    help="Cluster TSV file (barcode, cluster)"
)
parser.add_argument(
    "-o", "--output", required=True,
    help="Output TSV file for per-cell labels"
)
parser.add_argument(
    "--config", required=True,
    help="JSON config file for AnnotationAgent"
)
args = parser.parse_args()

with open(args.config) as f:
    cfg = json.load(f)

api = cfg.get("api", "openai")
key_env = API_KEY_ENVS.get(api)
if not key_env or not os.environ.get(key_env):
    sys.exit(
        f"API key for api '{api}' is not set. Set `envs.scagenttype.api_key` "
        f"or export the {key_env} environment variable."
    )
if not cfg.get("model"):
    model = DEFAULT_MODELS.get(api)
    if not model:
        sys.exit(f"No default model for api '{api}'; set `envs.scagenttype.model`")
    cfg["model"] = model

adata = sc.read_h5ad(args.input)
clusters = pd.read_csv(args.clusters, sep="\t", index_col=0)
# scAgentType requires the cluster_key column to be categorical
adata.obs["cluster"] = pd.Categorical(
    clusters["cluster"].reindex(adata.obs_names)
)
n_missing = int(adata.obs["cluster"].isna().sum())
if n_missing:
    sys.exit(
        f"{n_missing} h5ad barcodes are not found in the cluster table. "
        "The barcodes were likely changed during the h5ad conversion."
    )

from celltype_benchmark.agent import AnnotationAgent, LLMBackend  # noqa: E402

# Upstream crashes when the model puts prose strings in `actions` /
# `additional_actions` (only `response_format=json_object` is enforced, not
# the schema): the ReAct loop calls `a.get("tool")` on each item. Sanitize
# the parsed response — keep dict entries (executable tool calls), drop the
# rest with a warning, and let the loop synthesize from the evidence.
_complete_json_structured = LLMBackend.complete_json_structured


def _sanitized_complete_json_structured(self, messages, system):
    out = _complete_json_structured(self, messages, system)
    if not isinstance(out, dict):
        return out
    for key in ("actions", "additional_actions"):
        vals = out.get(key)
        if isinstance(vals, list) and any(not isinstance(v, dict) for v in vals):
            n = sum(not isinstance(v, dict) for v in vals)
            out[key] = [v for v in vals if isinstance(v, dict)]
            print(
                f"scAgentType: dropped {n} non-structured item(s) in "
                f"`{key}` of the model response (upstream expects dicts)",
                file=sys.stderr,
            )
    return out


LLMBackend.complete_json_structured = _sanitized_complete_json_structured

def make_agent():
    agent = AnnotationAgent(**{k: v for k, v in cfg.items() if v is not None})
    # Upstream default max_tokens=4096 truncates long responses on models
    # with heavy reasoning (observed empty content at the cap); give headroom.
    agent.llm.max_tokens = 8192
    return agent


agent = make_agent()
# The agent caches responses before parsing them, so a malformed response is
# replayed on retries; drop the response cache before each retry so the
# retried attempts are fresh. Also guards against a poisoned cache from a
# previous failed run in the same job dir.
result = None
for attempt in range(1, 4):
    try:
        result = agent.annotate(adata, cluster_key="cluster")
        break
    except Exception:
        if attempt == 3:
            raise
        print(
            f"scAgentType annotate() failed (attempt {attempt}/3); "
            "clearing the response cache and retrying ...",
            file=sys.stderr,
        )
        shutil.rmtree(agent.llm._cache_dir, ignore_errors=True)
        agent = make_agent()
assert result is not None
pd.DataFrame(
    {
        "barcode": adata.obs_names,
        "scagenttype_celltype": np.asarray(result["labels"]),
    }
).to_csv(args.output, sep="\t", index=False)
