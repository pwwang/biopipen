---
name: pipen-process
description: How to write a new biopipen pipen Process subclass from scratch — the declarative class in biopipen/ns/*.py, the R/Python script template(s) in biopipen/scripts/*/, the reporter-based HTML report, and the local test class. Use this whenever the user asks to add a new process, wrap an external tool or analysis (Seurat, hdWGCNA, cellranger, scanpy, any R/Python function) as a pipeline step, or port an existing R/Python script into a reusable biopipen process — even if they don't explicitly say "process" (e.g. "add a step that runs X", "make X a pipen job", "I need X in scrna.py"). The worked example is the HdWGCNA process (biopipen/ns/scrna.py:991 with biopipen/scripts/scrna/HdWGCNA*.R).
---

# Writing a pipen Process

A biopipen process is three pieces:

1. **Declarative class** in `biopipen/ns/<domain>.py` — declares I/O, envs (parameters), the script file, and plugin options. No `run()`, no `__init__`, no `__outdata__`.
2. **Script template(s)** in `biopipen/scripts/<domain>/` — the actual R/Python code, with liquid `{{ ... }}` interpolation of inputs/envs. Split into include files when sections are optional and repeatable (`HdWGCNA-core.R`, `HdWGCNA-plots.R`, ...).
3. **Report structure** — built from the script via `reporter$add2(...)` calls; rendered by `reports/common.svelte`.

**Before writing anything, read `HdWGCNA` (scrna.py ~line 991) plus its scripts.** It is the reference implementation: a heavy process with metacell/pseudobulk/consensus modes, 10+ optional analysis sections, and per-case plots/tables. New processes copy its skeleton and drop what isn't needed.

## Step 1 — Understand the tool you are wrapping

The envs of the process are a thin declarative mirror of the external function's arguments, so first enumerate:

- The exact function signatures (R: `args(<fn>)`; Python: `inspect.signature`). Every argument that a user may reasonably tune becomes an env key.
- What the function returns and mutates: new Seurat/AnnData fields, files written, tables returned.
- Which arguments are **control vars** — your process's own plumbing, NOT the tool's args: `devpars`, `more_formats`, `mode`, `descr`, and any user-facing flags that select between tool calls. These live in the env dict but must never be passed into the tool (see `extract_vars` in Step 3).
- Whether the tool writes its own files (e.g. `ModuleNetworkPlot`, `MotifOverlapBarPlot`) — then you pass it an absolute `outdir` under `plots/` and reference the files by path in the report, never `save_plot` its return value.

### Env naming rules (user directive, follows hdWGCNA)

- **Env keys = the tool's argument names**, e.g. `wgcna_name`, `gene_select`, `group_by`, `n_hubs`. Underscores pass through the `todot="-"` renderer unchanged.
- **Dotted R args use kebab-case**: `group.by.vars` → env key `group-by-vars` (the renderer turns the dash back into the dot).
- **Env sections that configure a single tool function are named after it**: `SetupForWGCNA`, `MetacellsByGroups`, `ConstructNetwork`, `GetHubGenes`, ...
- **Section envs bundling several functions get descriptive names**: `plots`, `dmes`, `enrich`, `module_preservations`.
- When an env value is an expression or a range that must reach R as raw code (e.g. `dims: "1:15"`, `features: "VariableFeatures(seurat_obj)"`), it's a string — the `r` filter passes it through verbatim. That's intentional.

## Step 2 — Write the class

Template (annotated):

```python
class MyProcess(Proc):                      # subclass an existing Proc from the same ns file
    """Short title.

    Input:
        srtobj: The Seurat object file (rds/qs).     # name: type: description
    Output:
        outdir: Directory of the analysis results.
    Envs:
        envs.name (type=...):                       # envs.<name> in the script
            Description, what it does, default if non-obvious. One per env.
            Keep the type tag: (type=float), (type=bool), (type=dict), ...
        # dotted-arg envs documented with their R name: group-by-vars (type=...):
    Requires:
        - check: {{proc.lang}} <(echo "library(<pkg>)")   # one per R package
    """
    # noqa: E501

    input = "srtobj:file"
    output = "outdir:dir:{{in.srtobj | stem}}.hdwgcn"   # {{in.x | stem}} pattern
    lang = config.lang.rscript                            # or config.lang.python
    envs_depth = 4        # deep-merge envs across subclass overrides; use 4 when
                          # envs contain nested per-case dicts (plots.NAME.devpars.res
                          # must survive subclass overrides). Use 1 for flat envs.
    script = "file://../scripts/scrna/MyProcess.R"        # relative to the ns file
    plugin_opts = {
        "report": "file://../reports/common.svelte",
        "report_paging": 2,   # heavy sections: 1-2 h1s per page
    }
```

Docstring rules:

- Sections in this exact order: `Input:`, `Output:`, `Envs:`, `Requires:`. Lines over 88 chars get `# noqa: E501` at the docstring end.
- `Requires:` lines use the `{{proc.lang}} <(echo "library(<pkg>)")` check so the environment is validated before the job runs.
- Link each env to the wrapped function's reference page where one exists (e.g. the hdWGCNA docs).

Class rules:

- Declarative only. If you need the output path pattern from the input, use `{{in.x | stem}}` style liquid in `output`.
- `envs_depth`: 4 for deep merges (the HdWGCNA/SeuratPreparing pattern), 1 for flat env dicts. Too shallow → subclass overrides clobber whole nested dicts instead of merging.
- `mutaters`, `subset`, `seed`, `ncores` are conventional shared envs on heavy processes (see HdWGCNA).
- Never hardcode a cache/TOM path; resolve at runtime: cache dir = `envs.cache %||% job.outdir` or an explicit subdir under it.

## Step 3 — Write the script template

R skeleton (Python scripts follow the same structure with their own idioms):

```r
library(<tool>)            # every package used, fully loaded
library(biopipen.utils)    # do_call, extract_vars, save_plot, write_table, read_obj, ...

log <- get_logger()
reporter <- get_reporter()

infile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
envs.x <- {{envs.x | r}}   # each env into a top-level variable; case maps with
                           # {{envs.plots | r: todot="-", skip=1}} (skip the outer key)
group_by_vars <- {{envs["group-by-vars"] | r}}   # dashed keys need jinja2
                           # SUBSCRIPT syntax — `{{envs.group-by-vars | r}}`
                           # (dot notation) does not render (MetaMarkers.R:24)
```

- **`do_call` every external call** — `do_call(fn, c(list(<fixed args>), Args))` from biopipen.utils. Fixed args (the object being mutated, the name) go in the list; everything from the env goes in the args variable.
- **Strip control vars before EVERY do_call — plot calls included**: `case <- extract_vars(case, "devpars", "more_formats", "mode", ...)`. Passing `mode`/`devpars`/nested flag dicts into a tool that doesn't accept them → "unused argument" errors. The plot loop is where this gets skipped: `SeuratClusterStats-dimplots.R:16` nulls `case$devpars` before `do_call(CellDimPlot, case)` — copy that discipline, or the plot silently dies in its `tryCatch`.
- **Optional sections in `tryCatch`** → `log$warn(glue("Failed to run <section> {name}: {conditionMessage(e)}"))`. Core pipeline failures propagate (stop), optional analysis never kills the job.
- **Plots**: `save_plot(p, prefix, auto_devpars(case$devpars, npanels, ncol = 4), formats = c("png", more_formats))`. `auto_devpars` computes pixel dims from panel count (`per=2.5`); `save_plot_sidefx` for functions that plot themselves (base R/WGCNA plotters). Width/height are in PIXELS (res-aware); 800×600 @ res=100 is the floor. **`auto_devpars` and `hs_section` are script-local helpers defined in HdWGCNA.R (lines 87/96), NOT biopipen.utils exports — copy them into your script** (the codebase's own pattern), or the job dies with "could not find function".
- **Tables**: `write_table(df, file.path(tables_dir, "x.tsv"), row.names = TRUE)` for rowname-bearing data.
- **Report items** (rendered by common.svelte):

```r
reporter$add2(list(name = "Plot", contents = list(
    list(kind = "descr", content = descr),
    reporter$image(info$prefix, case$more_formats, FALSE, kind = "image"),
    list(kind = "table", src = tbl_path, data = list(nrows = 100))
)), hs = hs_section(section, name), ui = "tabs")   # ui="tabs": plot + table side by side
```

- `hs_section(section, name)` gives each case its own h1/h2; per-case `descr` env overrides auto text.
- **Introduction-first**: `reporter$report` is a public named list, insertion-ordered. To put an intro section first, build everything, then `reporter$report <- c(reporter$report["Introduction"], reporter$report[setdiff(names(reporter$report), "Introduction")])` before `reporter$save(joboutdir)`.
- Save the mutated object at the end: `save_obj(srtobj, file.path(outdir, "<stem>.qs"))`.

## Step 4 — Tests

`tests/test_<domain>/<name>/test.py`:

- `run.env` with `ENVNAME="biopipen"` and `LOCAL_ONLY="true"` (test skipped in CI; run locally).
- One test class per use case, each subclassing the process and overriding envs (see the 9 HdWGCNA classes — one per tutorial).
- `testing(pipen)` asserts outputs: `pipen.run()` truthiness, output files exist, `qs`/`rds` load, table columns present.
- **Never run concurrent test invocations** — they share `tests/running/biopipen-tests-1` and corrupt each other (two jobs writing the same outdir produce "cannot open the connection" style failures). One at a time, from the venv python: `.venv/bin/python tests/test_<domain>/<name>/test.py`.

## Pitfalls (hard-won)

- **R factor comparison bug**: `mod_table$module != "grey"` on a factor compares *level codes*, not labels — with "grey" as the last level every gene counts as non-grey, `n_modules` = total genes, and `auto_devpars` computes a >50000px device → "One or both dimensions exceed the maximum (50000px)". Always `as.character(...) != "grey"` (or `as.character()` before any `!=` on a module/group column).
- **Duplicate formal args in do_call**: if the args variable already contains a name the fixed list also passes (e.g. an env dict leaked `wgcna_name`), R errors `formal argument "x" matched by multiple actual arguments`. When composing a call, verify the args variable doesn't carry names you pass explicitly — set them to `NULL` before `do_call`.
- **Dashed env keys don't render via dot notation**: `{{envs.group-by-vars | r}}` fails to render (jinja2 parses `envs.s-features` as `envs.s - features`). Use subscript syntax `{{envs["group-by-vars"] | r}}` — precedent: `MetaMarkers.R:24`.
- **Control vars leak** → "unused argument" from tools with strict signatures. `extract_vars` before every call; hdWGCNA-style tools don't take `...`.
- **List-returning tools**: `PlotModulePreservation` (list → `patchwork::wrap_plots`), `GetModuleTraitCorrelation` (list(cor, pval, fdr) → one table per element). Check the return type before writing the save/report code.
- **Self-writing plotters** return invisible objects — pass an absolute `outdir`, don't `save_plot` the return.
- **Caching semantics**: editing any included script or the test file invalidates ALL cached jobs for the process. A killed/failed run leaves `job.rc != 0` → rerun from scratch. TOM-like big intermediates go to the cache dir with an `overwrite_tom=FALSE`-style default so reruns reuse them.
- **`envs_depth` too shallow** → nested per-case env overrides clobber instead of deep-merging; the failure shows up as missing keys in a subclass that "should" inherit.
- **Path mangling**: when a tool returns a path, it may prefix the repo root (observed: `repo//abs/path`). Fix/verify paths before writing files.
- **Long sections**: split the script into include files by responsibility so each stays readable; report paging (`report_paging`) keeps heavy pages manageable.
