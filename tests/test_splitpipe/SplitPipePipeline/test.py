"""Test the wiring of the SplitPipePipeline ProcGroup

This test does not run `split-pipe` and does not need sequencing data, so it
can run in CI. It builds the process classes of the group and checks the
generated graph, the inputs and the environment variables.

Run from the repo root:

    .venv/bin/python tests/test_splitpipe/SplitPipePipeline/test.py
"""
from pathlib import Path

import pandas as pd
from pipen.template import get_template_engine
from pipen_args.parser_ import Parser

from biopipen.core.proc import Proc
from biopipen.ns.splitpipe_pipeline import SplitPipePipeline

REPO_ROOT = Path(__file__).parents[3]
ENGINE = get_template_engine("liquid")

GEX = dict(
    fq1=[["s1_R1.fastq.gz"], ["s2_R1.fastq.gz"]],
    fq2=[["s1_R2.fastq.gz"], ["s2_R2.fastq.gz"]],
    ids=["s1", "s2"],
    genome_dir="/genomes/test",
)

# The GEX outputs as generated before the tcr/bcr blocks were added
GEX_OUTPUTS = {
    "SplitPipeCount": (
        "outdir:dir:\n"
        "        {%- if in.id -%}\n"
        "            {{in.id}}\n"
        "        {%- else -%}\n"
        '            {%- set id = in.fq1 | first | basename |\n'
        '                regex_replace: "(_L\\d+)?(_R\\d+)?(_\\d+)?\\.f(ast)?q'
        '(\\.gz)?$", "" |\n'
        '                regex_replace: "_S\\d+$", "" -%}\n'
        "            {{- id -}}\n"
        "        {%- endif -%}\n"
        "    "
    ),
    "SplitPipeCombine": (
        "outdir:dir:\n"
        "        {%- if in.id -%}\n"
        "            {{in.id}}\n"
        "        {%- else -%}\n"
        '            {{- in.sublibs | first | basename | append: "-comb" -}}\n'
        "        {%- endif -%}\n"
        "    "
    ),
    "SplitPipeSummary": (
        "outdir:dir:{{in.indirs | first | stem | append: '-etc.splitpipe_summary'}}"
    ),
}


def build(**opts):
    """Build a fresh group instance

    `SplitPipePipeline` is a singleton ProcGroup and pipen_args uses a
    singleton argument parser, so both need to be reset for each case.
    """
    Parser._INST = None
    SplitPipePipeline._INST = None
    opts_ = dict(GEX)
    opts_.update(opts)
    return SplitPipePipeline(**opts_)


def render(template, data):
    """Render a Liquid template with biopipen's filters"""
    return ENGINE(template, **Proc.template_opts).render(data)


def expect_error(exc_type, fragment, **opts):
    """Expect the group construction to fail with a message fragment"""
    try:
        build(**opts)
    except exc_type as exc:
        assert fragment in str(exc), f"{fragment!r} not in {exc}"
        print(f"PASS: {exc_type.__name__}: {exc}")
    else:
        raise AssertionError(f"Expected {exc_type.__name__} containing {fragment!r}")


def test_no_immune():
    """Without tcr/bcr, the graph is the same as before"""
    group = build()
    assert sorted(group.procs) == [
        "SplitPipeCombine",
        "SplitPipeCount",
        "SplitPipeSummary",
    ], sorted(group.procs)
    for name, output in GEX_OUTPUTS.items():
        assert group.procs[name].output == output, name
    assert group.p_splitpipe_summary.requires == [group.p_splitpipe_combine]
    print("PASS: no tcr/bcr -> GEX proc names and outputs unchanged")


def test_tcr():
    """A tcr block builds a paired count/combine pair"""
    group = build(
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"], ["t2_R1.fastq.gz"]],
            fq2=[["t1_R2.fastq.gz"], ["t2_R2.fastq.gz"]],
            ids=["s1", "s2"],
        )
    )
    assert "SplitPipeTCRCount" in group.procs
    assert "SplitPipeTCRCombine" in group.procs

    count = group.p_splitpipe_tcr_count
    assert count.name == "SplitPipeTCRCount"
    assert count.input == "fq1:files, fq2:files, id, genome_dir:dir, parent_dir:dir"
    assert count.envs["tcr_analysis"] is True
    assert count.envs["bcr_analysis"] is False
    assert count.envs["immune_genome"] == "human"
    assert count.envs["use_imgt_db"] is False
    assert render(count.output, {"in": {"id": "s1"}}) == "outdir:dir:s1-tcr"

    gex_ch = pd.DataFrame(
        {"outdir": ["/out/SplitPipeCount/s1", "/out/SplitPipeCount/s2"]}
    )
    assert count.input_data(gex_ch) == [
        (
            ["t1_R1.fastq.gz"],
            ["t1_R2.fastq.gz"],
            "s1",
            "/genomes/test",
            "/out/SplitPipeCount/s1",
        ),
        (
            ["t2_R1.fastq.gz"],
            ["t2_R2.fastq.gz"],
            "s2",
            "/genomes/test",
            "/out/SplitPipeCount/s2",
        ),
    ], count.input_data(gex_ch)
    print("PASS: tcr count -> envs, output and GEX pairing")

    try:
        count.input_data(pd.DataFrame({"outdir": ["/out/SplitPipeCount/s2"]}))
    except ValueError as exc:
        assert "must equal its GEX sublibrary id" in str(exc), exc
    else:
        raise AssertionError("Unmatched tcr id did not raise")
    print("PASS: tcr count -> unmatched id raises at runtime")

    combine = group.p_splitpipe_tcr_combine
    assert combine.name == "SplitPipeTCRCombine"
    assert combine.requires == [count, group.p_splitpipe_combine]
    assert combine.envs["immune_genome"] == "human"
    assert combine.envs["use_imgt_db"] is False

    immune_ch = pd.DataFrame(
        {
            "outdir": [
                "/out/SplitPipeTCRCount/s1-tcr",
                "/out/SplitPipeTCRCount/s2-tcr",
            ]
        }
    )
    try:
        combine.input_data(immune_ch, pd.DataFrame({"outdir": [None]}))
    except ValueError as exc:
        assert "combined GEX output directory is missing" in str(exc), exc
    else:
        raise AssertionError("Missing GEX combine dir did not raise")

    gex_comb_ch = pd.DataFrame({"outdir": ["/out/SplitPipeCombine/s1-comb"]})
    assert combine.input_data(immune_ch, gex_comb_ch) == [
        (
            ["/out/SplitPipeTCRCount/s1-tcr", "/out/SplitPipeTCRCount/s2-tcr"],
            None,
            "/out/SplitPipeCombine/s1-comb",
        )
    ]
    assert combine.input == "sublibs:dirs, id, parent_dir:dir"
    assert (
        render(
            combine.output,
            {"in": {"id": None, "sublibs": ["/out/SplitPipeTCRCount/s1-tcr"]}},
        )
        == "outdir:dir:s1-tcr-comb"
    )
    print("PASS: tcr combine -> parent_dir, output and immune sublibrary dirs")

    summary = group.p_splitpipe_summary
    assert summary.requires == [group.p_splitpipe_combine, combine]
    assert summary.input_data(
        pd.DataFrame({"outdir": ["/out/SplitPipeCombine/s1-comb"]}),
        pd.DataFrame({"outdir": ["/out/SplitPipeTCRCombine/s1-tcr-comb"]}),
    ) == [["/out/SplitPipeCombine/s1-comb", "/out/SplitPipeTCRCombine/s1-tcr-comb"]]
    print("PASS: summary -> one job over the GEX and tcr final dirs")

    script = (REPO_ROOT / "biopipen/scripts/splitpipe/SplitPipeCount.py").read_text()
    rendered = render(
        script,
        {
            "in": {
                "fq1": ["/data/t1_R1.fastq.gz"],
                "fq2": None,
                "id": "s1",
                "genome_dir": "/genomes/test",
                "parent_dir": "/out/SplitPipeCount/s1",
            },
            "out": {"outdir": "s1-tcr"},
            "envs": count.envs,
            "job": {},
            "proc": None,
        },
    )
    assert "tcr_analysis: bool = True" in rendered
    assert "immune_genome = 'human'" in rendered
    assert 'command += ["--tcr_analysis"]' in rendered
    assert 'command += ["--immune_genome", immune_genome]' in rendered
    # the per-sublibrary parent arrives through the input channel, so it must
    # survive into the rendered command (envs alone would drop it)
    # a `dir`-typed input arrives as xqute's MountedLocalPath; the rendered line must
    # carry the plain path, never the class expression (which nothing imports)
    parent_line = next(
        line for line in rendered.splitlines() if line.startswith("parent_dir = ")
    )
    assert "/out/SplitPipeCount/s1" in parent_line, parent_line
    assert "MountedLocalPath" not in rendered, parent_line
    assert 'command += ["--parent_dir", str(parent_dir)]' in rendered

    combine_script = (
        REPO_ROOT / "biopipen/scripts/splitpipe/SplitPipeCombine.py"
    ).read_text()
    rendered_combine = render(
        combine_script,
        {
            "in": {
                "sublibs": ["/out/SplitPipeTCRCount/s1-tcr"],
                "id": None,
                "parent_dir": "/out/SplitPipeCombine/s1-comb",
            },
            "out": {"outdir": "s1-tcr-comb"},
            "envs": combine.envs,
            "job": {},
            "proc": None,
        },
    )
    assert "parent_dir = '/out/SplitPipeCombine/s1-comb'" in rendered_combine
    assert 'command += ["--parent_dir", str(parent_dir)]' in rendered_combine
    print("PASS: tcr flags render into the split-pipe command (incl. --parent_dir)")


def test_bcr():
    """A bcr block supports the transgenic_mouse immune genome"""
    group = build(
        bcr=dict(
            fq1=[["b1_R1.fastq.gz"]],
            ids=["s1"],
            immune_genome="transgenic_mouse",
        )
    )
    count = group.p_splitpipe_bcr_count
    assert count.name == "SplitPipeBCRCount"
    assert count.envs["bcr_analysis"] is True
    assert count.envs["tcr_analysis"] is False
    assert count.envs["immune_genome"] == "transgenic_mouse"
    assert render(count.output, {"in": {"id": "s1"}}) == "outdir:dir:s1-bcr"
    print("PASS: bcr block -> envs and output")


def test_explicit_parent_dir():
    """An explicit parent_dir skips id pairing and the GEX combine check"""
    group = build(
        combine=False,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"]],
            ids=["not-a-gex-id"],
            parent_dir="/parents/gex-s1",
        ),
    )
    count = group.p_splitpipe_tcr_count
    combine = group.p_splitpipe_tcr_combine
    assert count.envs["parent_dir"] == "/parents/gex-s1"
    assert count.input_data(pd.DataFrame()) == [
        (
            ["t1_R1.fastq.gz"],
            None,
            "not-a-gex-id",
            "/genomes/test",
            "/parents/gex-s1",
        )
    ]
    assert combine.envs["parent_dir"] == "/parents/gex-s1"
    assert combine.requires == [count]
    print("PASS: explicit tcr.parent_dir -> ids unpaired, parent used verbatim")


def test_existing_gex_parent_lists():
    """An empty GEX branch uses one existing parent per immune id"""
    group = build(
        fq1=[],
        combine=False,
        summary=False,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"], ["t2_R1.fastq.gz"]],
            ids=["s1", "s2"],
            parent_dir=["/parents/s1", "/parents/s2"],
            combine=False,
        ),
    )
    assert "SplitPipeTCRCount" in group.procs
    assert "SplitPipeCount" not in group.procs
    assert "SplitPipeCombine" not in group.procs
    count = group.p_splitpipe_tcr_count
    assert count.envs["parent_dir"] is None
    assert count.input_data == [
        (["t1_R1.fastq.gz"], None, "s1", "/genomes/test", "/parents/s1"),
        (
            ["t2_R1.fastq.gz"],
            None,
            "s2",
            "/genomes/test",
            "/parents/s2",
        ),
    ]
    print("PASS: empty GEX branch -> immune counts use per-id parent dirs")


def test_immune_only_with_gex_combine():
    """An immune-only run does not build a GEX combine process"""
    group = build(
        fq1=[],
        combine=True,
        summary=True,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"]],
            ids=["s1"],
            parent_dir="/parents/s1",
            parent_dir_comb="/combined/s1",
        ),
    )
    assert "SplitPipeCombine" not in group.procs
    assert group.p_splitpipe_combine is None
    assert group.p_splitpipe_summary.requires == [group.p_splitpipe_tcr_combine]
    print("PASS: immune-only run -> no GEX count/combine, summary over immune")


def test_per_id_parent_dir_comb():
    """A combined GEX parent list creates one immune combine row per id"""
    group = build(
        fq1=[],
        combine=False,
        summary=False,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"], ["t2_R1.fastq.gz"]],
            ids=["s1", "s2"],
            parent_dir=["/parents/s1", "/parents/s2"],
            parent_dir_comb=["/combined/s1", "/combined/s2"],
        ),
    )
    count = group.p_splitpipe_tcr_count
    combine = group.p_splitpipe_tcr_combine
    immune_ch = pd.DataFrame(
        {"outdir": ["/out/s1-tcr", "/out/s2-tcr"]}
    )
    assert combine.input_data(immune_ch) == [
        (["/out/s1-tcr"], None, "/combined/s1"),
        (["/out/s2-tcr"], None, "/combined/s2"),
    ]
    assert combine.envs["parent_dir"] is None
    assert combine.requires == [count]
    print("PASS: list parent_dir_comb -> one combine row per id")


def test_immune_parent_list_validation():
    """Parent directory lists must match the immune ids"""
    expect_error(
        ValueError,
        "2 ids but 1 parent_dir",
        fq1=[],
        combine=False,
        summary=False,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"], ["t2_R1.fastq.gz"]],
            ids=["s1", "s2"],
            parent_dir=["/parents/s1"],
            combine=False,
        ),
    )
    expect_error(
        ValueError,
        "2 ids but 1 parent_dir_comb",
        fq1=[],
        combine=False,
        summary=False,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"], ["t2_R1.fastq.gz"]],
            ids=["s1", "s2"],
            parent_dir=["/parents/s1", "/parents/s2"],
            parent_dir_comb=["/combined/s1"],
        ),
    )
    expect_error(
        ValueError,
        "parent_dir_comb",
        fq1=[],
        combine=False,
        summary=False,
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"], ["t2_R1.fastq.gz"]],
            ids=["s1", "s2"],
            parent_dir=["/parents/s1", "/parents/s2"],
        ),
    )


def test_empty_fq1_without_immune():
    """An empty top-level fq1 still requires an immune parent branch"""
    expect_error(
        TypeError,
        "empty `fq1` requires every immune block to set `parent_dir`",
        fq1=[],
    )


def test_immune_combine_disabled():
    """combine = false on a block summarizes the immune count dirs"""
    group = build(
        tcr=dict(
            fq1=[["t1_R1.fastq.gz"]],
            ids=["s1"],
            combine=False,
        )
    )
    assert group.p_splitpipe_tcr_combine is None
    assert group.p_splitpipe_summary.requires == [
        group.p_splitpipe_combine,
        group.p_splitpipe_tcr_count,
    ]
    print("PASS: tcr.combine = false -> no tcr combine, summary uses count dirs")


def test_validation():
    """Invalid tcr/bcr blocks are rejected with actionable messages"""
    expect_error(
        ValueError,
        "must equal its GEX sublibrary id",
        tcr=dict(fq1=[["t1_R1.fastq.gz"]], ids=["t3"]),
    )
    expect_error(
        ValueError,
        "immune_genome",
        tcr=dict(fq1=[["t1_R1.fastq.gz"]], ids=["s1"], immune_genome="rat"),
    )
    expect_error(
        ValueError,
        "requires `combine = true`",
        combine=False,
        tcr=dict(fq1=[["t1_R1.fastq.gz"]], ids=["s1"]),
    )
    expect_error(
        TypeError,
        "non-empty list of lists",
        tcr=dict(fq1=[], ids=[]),
    )
    expect_error(
        TypeError,
        "same length as `tcr.fq1`",
        tcr=dict(fq1=[["t1_R1.fastq.gz"]], ids=["s1", "s2"]),
    )
    expect_error(
        TypeError,
        "`ids` of `SplitPipePipeline` is required",
        ids=None,
        tcr=dict(fq1=[["t1_R1.fastq.gz"]], ids=["s1"]),
    )


if __name__ == "__main__":
    test_no_immune()
    test_tcr()
    test_bcr()
    test_explicit_parent_dir()
    test_existing_gex_parent_lists()
    test_immune_only_with_gex_combine()
    test_per_id_parent_dir_comb()
    test_immune_parent_list_validation()
    test_empty_fq1_without_immune()
    test_immune_combine_disabled()
    test_validation()
    print("ALL PASS")
