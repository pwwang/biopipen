from typing import Any

import subprocess as sp
from pathlib import Path

from biopipen.core.filters import r


R_FILE = Path(__file__).parent.parent.parent.parent.joinpath(
    "biopipen", "utils", "mutate_helpers.R"
)


def run_fun(
    fun: str,
    fun_args: dict,
    df_values: dict,
    mutate: bool,
    if_else: bool = False,
):
    df_r_repr = ", ".join(
        [f"{col} = {r(values)}" for col, values in df_values.items()]
    )
    fun_r_repr = ", ".join(
        [
            f"{arg.replace('_', '.')} = {value}"
            for arg, value in fun_args.items()
        ]
    )
    if mutate:
        rcode = f"""
            source("{R_FILE}")

            df <- tibble({df_r_repr})
            out <- df %>% mutate(out = {fun}(., {fun_r_repr})) %>% pull(out)
            cat(out)
        """
    elif if_else:
        rcode = f"""
            source("{R_FILE}")

            df <- tibble({df_r_repr})
            out <- df %>%
                mutate(
                    out = if_else(CDR3.aa %in% {fun}(., {fun_r_repr}), "X", "Y")
                ) %>%
                pull(out)
            cat(out)
        """
    else:
        rcode = f"""
            source("{R_FILE}")

            df <- tibble({df_r_repr})
            out <- {fun}(df, {fun_r_repr})
            cat(out)
        """

    out = sp.check_output(["Rscript", "-e", rcode])
    return out.decode().strip()


def run(
    df_values: dict,
    expected: str,
    fun: str,
    mutate: bool = False,
    if_else: bool = False,
    **fun_args: Any,
):
    print(f">>> TESTING {fun} - {fun_args}")
    out = run_fun(fun, fun_args, df_values, mutate=mutate, if_else=if_else)
    assert out == expected, f"{out} != {expected}"
    print("    PASSED")
    print("")


if __name__ == "__main__":
    df_values = {
        "Clones": [10, 8, 1, 5, 9, 2, 3, 7, 6, 4, 9, 9],
        "Source": [
            "Tumor", "Normal", "Normal", "Normal", "Tumor", "Tumor",
            "Tumor", "Normal", "Normal", "Normal", "NA", "X"
        ],
        "CDR3.aa": [
            "A", "C", "B", "E", "D", "E", "E", "B", "B", "B", "A", "A"
        ]
    }
    run(
        df_values=df_values,
        expected="A D",
        fun="expanded",
        group_by="Source",
        idents='c("Tumor", "Normal")',
    )
    run(
        df_values=df_values,
        expected="A NA NA NA D NA NA NA NA NA A A",
        fun="expanded",
        group_by="Source",
        idents='c("Tumor", "Normal")',
        uniq='FALSE',
    )
    run(
        df_values=df_values,
        expected="A NA NA NA D NA NA NA NA NA A A",
        fun="expanded",
        mutate=True,
        group_by="Source",
        idents='c("Tumor", "Normal")',
        uniq='FALSE',
    )
    run(
        df_values=df_values,
        expected="B E C",
        fun="collapsed",
        group_by="Source",
        idents='c("Tumor", "Normal")',
    )
    run(
        df_values=df_values,
        expected="NA C B E NA E E B B B NA NA",
        fun="collapsed",
        group_by="Source",
        idents='c("Tumor", "Normal")',
        uniq='FALSE',
    )
    run(
        df_values=df_values,
        expected="A D",
        fun="emerged",
        group_by="Source",
        idents='c("Tumor", "Normal")',
    )
    run(
        df_values=df_values,
        expected="A NA NA NA D NA NA NA NA NA A A",
        fun="emerged",
        group_by="Source",
        idents='c("Tumor", "Normal")',
        uniq='FALSE',
    )
    run(
        df_values=df_values,
        expected="B C",
        fun="vanished",
        group_by="Source",
        idents='c("Tumor", "Normal")',
    )
    run(
        df_values=df_values,
        expected="NA C B NA NA NA NA B B B NA NA",
        fun="vanished",
        group_by="Source",
        idents='c("Tumor", "Normal")',
        uniq='FALSE',
    )
    # what if quoted group-by
    run(
        df_values=df_values,
        expected="B C",
        fun="vanished",
        group_by='"Source"',
        idents='c("Tumor", "Normal")',
    )
    run(
        df_values=df_values,
        expected="Y X X Y Y Y Y X X X Y Y",
        fun="vanished",
        group_by="Source",
        idents='c("Tumor", "Normal")',
        uniq='FALSE',
        if_else=True,
    )
