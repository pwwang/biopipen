from unittest import TestCase, main

from pathlib import Path

from biopipen.core.filters import r
from biopipen.core.testing import r_test


def _get_size_rcode(
    fun: str,
    fun_args: dict,
    df_values: dict,
    expected: str,
    mutate: bool = False,
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

    if fun == "expanded+":
        fun = "expanded"
        fun_r_repr = f"{fun_r_repr}, include_emerged = TRUE"
    elif fun == "collapsed+":
        fun = "collapsed"
        fun_r_repr = f"{fun_r_repr}, include_vanished = TRUE"

    expected = ", ".join([f"'{x}'" for x in expected.split(" ")])
    expected = f"c({expected})"

    if mutate:
        rcode = f"""
            df <- tibble({df_r_repr})
            out <- df %>% mutate(out = {fun}(., {fun_r_repr})) %>% pull(out)
        """
    elif if_else:
        rcode = f"""
            df <- tibble({df_r_repr})
            out <- df %>%
                mutate(
                    out = if_else(CDR3.aa %in% {fun}(., {fun_r_repr}), "X", "Y")
                ) %>%
                pull(out)
        """
    else:
        rcode = f"""
            df <- tibble({df_r_repr})
            out <- {fun}(df, {fun_r_repr})
        """

    rcode = f"""
        {rcode}

        expected <- {expected}
        expected[expected == "NA"] <- NA
        expect(identical(out, expected), paste(out, collapse = " "))
    """

    return rcode


def _get_paired_rcode(args: dict, df_values: dict, expected: str) -> str:
    df_r_repr = ", ".join(
        [f"{col} = {r(values)}" for col, values in df_values.items()]
    )
    args_r_repr = ", ".join(
        [f"{arg} = {value}" for arg, value in args.items()]
    )

    expected = ", ".join([f"'{x}'" for x in expected.split(" ")])
    expected = f"c({expected})"

    return f"""
        df <- tibble({df_r_repr})
        out <- paired(df, {args_r_repr})

        expected <- {expected}
        expected[expected == "NA"] <- NA
        expect(identical(out, expected), paste(out, collapse = " "))
    """


class TestFiltersR(TestCase):

    SOURCE_FILE = Path(__file__).parent.parent.parent.parent.joinpath(
        "biopipen", "utils", "mutate_helpers.R"
    )

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

    paired_df_values = {
        "id": ["A", "A", "B", "B", "C", "C", "D", "D"],
        "compare": [1, 2, 1, 1, 1, 2, 1, 2],
    }

    @r_test
    def test_size1(self):
        df_values = self.df_values
        expected = "A D"
        fun = "expanded+"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size2(self):
        df_values = self.df_values
        expected = "A NA NA NA D NA NA NA NA NA NA NA"
        fun = "expanded+"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        uniq = 'FALSE'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
                "uniq": uniq,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size3(self):
        df_values = self.df_values
        expected = "A NA NA E D E E NA NA NA NA NA"
        fun = "expanded+"
        mutate = True
        group_by = "Source"
        idents = 'c("Tumor", "Normal")'
        uniq = 'FALSE'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "idents": idents,
                "uniq": uniq,
            },
            df_values=df_values,
            expected=expected,
            mutate=mutate,
        )

    @r_test
    def test_size4(self):
        df_values = self.df_values
        expected = "B E C"
        fun = "collapsed+"
        group_by = "Source"
        compare = "Clones"
        order = "desc(.diff)"
        idents = 'c("Tumor", "Normal")'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "order": order,
                "idents": idents,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size5(self):
        df_values = self.df_values
        expected = "NA C B E NA E E B B B NA NA"
        fun = "collapsed+"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        uniq = 'FALSE'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
                "uniq": uniq,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size6(self):
        df_values = self.df_values
        expected = "A D"
        fun = "emerged"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size7(self):
        df_values = self.df_values
        expected = "A NA NA NA D NA NA NA NA NA NA NA"
        fun = "emerged"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        uniq = 'FALSE'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
                "uniq": uniq,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size8(self):
        df_values = self.df_values
        expected = "B C"
        fun = "vanished"
        group_by = "Source"
        compare = "Clones"
        order = "desc(.diff)"
        idents = 'c("Tumor", "Normal")'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "order": order,
                "idents": idents,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size9(self):
        df_values = self.df_values
        expected = "NA C B NA NA NA NA B B B NA NA"
        fun = "vanished"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        uniq = 'FALSE'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
                "uniq": uniq,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size10(self):
        df_values = self.df_values
        expected = "B C"
        fun = "vanished"
        group_by = '"Source"'
        compare = "Clones"
        order = "desc(.diff)"
        idents = 'c("Tumor", "Normal")'
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "order": order,
                "idents": idents,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_size11(self):
        df_values = self.df_values
        expected = "Y X X Y Y Y Y X X X Y Y"
        fun = "vanished"
        group_by = "Source"
        compare = "Clones"
        idents = 'c("Tumor", "Normal")'
        uniq = 'FALSE'
        if_else = True
        return _get_size_rcode(
            fun=fun,
            fun_args={
                "group_by": group_by,
                "compare": compare,
                "idents": idents,
                "uniq": uniq,
            },
            df_values=df_values,
            expected=expected,
            if_else=if_else,
        )

    @r_test
    def test_paired1(self):
        df_values = self.paired_df_values
        expected = "A B C D"
        return _get_paired_rcode(
            args={
                "id": "id",
                "compare": "compare",
                "idents": 2,
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_paired2(self):
        df_values = self.paired_df_values
        expected = "A C D"
        return _get_paired_rcode(
            args={
                "id": "id",
                "compare": "compare",
                "idents": "1:2",
            },
            df_values=df_values,
            expected=expected,
        )

    @r_test
    def test_paired3(self):
        df_values = self.paired_df_values
        expected = "A A NA NA C C D D"
        return _get_paired_rcode(
            args={
                "id": "id",
                "compare": "compare",
                "idents": "1:2",
                "uniq": "FALSE",
            },
            df_values=df_values,
            expected=expected,
        )


# def run_paired(
#     df_values: dict,
#     expected: str,
#     **args: Any,
# ):
#     print(f">>> TESTING paired - {args}")
#     out = _run_paired_fun(args, df_values)
#     assert out == expected, f"{out} != {expected}"
#     print("    PASSED")
#     print("")


if __name__ == "__main__":
    main(verbosity=2)

    # # paired
    # paired_df_values = {
    #     "id": ["A", "A", "B", "B", "C", "C", "D", "D"],
    #     "compare": [1, 2, 1, 1, 1, 2, 1, 2],
    # }
    # run_paired(
    #     df_values=paired_df_values,
    #     expected="A B C D",
    #     id="id",
    #     compare="compare",
    #     idents=2,
    # )
    # run_paired(
    #     df_values=paired_df_values,
    #     expected="A C D",
    #     id="id",
    #     compare="compare",
    #     idents="1:2",
    # )
    # run_paired(
    #     df_values=paired_df_values,
    #     expected="A A NA NA C C D D",
    #     id="id",
    #     compare="compare",
    #     idents="1:2",
    #     uniq="FALSE",
    # )
