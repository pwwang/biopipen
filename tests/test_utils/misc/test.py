from pathlib import Path
from unittest import TestCase, main
from biopipen.core.testing import r_test


class TestUtilsMisc(TestCase):

    SOURCE_FILE = Path(__file__).parent.parent.parent.parent.joinpath(
        "biopipen", "utils", "misc.R"
    )

    @r_test
    def test_expand_cases(self):
        return """
            defaults <- list(a = 0, b = 0, x = list(c = 0, d = 0))
            cases <- list(
                case1 = list(a = 1),
                case2 = list(b = 2),
                case3 = list(x = list(d = 3)),
                case4 = list()
            )
            expand_each <- function(name, case) {{
                if (name == "case4") {{
                    return(list("case4_1" = case, "case4_2" = case))
                }} else {{
                    out <- list(case)
                    names(out) <- name
                    return(out)
                }}
            }}
            excases <- expand_cases(cases, defaults, expand_each)

            expect(identical(excases$case1$a, 1), "excases$case1$a ==", excases$case1$a)
            expect(identical(excases$case1$b, 0), "excases$case1$b ==", excases$case1$b)
            expect(identical(excases$case1$x$c, 0), "excases$case1$x$c ==", excases$case1$x$c)
            expect(identical(excases$case1$x$d, 0), "excases$case1$x$d ==", excases$case1$x$d)

            expect(identical(excases$case2$a, 0), "excases$case2$a ==", excases$case2$a)
            expect(identical(excases$case2$b, 2), "excases$case2$b ==", excases$case2$b)
            expect(identical(excases$case2$x$c, 0), "excases$case2$x$c ==", excases$case2$x$c)
            expect(identical(excases$case2$x$d, 0), "excases$case2$x$d ==", excases$case2$x$d)

            expect(identical(excases$case3$a, 0), "excases$case3$a ==", excases$case3$a)
            expect(identical(excases$case3$b, 0), "excases$case3$b ==", excases$case3$b)
            expect(identical(excases$case3$x$c, 0), "excases$case3$x$c ==", excases$case3$x$c)
            expect(identical(excases$case3$x$d, 3), "excases$case3$x$d ==", excases$case3$x$d)

            expect(identical(excases$case4_1$a, 0), "excases$case4_1$a ==", excases$case4_1$a)
            expect(identical(excases$case4_1$b, 0), "excases$case4_1$b ==", excases$case4_1$b)
            expect(identical(excases$case4_1$x$c, 0), "excases$case4_1$x$c ==", excases$case4_1$x$c)
            expect(identical(excases$case4_1$x$d, 0), "excases$case4_1$x$d ==", excases$case4_1$x$d)

            expect(identical(excases$case4_2$a, 0), "excases$case4_2$a ==", excases$case4_2$a)
            expect(identical(excases$case4_2$b, 0), "excases$case4_2$b ==", excases$case4_2$b)
            expect(identical(excases$case4_2$x$c, 0), "excases$case4_2$x$c ==", excases$case4_2$x$c)
            expect(identical(excases$case4_2$x$d, 0), "excases$case4_2$x$d ==", excases$case4_2$x$d)
        """

    @r_test
    def test_casename_info_no_sections(self):
        return """
            cases <- list(
                case1 = list(a = 1),
                case2 = list(b = 2)
            )
            info <- casename_info("case1", cases, "/tmp", create = FALSE)

            expect(identical(info$casename, "case1"), "info$casename ==", info$casename)
            expect(identical(info$case, "case1"), "info$case ==", info$case)
            expect(identical(info$case_slug, "case1"), "info$case_slug ==", info$case_slug)
            expect(identical(info$section, "DEFAULT"), "info$section ==", info$section)
            expect(identical(info$section_slug, "DEFAULT"), "info$section_slug ==", info$section_slug)
            expect(identical(info$h1, "case1"), "info$h1 ==", info$h1)
            expect(identical(info$h2, "#"), "info$h2 ==", info$h2)
            expect(identical(info$casedir, "/tmp/case1"), "info$casedir ==", info$casedir)
        """

    @r_test
    def test_casename_info_mixed_section(self):
        return """
            cases <- list(
                case1 = list(a = 1, section = "section1"),
                case2 = list(b = 2)
            )

            info1 <- casename_info("section1::case1", cases, "/tmp", create = FALSE)

            expect(identical(info1$casename, "section1::case1"), "info1$casename ==", info1$casename)
            expect(identical(info1$case, "case1"), "info1$case ==", info1$case)
            expect(identical(info1$case_slug, "case1"), "info1$case_slug ==", info1$case_slug)
            expect(identical(info1$section, "section1"), "info1$section ==", info1$section)
            expect(identical(info1$section_slug, "section1"), "info1$section_slug ==", info1$section_slug)
            expect(identical(info1$h1, "section1"), "info1$h1 ==", info1$h1)
            expect(identical(info1$h2, "case1"), "info1$h2 ==", info1$h2)
            expect(identical(info1$casedir, "/tmp/section1/case1"), "info1$casedir ==", info1$casedir)

            info2 <- casename_info("case2", cases, "/tmp", create = FALSE)

            expect(identical(info2$casename, "case2"), "info2$casename ==", info2$casename)
            expect(identical(info2$case, "case2"), "info2$case ==", info2$case)
            expect(identical(info2$case_slug, "case2"), "info2$case_slug ==", info2$case_slug)
            expect(identical(info2$section, "DEFAULT"), "info2$section ==", info2$section)
            expect(identical(info2$section_slug, "DEFAULT"), "info2$section_slug ==", info2$section_slug)
            expect(identical(info2$h1, "DEFAULT"), "info2$h1 ==", info2$h1)
            expect(identical(info2$h2, "case2"), "info2$h2 ==", info2$h2)
            expect(identical(info2$casedir, "/tmp/DEFAULT/case2"), "info2$casedir ==", info2$casedir)
        """

    @r_test
    def test_casename_info_single_section(self):
        return """
            cases <- list(
                case1 = list(a = 1, section = "section1"),
                case2 = list(b = 2, section = "section1")
            )

            info1 <- casename_info("section1::case1", cases, "/tmp", create = FALSE)

            expect(identical(info1$casename, "section1::case1"), "info1$casename ==", info1$casename)
            expect(identical(info1$case, "case1"), "info1$case ==", info1$case)
            expect(identical(info1$case_slug, "case1"), "info1$case_slug ==", info1$case_slug)
            expect(identical(info1$section, "section1"), "info1$section ==", info1$section)
            expect(identical(info1$section_slug, "section1"), "info1$section_slug ==", info1$section_slug)
            expect(identical(info1$h1, "section1: case1"), "info1$h1 ==", info1$h1)
            expect(identical(info1$h2, "#"), "info1$h2 ==", info1$h2)
            expect(identical(info1$casedir, "/tmp/section1/case1"), "info1$casedir ==", info1$casedir)

            info2 <- casename_info("section1::case2", cases, "/tmp", create = FALSE)

            expect(identical(info2$casename, "section1::case2"), "info2$casename ==", info2$casename)
            expect(identical(info2$case, "case2"), "info2$case ==", info2$case)
            expect(identical(info2$case_slug, "case2"), "info2$case_slug ==", info2$case_slug)
            expect(identical(info2$section, "section1"), "info2$section ==", info2$section)
            expect(identical(info2$section_slug, "section1"), "info2$section_slug ==", info2$section_slug)
            expect(identical(info2$h1, "section1: case2"), "info2$h1 ==", info2$h1)
            expect(identical(info2$h2, "#"), "info2$h2 ==", info2$h2)
            expect(identical(info2$casedir, "/tmp/section1/case2"), "info2$casedir ==", info2$casedir)
        """

    @r_test
    def test_list_update(self):
        return """
            list1 <- list(a = 1, b = 2)
            list2 <- list(b = 3, c = 4)
            list3 <- list_update(list1, list2)

            expect(identical(list3$a, 1), "list3$a ==", list3$a)
            expect(identical(list3$b, 3), "list3$b ==", list3$b)
            expect(identical(list3$c, 4), "list3$c ==", list3$c)
        """

    @r_test
    def test_list_update_recursively(self):
        return """
            list1 <- list(a = 1, b = list(c = 2, d = list(e = 3, f = 4)))
            list2 <- list(b = list(c = 4, d = list(e = 5)))
            list3 <- list_update(list1, list2)

            expect(identical(list3$a, 1), "list3$a ==", list3$a)
            expect(identical(list3$b$c, 4), "list3$b$c ==", list3$b$c)
            expect(identical(list3$b$d$e, 5), "list3$b$d$e ==", list3$b$d$e)
            expect(identical(list3$b$d$f, 4), "list3$b$d$f ==", list3$b$d$f)
        """

    @r_test
    def test_list_update_recursively_with_depth(self):
        return """
            list1 <- list(a = 1, b = list(c = 2, d = list(e = 3, f = 4)))
            list2 <- list(b = list(c = 4, d = list(e = 5)))
            list3 <- list_update(list1, list2, depth = 1)

            expect(identical(list3$a, 1), "list3$a ==", list3$a)
            expect(identical(list3$b$c, 4), "list3$b$c ==", list3$b$c)
            expect(identical(list3$b$d$e, 5), "list3$b$d$e ==", list3$b$d$e)
            expect(is.null(list3$b$d$f), "list3$b$d$f ==", list3$b$d$f)
        """

    @r_test
    def test_list_update_with_null(self):
        return """
            list1 <- list(a = 1, b = 2)
            list2 <- list(b = NULL, c = 4)
            list3 <- list_update(list1, list2)

            expect(identical(list3$a, 1), "list3$a ==", list3$a)
            expect(is.null(list3$b), "list3$b ==", list3$b)
            expect(identical(list3$c, 4), "list3$c ==", list3$c)
            # "b" is not removed, but set to NULL
            expect(
                "b" %in% names(list3),
                "'b' not in ", paste(names(list3), collapse = ',')
            )
        """

    @r_test
    def test_list_setdefault(self):
        return """
            list1 <- list(a = 1)
            list1 <- list_setdefault(list1, a = 2)
            list1 <- list_setdefault(list1, b = 3)
            list1 <- list_setdefault(list1, c = NULL)

            expect(identical(list1$a, 1), "list1$a ==", list1$a)
            expect(identical(list1$b, 3), "list1$b ==", list1$b)
            expect(is.null(list1$c), "list1$c ==", list1$c)
            expect(
                "c" %in% names(list1),
                "'c' not in ", paste(names(list1), collapse = ',')
            )
        """

    @r_test
    def test_bquote(self):
        return """
            expect(isFALSE(.isBQuoted("a")))
            expect(isFALSE(.isBQuoted("`")))
            expect(isTRUE(.isBQuoted("``")))
            expect(isTRUE(.isBQuoted("`a`")))
            expect(identical(bQuote("a"), "`a`"))
            expect(identical(bQuote("`a`"), "`a`"))
        """

    @r_test
    def test_slugify(self):
        return """
            # create one including greek letters
            x <- "Title: α1 - β2"

            s1 <- slugify(x)
            expect(identical(s1, "Title-a1-b2"), "s1 ==", s1)

            s2 <- slugify(x, non_alphanum_replace = ".")
            expect(identical(s2, "Title.a1.b2"), "s2 ==", s2)

            s3 <- slugify(x, non_alphanum_replace = ".", tolower = TRUE)
            expect(identical(s3, "title.a1.b2"), "s3 ==", s3)

            s4 <- slugify(x, non_alphanum_replace = ".", collapse_replace = FALSE)
            expect(identical(s4, "Title..a1...b2"), "s4 ==", s4)
        """


if __name__ == "__main__":
    main(verbosity=2)
