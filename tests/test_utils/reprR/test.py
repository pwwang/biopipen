from pathlib import Path
from unittest import TestCase, main
from biopipen.core.testing import r_test


class TestUtilsReprR(TestCase):
    SOURCE_FILE = [
        Path(__file__).parent.parent.parent.parent.joinpath(
            "biopipen", "utils", "misc.R"
        ),
        Path(__file__).parent.parent.parent.parent.joinpath(
            "biopipen", "utils", "repr.R"
        ),
    ]

    @r_test
    def test_reprR_default_dot_repr_protocol(self):
        return """
            K <- setRefClass("K",
                fields = list(x = "numeric"),
                methods = list(.repr = function(newline) { x + 1 }))

            k <- K$new(x = 1)
            expect(identical(repr(k), 2))
        """

    @r_test
    def test_reprR_default_fallback(self):
        return """
            k <- 1:3
            class(k) <- c("A", "B")
            expect(identical(repr(k), "<A/B: k>"))
        """

    @r_test
    def test_reprR_numeric(self):
        return """
            k <- 1
            expect(identical(repr(k), "1"))
        """

    @r_test
    def test_reprR_character(self):
        return """
            k <- "a"
            expect(identical(repr(k), '"a"'))
        """

    @r_test
    def test_reprR_logical(self):
        return """
            k <- TRUE
            expect(identical(repr(k), "TRUE"))
        """

    @r_test
    def test_reprR_factor(self):
        return """
            k <- factor("a")
            expect(identical(repr(k), 'factor("a", levels = "a")'))
        """

    @r_test
    def test_reprR_list(self):
        return """
            k <- list(a = 1, b = 2)
            expect(
                identical(repr(k, newline = TRUE), 'list(\n  `a` = 1,\n  `b` = 2\n)'),
                repr(k, newline = TRUE)
            )
        """

    @r_test
    def test_reprR_data_frame(self):
        return """
            k <- data.frame(a = 1, b = 2)
            expect(
                identical(repr(k, newline = TRUE), 'data.frame(\n  `a` = 1,\n  `b` = 2\n)'),
                repr(k, newline = TRUE)
            )
        """

    @r_test
    def test_reprR_NULL(self):
        return """
            k <- NULL
            expect(identical(repr(k), "NULL"))
        """

    @r_test
    def test_reprR_environment(self):
        return """
            k <- rlang::env(a = 1)
            expect(identical(repr(k), "rlang::env( `a` = 1 )"), repr(k))
        """


if __name__ == "__main__":
    main(verbosity=2)
