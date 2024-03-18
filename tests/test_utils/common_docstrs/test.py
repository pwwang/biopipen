from unittest import TestCase, main
from biopipen.utils.common_docstrs import indent_docstr, format_placeholder


class TestCommonDocstrs(TestCase):

    def test_indent_docstr(self):
        self.assertEqual(indent_docstr("1\n2\n3", "    "), "1\n    2\n    3")
        self.assertEqual(indent_docstr("  \n1\n2\n3\n  ", " "), "1\n 2\n 3")

    def test_format_placeholder(self):
        @format_placeholder(test="test")
        class Test:
            """%(test)s"""
        self.assertEqual(Test.__doc__, "test")

        @format_placeholder(test="1\n2")
        class Test:
            """%(test)s"""
        self.assertEqual(Test.__doc__, "1\n2")


if __name__ == "__main__":
    main(verbosity=2)
