from biopipen.utils.common_docstrs import indent_docstr, format_placeholder

TEST_NO = 1


def run_indent_docstr(in_, prefix, out):
    global TEST_NO
    print(f">>> TESTING indent_docstr #{TEST_NO}")
    TEST_NO += 1
    assert indent_docstr(in_, prefix) == out
    print(">>> PASSED")
    print(">>> ")


def run_format_placeholder(in_, out):
    global TEST_NO
    print(f">>> TESTING format_placeholder #{TEST_NO}")
    TEST_NO += 1

    @format_placeholder(**in_)
    class Test:
        """%(test)s"""

    assert Test.__doc__ == out
    print(">>> PASSED")
    print(">>> ")


if __name__ == "__main__":
    run_indent_docstr("1\n2\n3", "    ", "1\n    2\n    3")
    run_indent_docstr("  \n1\n2\n3\n  ", " ", "1\n 2\n 3")
    run_format_placeholder({"test": "test"}, "test")
    run_format_placeholder({"test": "1\n2"}, "1\n2")
