from unittest import TestCase, main
from pathlib import Path
from biopipen.core.filters import dict_to_cli_args, r


class TestFilters(TestCase):

    def test_dict_to_cli_args(self):
        # Test empty dictionary
        self.assertEqual(dict_to_cli_args({}), [])

        # Test single key-value pair
        self.assertEqual(dict_to_cli_args({"a": 1}), ["-a", "1"])

        # Test multiple key-value pairs
        self.assertEqual(dict_to_cli_args({"a": 1, "b": 2}), ["-a", "1", "-b", "2"])

        # Test boolean values
        self.assertEqual(dict_to_cli_args({"a": True, "b": False}), ["-a"])

        # Test list values
        self.assertEqual(
            dict_to_cli_args({"a": [1, 2], "b": [3]}),
            ["-a", "1", "-a", "2", "-b", "3"],
        )

        # Test list values with duplicate key
        self.assertEqual(
            dict_to_cli_args({"a": [1, 2], "b": [3]}, dup_key=True),
            ["-a", "1", "-a", "2", "-b", "3"],
        )

        # Test list values with separator
        self.assertEqual(
            dict_to_cli_args({"a": [1, 2]}, sep="=", dup_key=True),
            ["-a=1", "-a=2"],
        )

        # Test nested dictionary
        self.assertEqual(
            dict_to_cli_args({"a": {"b": 1, "c": 2}}),
            ['-a', "{'b': 1, 'c': 2}"],
        )

        # Test dictionary with non-string keys
        self.assertEqual(
            dict_to_cli_args({1: "a", True: "b"}),  # noqa: F601
            ["-1", "b"],
        )

        # Test dictionary with non-string values
        self.assertEqual(
            dict_to_cli_args({"a": 1, "b": True}),
            ["-a", "1", "-b"],
        )

        # Test dictionary with non-hashable keys
        ex = None
        try:
            dict_to_cli_args({[1, 2]: "a"})
        except TypeError as e:
            ex = e
        assert isinstance(ex, TypeError), "Expected TypeError"

        # Test sep='=' and dup_key=True
        ex = None
        try:
            dict_to_cli_args({"a": [1, 2]}, sep="=", dup_key=False)
        except ValueError as e:
            ex = e
        assert isinstance(ex, ValueError), "Expected ValueError"

        # Test prefix
        self.assertEqual(
            dict_to_cli_args({"a": 1, "b": 2}, prefix="--"),
            ["--a", "1", "--b", "2"]
        )

        # Test start
        self.assertEqual(
            dict_to_cli_args({"a": 1, "_": 3, "": 2}),
            ["2", "-a", "1", "3"]
        )

        # Test dashify
        self.assertEqual(
            dict_to_cli_args({"a": 1, "_": 3, "a_b": 4, "": 2}, dashify=True),
            ["2", "-a", "1", "--a-b", "4", "3"]
        )

    def test_r(self):
        self.assertEqual(r(True), "TRUE")
        self.assertEqual(r(False), "FALSE")
        self.assertEqual(r(None), "NULL")
        self.assertEqual(r("+INF"), "Inf")
        self.assertEqual(r("-INF"), "-Inf")
        self.assertEqual(r("TRUE"), "TRUE")
        self.assertEqual(r("FALSE"), "FALSE")
        self.assertEqual(r("na"), "NA")
        self.assertEqual(r("null"), "NULL")
        self.assertEqual(r("r:c(1, 2, 3)"), "c(1, 2, 3)")
        self.assertEqual(r("R:abc"), "abc")
        self.assertEqual(r(Path("/x/y/z")), "'/x/y/z'")
        self.assertEqual(r([1, 2, 3]), "c(1, 2, 3)")
        self.assertEqual(r(({"a": 1}, 2, 3)), "list(list(`a`=1), 2, 3)")
        self.assertEqual(r(({"a": 1}, 2, 3)), "list(list(`a`=1), 2, 3)")
        self.assertEqual(
            r({2: 0, 0: 1, 1: 2}, ignoreintkey=False, sortkeys=True),
            "list(`0`=1, `1`=2, `2`=0)",
        )
        self.assertEqual(
            r({2: 0, 0: 1, 1: 2}, ignoreintkey=True, sortkeys=True),
            "list(1, 2, 0)",
        )
        self.assertEqual(
            r({2: 0, 0: 1, 1: 2}, ignoreintkey=False, sortkeys=False),
            "list(`2`=0, `0`=1, `1`=2)",
        )
        self.assertEqual(
            r({2: 0, 0: 1, 1: 2}, ignoreintkey=True, sortkeys=False),
            "list(0, 1, 2)",
        )
        self.assertEqual(
            r({"a-b": {"c-d": 1}}, todot="-", skip=0),
            "list(`a.b`=list(`c.d`=1))",
        )
        self.assertEqual(
            r({"a-b": {"c-d": 1}}, todot="-", skip=1),
            "list(`a-b`=list(`c.d`=1))",
        )


if __name__ == "__main__":
    main(verbosity=2)
