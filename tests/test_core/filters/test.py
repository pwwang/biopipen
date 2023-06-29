from biopipen.core.filters import dict_to_cli_args


def assert_equal(expected, result):
    assert expected == result, f"Expected {expected}, but got {result}"


def test_dict_to_cli_args():
    # Test empty dictionary
    assert_equal(dict_to_cli_args({}), [])

    # Test single key-value pair
    assert_equal(dict_to_cli_args({"a": 1}), ["-a", "1"])

    # Test multiple key-value pairs
    assert_equal(dict_to_cli_args({"a": 1, "b": 2}), ["-a", "1", "-b", "2"])

    # Test boolean values
    assert_equal(dict_to_cli_args({"a": True, "b": False}), ["-a"])

    # Test list values
    assert_equal(
        dict_to_cli_args({"a": [1, 2], "b": [3]}),
        ["-a", "1", "-a", "2", "-b", "3"],
    )

    # Test list values with duplicate key
    assert_equal(
        dict_to_cli_args({"a": [1, 2], "b": [3]}, dup_key=True),
        ["-a", "1", "-a", "2", "-b", "3"],
    )

    # Test list values with separator
    assert_equal(
        dict_to_cli_args({"a": [1, 2]}, sep="=", dup_key=True),
        ["-a=1", "-a=2"],
    )

    # Test nested dictionary
    assert_equal(
        dict_to_cli_args({"a": {"b": 1, "c": 2}}),
        ['-a', "{'b': 1, 'c': 2}"],
    )

    # Test dictionary with non-string keys
    assert_equal(
        dict_to_cli_args({1: "a", True: "b"}),
        ["-1", "b"],
    )

    # Test dictionary with non-string values
    assert_equal(
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
    assert_equal(
        dict_to_cli_args({"a": 1, "b": 2}, prefix="--"),
        ["--a", "1", "--b", "2"]
    )

    # Test start
    assert_equal(
        dict_to_cli_args({"a": 1, "_": 3, "": 2}),
        ["2", "-a", "1", "3"]
    )

    # Test dashify
    assert_equal(
        dict_to_cli_args({"a": 1, "_": 3, "a_b": 4, "": 2}, dashify=True),
        ["2", "-a", "1", "--a-b", "4", "3"]
    )


if __name__ == "__main__":
    test_dict_to_cli_args()
