from unittest import TestCase, main
from pathlib import Path
from io import StringIO
import tempfile
import sys
from biopipen.utils.misc import require_package, run_command, write_table, read_table


class TestUtilsMisc(TestCase):

    def test_require_package(self):
        # Test with an installed package without version check
        require_package("pip")

        # Test with an installed package with satisfying version
        require_package("pip", ">=0.0.1")

        # Test with an installed package with non-satisfying version
        with self.assertRaises(ImportError):
            require_package("pip", ">=9999.0.0")

        # Test with a non-installed package
        with self.assertRaises(ImportError):
            require_package("non_existent_package_12345")

    def test_require_package_with_python(self):
        # Test with current Python interpreter
        python_exe = sys.executable

        # Test with an installed package without version check
        require_package("pip", python=python_exe)

        # Test with an installed package with satisfying version
        require_package("pip", version=">=0.0.1", python=python_exe)

        # Test with an installed package with non-satisfying version
        with self.assertRaises(ImportError) as cm:
            require_package("pip", version=">=9999.0.0", python=python_exe)
        self.assertIn("does not satisfy the requirement", str(cm.exception))

        # Test with a non-installed package
        with self.assertRaises(ImportError) as cm:
            require_package(
                "non_existent_package_12345",
                python=python_exe
            )
        self.assertIn("is required but not installed", str(cm.exception))

        # Test with a non-existent Python interpreter
        with self.assertRaises(ImportError) as cm:
            require_package(
                "pip",
                python="/nonexistent/path/to/python"
            )
        self.assertIn("Python interpreter", str(cm.exception))
        self.assertIn("not found", str(cm.exception))

    def test_run_command_string(self):
        # Test running a simple string command
        result = run_command("echo hello", wait=True, print_command=False)
        self.assertEqual(result.returncode, 0)

    def test_run_command_list(self):
        # Test running a command as a list
        result = run_command(["echo", "hello"], wait=True, print_command=False)
        self.assertEqual(result.returncode, 0)

    def test_run_command_return_stdout(self):
        # Test returning stdout
        result = run_command("echo hello", stdout="RETURN", print_command=False)
        self.assertEqual(result.strip(), "hello")

        # Test with lowercase "return"
        result = run_command("echo world", stdout="return", print_command=False)
        self.assertEqual(result.strip(), "world")

    def test_run_command_stdout_pipe(self):
        # Test stdout=True (PIPE)
        result = run_command(["echo", "test"], stdout=True, wait=True, print_command=False)
        try:
            output = result.stdout.read().decode().strip()
            self.assertEqual(output, "test")
        finally:
            result.stdout.close()

    def test_run_command_stdout_file(self):
        # Test redirecting stdout to a file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            temp_file = Path(f.name)

        try:
            run_command("echo file_output", stdout=temp_file, wait=True, print_command=False)
            content = temp_file.read_text().strip()
            self.assertEqual(content, "file_output")
        finally:
            temp_file.unlink()

    def test_run_command_stdout_file_str(self):
        # Test redirecting stdout to a file (string path)
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            temp_file = f.name

        try:
            run_command("echo string_path", stdout=temp_file, wait=True, print_command=False)
            content = Path(temp_file).read_text().strip()
            self.assertEqual(content, "string_path")
        finally:
            Path(temp_file).unlink()

    def test_run_command_stderr_pipe(self):
        # Test stderr=True (PIPE)
        result = run_command(
            "python -c 'import sys; sys.stderr.write(\"error\")'",
            stderr=True,
            wait=True,
            print_command=False
        )
        try:
            error = result.stderr.read().decode().strip()
            self.assertEqual(error, "error")
        finally:
            result.stderr.close()

    def test_run_command_stderr_stdout(self):
        # Test redirecting stderr to stdout
        result = run_command(
            "python -c 'import sys; sys.stderr.write(\"error\")'",
            stdout=True,
            stderr="STDOUT",
            wait=True,
            print_command=False
        )
        try:
            output = result.stdout.read().decode().strip()
            self.assertEqual(output, "error")
        finally:
            result.stdout.close()

    def test_run_command_env(self):
        # Test passing environment variables
        result = run_command(
            "echo $TEST_VAR",
            env={"TEST_VAR": "test_value"},
            stdout="RETURN",
            print_command=False
        )
        self.assertEqual(result.strip(), "test_value")

    def test_run_command_failed(self):
        # Test that a failed command raises RuntimeError
        with self.assertRaises(RuntimeError) as cm:
            run_command("exit 1", wait=True, print_command=False)
        self.assertIn("Failed to run command: rc=1", str(cm.exception))

    def test_run_command_failed_return_stdout(self):
        # Test that a failed command with return stdout raises RuntimeError
        with self.assertRaises(RuntimeError) as cm:
            run_command("exit 1", stdout="RETURN", print_command=False)
        self.assertIn("Failed to run command: rc=1", str(cm.exception))

    def test_run_command_no_wait(self):
        # Test running a command without waiting
        result = run_command("sleep 0.1", wait=False, print_command=False)
        # Process should still be running or just finished
        self.assertIsNotNone(result)
        result.wait()  # Clean up
        self.assertEqual(result.returncode, 0)

    def test_run_command_print_command(self):
        # Test that the command is printed
        captured_output = []

        def capture_print(*args):
            captured_output.append(str(args[0]) if args else "")

        run_command(
            "echo test",
            wait=True,
            print_command=True,
            print_command_handler=capture_print
        )

        self.assertIn("RUNNING COMMAND:", captured_output)
        self.assertTrue(any("echo test" in line for line in captured_output))

    def test_run_command_no_print(self):
        # Test disabling command printing
        captured_output = []

        def capture_print(*args):
            captured_output.append(str(args[0]) if args else "")

        run_command(
            "echo test",
            wait=True,
            print_command=False,
            print_command_handler=capture_print
        )

        self.assertEqual(len(captured_output), 0)

    def test_run_command_fg_error(self):
        # Test that fg=True with stdout/stderr redirection raises ValueError
        with self.assertRaises(ValueError) as cm:
            run_command("echo test", fg=True, stdout=True, print_command=False)
        self.assertIn("Cannot redirect stdout or stderr when running in foreground", str(cm.exception))

        with self.assertRaises(ValueError) as cm:
            run_command("echo test", fg=True, stderr=True, print_command=False)
        self.assertIn("Cannot redirect stdout or stderr when running in foreground", str(cm.exception))

    def test_run_command_invalid_command(self):
        # Test that an invalid command raises RuntimeError
        with self.assertRaises(RuntimeError) as cm:
            run_command("nonexistent_command_12345", wait=True, print_command=False)
        self.assertIn("Failed to run command", str(cm.exception))

    def test_run_command_list_with_numbers(self):
        # Test that list elements are converted to strings
        result = run_command(["echo", 123, 456], stdout="RETURN", print_command=False)
        self.assertEqual(result.strip(), "123 456")

    def setUp(self):
        import pandas as pd
        self.df = pd.DataFrame({
            'A': pd.Categorical(['a', 'b', 'a'], categories=['a', 'b', 'c']),
            'B': [1, 2, 3]
        })

    def test_write_table_str(self):
        # No factor levels, no path -> no leading newline
        result = write_table(self.df[['B']], index=False)
        self.assertEqual(result, "B\n1\n2\n3\n")
        # With factor levels
        result = write_table(self.df, index=False)
        self.assertEqual(
            result,
            "# factor-levels: A=a|b|c\nA,B\na,1\nb,2\na,3\n"
        )

    def test_write_table_path(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            path = f.name
        try:
            write_table(self.df, path, index=False)
            self.assertEqual(
                Path(path).read_text(),
                "# factor-levels: A=a|b|c\nA,B\na,1\nb,2\na,3\n"
            )
        finally:
            Path(path).unlink()

    def test_write_table_buffer(self):
        buf = StringIO()
        write_table(self.df, buf, index=False)
        buf.seek(0)
        self.assertEqual(
            buf.getvalue(),
            "# factor-levels: A=a|b|c\nA,B\na,1\nb,2\na,3\n"
        )

    def test_read_table_path(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
            path = f.name
        try:
            write_table(self.df, path)
            df = read_table(path)
            self.assertEqual(
                df['A'].cat.categories.tolist(), ['a', 'b', 'c']
            )
            self.assertEqual(df['B'].tolist(), [1, 2, 3])
        finally:
            Path(path).unlink()

    def test_read_table_buffer(self):
        buf = StringIO()
        write_table(self.df, buf)
        buf.seek(0)
        df = read_table(buf)
        self.assertEqual(df['A'].cat.categories.tolist(), ['a', 'b', 'c'])
        self.assertEqual(df['B'].tolist(), [1, 2, 3])

    def test_read_table_no_factor_levels(self):
        # Plain csv with no factor-level lines
        buf = StringIO("B\n1\n2\n3\n")
        df = read_table(buf)
        self.assertEqual(df['B'].tolist(), [1, 2, 3])


if __name__ == "__main__":
    main(verbosity=2)
