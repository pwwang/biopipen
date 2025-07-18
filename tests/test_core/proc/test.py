from unittest import TestCase, main
from biopipen.core.proc import Proc


class TestProc(TestCase):

    def test_proc(self):

        class P1(Proc):
            pass

        # plugin_opts should be inherited
        self.assertIn("poplog_pattern", P1.plugin_opts)


if __name__ == "__main__":
    main(verbosity=2)
