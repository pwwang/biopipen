from biopipen.core.proc import Proc
from pipen_cli_run import Pipeline


class MyPipeline(Pipeline):

    def build(self) -> None:

        class MyProc(Proc):
            input = "var:var"
            output = "outfile:file:{{in.var}}.out"
            script = f"""echo '{self.options}' > {{{{out.outfile}}}}"""

        self.starts.append(MyProc)
        self.procs.MyProc = MyProc


pipe = MyPipeline(options={"a": "b"}).run(["var"], plugins=["no:report"])

outfile = pipe.workdir.joinpath("myproc", "0", "output", "var.out")

assert outfile.exists()
assert outfile.read_text().strip() == "{a: b}"
