from biopipen.core.proc import Proc, Pipeline


class MyPipeline(Pipeline):

    def build(self) -> None:

        class MyProc(Proc):
            input = "var:var"
            output = "outfile:file:{{in.var}}.out"
            script = f"""echo '{self.options}' > {{{{out.outfile}}}}"""

        self.starts.append(MyProc)
        self.procs.MyProc = Proc


pipe = MyPipeline(options={"a": "b"}).run(["var"])

outfile = pipe.workdir.joinpath("myproc", "0", "output", "var.out")

assert outfile.exists()
assert outfile.read_text().strip() == "{a: b}"
