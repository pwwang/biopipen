
infile = {{in.infile | r}}
outprefix = {{out.outfile | prefix | replace: ".fancyvj.wt", "" | r}}
vdjtools = {{ envs.vdjtools | r }}
vdjtools_patch = {{ envs.vdjtools_patch | r }}
joboutdir = {{job.outdir | r}}

command = sprintf(
    "cd %s && bash %s %s PlotFancyVJUsage --plot-type png %s %s",
    joboutdir, vdjtools_patch, vdjtools, infile, outprefix
)

print("Running command:")
print(command)

system(command)
