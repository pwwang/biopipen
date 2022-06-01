
infile = {{in.infile | quote}}
outprefix = {{out.outfile | prefix | replace: ".fancyvj.wt", "" | quote}}
vdjtools = {{ envs.vdjtools | quote }}
vdjtools_patch = {{ envs.vdjtools_patch | quote }}
joboutdir = {{job.outdir | quote}}

command = sprintf(
    "cd %s && bash %s %s PlotFancyVJUsage --plot-type png %s %s",
    joboutdir, vdjtools_patch, vdjtools, infile, outprefix
)

print("Running command:")
print(command)

system(command)
