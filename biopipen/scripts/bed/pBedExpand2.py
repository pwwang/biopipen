from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
base = {{args.base | int}}
binsize = {{args.binsize | repr}}
datacols = {{args.datacols | ?isinstance: str | =.split: ","
             			   | @map: _, int | $repr}}
name = {{args.name | repr}}
origcols = {{args.origcols | quote}}
aggrs = {{args.aggrs | ?isinstance: str | =:[_]
                     | $repr }}
for i, aggr in enumerate(aggrs):
    if aggr == "raw":
        aggr = lambda data: ",".join(str(dat) for dat in data)
    elif aggr == "mean":
        aggr = lambda data: sum(float(dat) for dat in data) / len(data)
    elif aggr == "sum":
        aggr = lambda data: sum(float(dat) for dat in data)
    elif aggr == "min":
        aggr = min
    elif aggr == "max":
        aggr = max
    elif aggr == "count":
        aggr = len
    elif aggr.startswith("lambda"):
        aggr = eval(aggr)
    aggrs[i] = aggr

if len(aggrs) == 1:
    aggrs = aggrs * len(datacols)

# required arguments
if not binsize:
    raise ValueError("Require argument `args.binsize`")
if not datacols:
    raise ValueError("Require argument `args.datacols`")

reader = TsvReader(infile, cnames=False)
writer = TsvWriter(outfile)
for ridx, row in enumerate(reader):
    chrom = row[0]
    start = int(row[1]) if base == 0 else int(row[1]) - 1
    end = int(row[2])
    origname = row[3] if len(row) > 3 else f"src{ridx+1}"
    datas = [[dat.strip() for dat in row[dcol-1].split(",")]
             for dcol in datacols]

    # length checking
    if (end - start) % binsize != 0:
        raise ValueError(f"Length of region {row} cannot be divided by binsize")
    nbins = (end - start) // binsize
    for i, data in enumerate(datas):
        if len(data) % nbins != 0:
            raise ValueError(f"Number of data {len(data)} "
                             f"at column {datacols[i]} cannot be divided by "
                             f"number of bins for region {row[:3]}")
        chunk_size = len(data) // nbins
        data = [data[n*chunk_size: (n+1)*chunk_size] for n in range(nbins)]
        # datas: mark_index => bin_index => [data]
        datas[i] = data

    for binidx in range(nbins):
        outs = [chrom, start + binidx * binsize, start + (binidx + 1) * binsize]
        if name == 'src':
            outs.append(origname)
        elif name == 'srcbin':
            outs.append(f"{origname}_bin{binidx+1}")
        if origcols == 'keep':
            outs.extend(row[4:])
        # add aggregation value of each mark
        for i, data in enumerate(datas):
            outs.append(aggrs[i](data[binidx]))
        writer.write(outs)

reader.close()
writer.close()
