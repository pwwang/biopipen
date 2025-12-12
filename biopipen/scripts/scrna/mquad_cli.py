# CLI for mitoMut

import os
import sys
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy.io import mmread
from scipy.io import mmwrite

# from scipy import sparse
from scipy.sparse import csc_matrix, csr_matrix
from optparse import OptionParser, OptionGroup

from mquad.version import __version__
from vireoSNP.utils.io_utils import read_cellSNP, read_sparse_GeneINFO
from vireoSNP.utils.vcf_utils import load_VCF

from mquad.mquad import Mquad
from mquad.mquad_utils import findKnee
from mquad.mquad_batch_mixbin import MquadSparseMixBin, fit_batch

START_TIME = time.time()


def fit_batch2(*args, **kwargs):
    try:
        return fit_batch(*args, **kwargs)
    except Exception as e:
        print(f"Warning: Batch failed with error: {e}")
        return None


def fit_deltaBIC(
    self,
    out_dir,
    minDP=10,
    minAD=1,
    export_csv=True,
    nproc=30,
    batch_size=128,
):
    # here we fit and choose model based on deltaBIC
    n_variants = self.dp.shape[0]

    # adjust batch size to avoid unused processes
    adj_bs = min(n_variants // nproc, batch_size)
    print(
        "CPUs used: {}, batch size: {} {}".format(
            nproc,
            adj_bs,
            "" if adj_bs == batch_size else "(adjusted to avoid idle processes)",
        )
    )
    print("Fitting in sparse mode...")
    t0 = time.time()

    print(
        "Initializing fit(mode: deltaBIC) on " + str(self.ad.shape[0]) + " variants..."
    )

    dp_row, dp_col = np.nonzero(self.dp >= minDP)
    ad_row, _ = np.nonzero(self.ad >= minAD)

    # only variant with at leat one valid dp and ad records are included
    valid_rows = np.intersect1d(dp_row, ad_row)
    # filter invalid variants
    x, y = zip(*[(r, c) for r, c in zip(dp_row, dp_col) if r in valid_rows])

    # split batch
    x, y = np.array(x), np.array(y)
    valid_row_sizes = Counter(x)
    valid_row_sizes = np.array([valid_row_sizes[r_idx] for r_idx in valid_rows])

    assert np.sum(valid_row_sizes) == x.shape[0]

    with mp.Pool(processes=nproc) as pool:
        results = pool.starmap_async(
            fit_batch2, self._batchify(batch_size, x, y, valid_row_sizes)
        ).get()

    # Filter out None results and track successful indices
    successful_results = []
    successful_variant_indices = []
    current_batch_start = 0

    for batch_idx, result in enumerate(results):
        if result is not None:
            batch_variants = len(result["deltaBIC"])
            successful_results.append(pd.DataFrame(result))
            # Track which variants from valid_rows succeeded
            successful_variant_indices.extend(
                range(current_batch_start, current_batch_start + batch_variants)
            )
            current_batch_start += batch_variants
        else:
            # Skip failed batch
            print(f"Warning: Batch {batch_idx} failed and will be skipped")

    if not successful_results:
        raise RuntimeError("All batches failed during fitting!")

    self.df = pd.concat(successful_results, axis=0, ignore_index=True)

    # Store successful variant indices for later use
    self.successful_variant_indices = np.array(successful_variant_indices)

    t1 = time.time()
    print(
        "deltaBIC was calculated for "
        + str(self.ad.shape[0])
        + " variants and took:%.2f minutes" % ((t1 - t0) / 60)
    )

    #         self.df = pd.DataFrame()
    #         for col, res in results.items():
    #             self.df[col] = res.tolist()

    if self.variants is not None:
        # self._addVariantNames(valid_rows)
        self._addVariantNames(valid_rows[self.successful_variant_indices])

    # sort df but keep index
    self.sorted_df = self.df.sort_values(by=["deltaBIC"], ascending=False)

    if export_csv is True:
        if self._check_outdir_exist(out_dir) is True:
            self.sorted_df.to_csv(out_dir + "/BIC_params.csv", index=False)
        else:
            self.sorted_df.to_csv("BIC_params.csv", index=False)

    self.df.to_csv(out_dir + "/debug_unsorted_BIC_params.csv", index=False)
    # return df of all metrics
    return self.df


def selectInformativeVariants(
    self,
    min_cells=2,
    export_heatmap=True,
    export_mtx=True,
    out_dir=None,
    existing_df=None,
    tenx_cutoff=None,
):
    # takes self.df, return best_ad and best_dp as array

    if self.df is None:
        print("Fitted model not found! Have you run fit_deltaBIC/fit_logLik yet?")
    else:
        if out_dir is not None:
            if os.path.exists(out_dir) is not True:
                try:
                    os.mkdir(out_dir)
                except Exception:
                    print("Can't make directory, do you have permission?")
        else:
            print("Out directory already exists, overwriting content inside...")

        if tenx_cutoff is None:
            x, y, knee, cutoff = findKnee(self.df.deltaBIC)

            plt.plot(x, y)
            plt.axvline(
                x=knee,  # type: ignore
                color="black",
                linestyle="--",
                label="cutoff",
            )
            plt.legend()
            plt.ylabel("\u0394BIC")
            plt.xlabel("Cumulative probability")
            plt.savefig(out_dir + "/deltaBIC_cdf.pdf")

            # make a PASS/FAIL column in self.df for easier subsetting
            print("deltaBIC cutoff = ", cutoff)
            # self.sorted_df['VALID'] = self.validateSNP(self.sorted_df.variant_name)
            self.sorted_df["PASS_KP"] = self.sorted_df.deltaBIC.apply(
                lambda x: True if x >= cutoff else False
            )
            self.sorted_df["PASS_MINCELLS"] = self.sorted_df.num_cells_minor_cpt.apply(
                lambda x: True if x >= min_cells else False
            )

            self.final_df = self.sorted_df[
                (self.sorted_df.PASS_KP == True)
                & (self.sorted_df.PASS_MINCELLS == True)
            ]
            # print(self.final_df.head())

            # will deprecate in later versions
            # self.final_df = self.sorted_df[0:int(len(y) * (1 - knee))]
            # self.final_df = self.final_df[
            #   self.sorted_df.num_cells_minor_cpt >= min_cells]

            print(
                "Number of variants passing threshold: "
                + str(len(self.final_df["variant_name"]))
            )

        else:
            print("[MQuad] Tenx mode used with cutoff = " + str(tenx_cutoff))
            self.final_df = self.sorted_df[
                self.sorted_df.deltaBIC >= float(tenx_cutoff)
            ]
            self.final_df = self.final_df[
                self.sorted_df.num_cells_minor_cpt >= min_cells
            ]
            print(
                "Number of variants passing threshold: "
                + str(len(self.final_df["variant_name"]))
            )

        if len(self.final_df["variant_name"]) != 0:
            passed_variants = self.final_df["variant_name"]
            idx = [self.variants.index(i) for i in passed_variants]

            best_ad = self.ad[idx]
            best_dp = self.dp[idx]
        else:
            print(
                "No informative variants detected! If you are using 10x data, "
                "try setting --minDP to a smaller number."
            )
            best_ad = self.ad[[]]
            best_dp = self.dp[[]]
            passed_variants = []

    self.sorted_df.to_csv(out_dir + "/BIC_params.csv", index=False)
    # fname = by + '_' + str(threshold) + '_'

    if self.variants is not None:
        # best_vars = np.array(self.variants)[idx]
        renamed_vars = []
        for var in passed_variants:
            renamed_vars.append(
                (var.split("_")[1] + var.split("_")[2] + ">" + var.split("_")[3])
            )

        with open(out_dir + "/passed_variant_names.txt", "w+") as var_file:
            var_file.write("\n".join(str(var) for var in renamed_vars))

    if export_heatmap is True:
        af = best_ad / best_dp
        # print(af.shape)
        # af = af.fillna(0)
        fig, ax = plt.subplots(figsize=(15, 10))
        plt.title("Allele frequency of top variants")
        plt.style.use("seaborn-v0_8-dark")
        if self.variants and renamed_vars:
            sns.heatmap(af, cmap="Greens", yticklabels=renamed_vars)
        else:
            # sns.heatmap(af, cmap="Greens")
            # generate a plot with text "No variants passed the threshold"
            plt.text(
                0.5,
                0.5,
                "No variants passed the threshold",
                ha="center",
                va="center",
                fontsize=20,
            )
            plt.axis("off")

        plt.savefig(out_dir + "/top variants heatmap.pdf")

    # export ad dp mtx out for vireo
    if export_mtx is True:
        mmwrite(out_dir + "/passed_ad.mtx", csr_matrix(best_ad))
        mmwrite(out_dir + "/passed_dp.mtx", csr_matrix(best_dp))

    return best_ad, best_dp


MquadSparseMixBin.fit_deltaBIC = fit_deltaBIC  # type: ignore
MquadSparseMixBin.selectInformativeVariants = selectInformativeVariants  # type: ignore


def main():
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option(
        "--cellData",
        "-c",
        dest="cell_data",
        default=None,
        help=("The cellSNP folder with AD and DP sparse matrices."),
    )
    parser.add_option(
        "--outDir",
        "-o",
        dest="out_dir",
        default=None,
        help=("Dirtectory for output files [default: $cellData/mitoMut]"),
    )

    group0 = OptionGroup(parser, "Alternative input formats")
    group0.add_option(
        "--mtxData",
        "-m",
        dest="mtx_data",
        default=None,
        help=("The two mtx files for AD and DP matrices, comma separated"),
    )
    parser.add_option(
        "--vcfData",
        dest="vcf_data",
        default=None,
        help=("The cell genotype file in VCF format"),
    )
    group0.add_option(
        "--BICparams",
        "--b",
        dest="BIC_params",
        default=None,
        help=("Existing unsorted_debug_BIC_params.csv"),
    )
    group0.add_option(
        "--tenx",
        "--t",
        dest="cutoff",
        default=None,
        help=("User-defined deltaBIC cutoff mainly for low-depth data"),
    )

    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option(
        "--randSeed",
        type="int",
        dest="rand_seed",
        default=None,
        help="Seed for random initialization [default: %default]",
    )
    group1.add_option(
        "--nproc",
        "-p",
        type="int",
        dest="nproc",
        default=1,
        help="Number of subprocesses [default: %default]",
    )
    group1.add_option(
        "--minDP",
        type="int",
        dest="minDP",
        default=10,
        help="Minimum DP to include for modelling [default: 10]",
    )
    group1.add_option(
        "--minCell",
        type="int",
        dest="minCell",
        default=2,
        help=("Minimum no. of cells in minor component [default: 2]"),
    )
    group1.add_option(
        "--batchFit",
        type="int",
        dest="batch_fit",
        default=1,
        help=("1 if fit MixBin model using batch mode, 0 else [default: 1]"),
    )
    group1.add_option(
        "--batchSize",
        type="int",
        dest="batch_size",
        default=128,
        help=(
            "Number of variants in one batch, cooperate with "
            "--nproc for speeding up [default: 128]"
        ),
    )
    group1.add_option(
        "--beta",
        type="int",
        dest="beta_mode",
        default=0,
        help=("Use betabinomial model if True [default: 0]"),
    )

    parser.add_option_group(group0)
    parser.add_option_group(group1)
    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        print("Welcome to MQuad v%s!\n" % (__version__))
        print("use -h or --help for help on argument.")
        sys.exit(1)

    # out directory
    if options.out_dir is None:
        print("Warning: no outDir provided, we use $cellFilePath/mquad")
        out_dir = os.path.dirname(os.path.abspath(options.cell_file)) + "/mquad"
    elif os.path.dirname(options.out_dir) == "":
        out_dir = "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # input data (VCF.gz or a folder with sparse matrices)
    if (
        options.cell_data is None
        and options.mtx_data is None
        and options.vcf_data is None
    ):
        print(
            "Error: need cell data in vcf file, or cellSNP output folder, or "
            "matrix files: AD.mtx,DP.mtx"
        )
        sys.exit(1)
    elif options.mtx_data is not None:
        print("[MQuad] Loading matrix files ...")
        matrix_files = options.mtx_data.split(",")
        if len(matrix_files) != 2:
            print("Error: mtxData requires 2 comma separated files")
            sys.exit(1)
        cell_dat = {}
        cell_dat["AD"] = csc_matrix(mmread(matrix_files[0]))
        cell_dat["DP"] = csc_matrix(mmread(matrix_files[1]))
        cell_dat["variants"] = [
            "SNP%d" % (x + 1) for x in range(cell_dat["AD"].shape[0])  # type: ignore
        ]

    elif options.vcf_data is not None:
        print("[MQuad] Loading cell VCF file ...")
        cell_vcf = load_VCF(options.vcf_data, biallelic_only=True)
        cell_dat = read_sparse_GeneINFO(cell_vcf["GenoINFO"], keys=["AD", "DP"])
        for _key in ["samples", "variants", "FixedINFO", "contigs", "comments"]:
            cell_dat[_key] = cell_vcf[_key]

    else:
        # os.path.isdir(os.path.abspath(options.cell_data)):
        print("[MQuad] Loading cell folder ...")
        cell_dat = read_cellSNP(options.cell_data)

    # More options
    nproc = options.nproc
    minDP = options.minDP
    batch_size = options.batch_size
    cutoff = options.cutoff
    minCell = options.minCell
    if options.beta_mode == 1:
        beta = True
    elif options.beta_mode == 0:
        beta = False

    # Main functions
    if options.BIC_params is not None:
        mdphd = Mquad(
            AD=cell_dat["AD"], DP=cell_dat["DP"], variant_names=cell_dat["variants"]
        )
        print("[MQuad] Using existing BIC params to filter variants only...")
        best_ad, best_dp = mdphd.selectInformativeVariants(
            min_cells=minCell,
            out_dir=out_dir,
            existing_df=options.BIC_params,
            tenx_cutoff=cutoff,
        )
    else:
        if options.batch_fit == 0:
            mdphd = Mquad(
                AD=cell_dat["AD"], DP=cell_dat["DP"], variant_names=cell_dat["variants"]
            )
            df = mdphd.fit_deltaBIC(
                out_dir=out_dir, nproc=nproc, minDP=minDP, beta_mode=beta  # type: ignore
            )
            best_ad, best_dp = mdphd.selectInformativeVariants(
                min_cells=minCell, out_dir=out_dir, tenx_cutoff=cutoff
            )
        else:
            # use sparse mode for faster performance
            # default to sparse mode in v0.1.6
            mdphd = MquadSparseMixBin(
                AD=cell_dat["AD"], DP=cell_dat["DP"], variant_names=cell_dat["variants"]
            )
            df = mdphd.fit_deltaBIC(  # noqa
                out_dir=out_dir, minDP=minDP, nproc=nproc, batch_size=batch_size
            )
            best_ad, best_dp = mdphd.selectInformativeVariants(
                min_cells=minCell, out_dir=out_dir, tenx_cutoff=cutoff
            )

    run_time = time.time() - START_TIME
    print("[MQuad] All done: %d min %.1f sec" % (int(run_time / 60), run_time % 60))
    print()


if __name__ == "__main__":
    main()
