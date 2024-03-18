from unittest import TestCase, main

from pathlib import Path
from tempfile import gettempdir
from biopipen.core.testing import r_test

EXDATA_FILE = Path(gettempdir()).joinpath("biopipen-test_utils-single_cell-exdata.txt")
FLDATA_FILE = Path(gettempdir()).joinpath("biopipen-test_utils-single_cell-filtered.txt")
REDATA_FILE = Path(gettempdir()).joinpath("biopipen-test_utils-single_cell-recovered.txt")


class TestSingleCell(TestCase):

    SOURCE_FILE = Path(__file__).parent.parent.parent.parent.joinpath(
        "biopipen", "utils", "single_cell.R"
    )

    @r_test
    def test_expand_immdata(self):
        return f"""
            samples <- names(immdata$data)
            all_cells <- list()
            df <- lapply(samples, function(x) {{
                all_cells[[x]] <<- sum(immdata$data[[x]]$Clones)
                out <- apply(immdata$data[[x]], 1, function(y) {{
                    y['Barcode'] <- paste0(
                        x, "_", y['CDR3.aa'] ,"_", 1:y['Clones'], collapse = ";"
                    )
                    y
                }})
                as.data.frame(t(out))
            }})
            names(df) <- samples
            immdata$data <- df
            exdata <- expand_immdata(immdata)
            # Save it for later tests
            write.table(exdata, "{EXDATA_FILE}", sep="\\t", quote=F, row.names=F)

            for (sam in samples) {{
                d <- exdata[exdata$Sample == sam,]
                # if (sam == "A2-i129")
                #    write.table(d, "d.txt", sep="\\t", quote=F, row.names=F)
                if (all_cells[[sam]] != nrow(d)) {{
                    stop(paste0(
                        "# cells (", sam, ") in expanded data: ",
                        nrow(d), " != Original: ", all_cells[[sam]]
                    ))
                }}
            }}
        """

    @r_test
    def test_filter_expanded_immdata(self):
        return f"""
            exdata <- read.table("{EXDATA_FILE}", sep="\\t", header=T, stringsAsFactors=F)
            filtered <- filter_expanded_immdata(
                exdata,
                'as.numeric(sub("^.+_", "", Barcode)) %% 2 == 0'
            )
            # Save it for later tests
            write.table(filtered, "{FLDATA_FILE}", sep="\\t", quote=F, row.names=F)
            # We got a lot of singletons, which were filtered out
            if (nrow(filtered) >= nrow(exdata) / 2) {{
                stop(paste0(
                    "Filtered rows (", nrow(filtered), ") !< Original rows (",
                    nrow(exdata), ") / 2"
                ))
            }}
        """

    @r_test
    def test_filter_expanded_immdata2(self):
        return f"""
            exdata <- read.table("{EXDATA_FILE}", sep="\\t", header=T, stringsAsFactors=F)
            filtered <- filter_expanded_immdata(
                exdata,
                'as.numeric(sub("^.+_", "", Barcode)) %% 2 == 0',
                update_clones = TRUE
            )
            # Save it for later tests
            write.table(filtered, "{FLDATA_FILE}2", sep="\\t", quote=F, row.names=F)
            if (filtered$Clones[1] != 86) {{
                stop(paste0("Clones (", filtered$Clones[1], ") != 86"))
            }}
            if (filtered$Clones[1] != 86 || abs(filtered$Proportion[1] - 0.0711331679073615) > 1e-2) {{
                stop(
                    paste0(
                        "Clones (",
                        filtered$Clones[1],
                        ") != 86 or Proportion (",
                        filtered$Proportion[1],
                        ") != 0.0711331679073615"
                    )
                )
            }}
        """

    @r_test
    def test_immdata_from_expanded(self):
        return f"""
            fldata <- read.table("{FLDATA_FILE}", sep="\\t", header=T, stringsAsFactors=F)
            immdata <- immdata_from_expanded(fldata)
            write.table(immdata$data[["A2-i129"]], "{REDATA_FILE}", sep="\\t", quote=F, row.names=F)
            fl_samples <- unique(fldata$Sample)
            immdata_samples <- names(immdata$data)
            if (!setequal(fl_samples, immdata_samples)) {{
                stop(paste0(
                    "Filtered samples (",
                    paste0(fl_samples, collapse = ", "),
                    ") != immdata samples (",
                    paste0(immdata_samples, collapse = ", "), ")"
                ))
            }}
            # check if chones recalculation is correct
            d <- immdata$data[["A2-i129"]]
            if (d$Clones[1] != 86 || abs(d$Proportion[1] - 0.0711331679073615) > 1e-2) {{
                stop(
                    paste0(
                        "Clones (",
                        d$Clones[1],
                        ") != 86 or Proportion (",
                        d$Proportion[1],
                        ") != 0.0711331679073615"
                    )
                )
            }}
        """


if __name__ == "__main__":
    main(verbosity=2)
