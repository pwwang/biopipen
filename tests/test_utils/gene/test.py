from unittest import TestCase, main
from biopipen.utils.gene import gene_name_conversion

from pathlib import Path

# from biopipen.core.filters import r
from biopipen.core.testing import r_test


class TestGenePy(TestCase):

    def test_gene_name_conversion(self):
        # gene_name_conversion(genes, species, infmt, outfmt, notfound)
        out = gene_name_conversion(
            [
                "ENSG00000230373",
                "ENSG00000236269",
                "ENSG00000227232.5",
                "ENSG00000227232.5",
            ],
            species="human",
            infmt="ensg",
            outfmt="symbol",
            suppress_messages=True,
        )
        self.assertListEqual(
            out["symbol"].tolist(),
            ['GOLGA6L3P', 'NA', 'WASH7P', 'WASH7P'],
        )

    def test_gene_name_conversion1(self):
        out = gene_name_conversion(
            ['NM_003466', 'CDK2', "695", '1320_at', 'Q08345'],
            species="human",
            infmt="refseq,symbol,entrezgene,reporter,uniprot",
            outfmt="symbol",
            notfound="ignore",
            suppress_messages=True,
        )
        self.assertListEqual(
            out["symbol"].tolist(),
            ["PAX8", "CDK2", "BTK", "PTPN21", "DDR1"]
        )

    def test_gene_name_conversion_keep_dup(self):
        out = gene_name_conversion(
            [
                "ENSG00000230373",
                "ENSG00000236269",
                "ENSG00000227232.5",
            ],
            species="human",
            infmt="ensg",
            outfmt="symbol",
            notfound="ignore",
            dup=";",
            suppress_messages=True,
        )
        self.assertIn("GOLGA6L3P;GOLGA6L17P", out["symbol"].tolist())

    def test_gene_name_conversion_skip(self):
        out = gene_name_conversion(
            [
                "ENSG00000230373",
                "ENSG00000236269",
                "ENSG00000227232.5",
            ],
            species="human",
            infmt="ensg",
            outfmt="symbol",
            notfound="skip",
            suppress_messages=True,
        )
        self.assertNotIn("NA", out["symbol"].tolist())


class TestGeneR(TestCase):

    SOURCE_FILE = Path(__file__).parent.parent.parent.parent.joinpath(
        "biopipen", "utils", "gene.R"
    )

    @r_test
    def test_gene_name_conversion(self):
        return """
            ensgs = c(
                "ENSG00000230373",
                "ENSG00000236269",
                "ENSG00000227232.5",
                "ENSG00000227232.5"
            )
            symbols = gene_name_conversion(
                ensgs, "ensg", "symbol", suppress_messages = TRUE
            )
            expect(length(symbols$symbol) == 4, "Wrong number of symbols returned")
            expect(
                symbols$symbol[1] == "GOLGA6L3P",
                "Expected GOLGA6L3P as 1st symbol", "Got: ", symbols$symbol[1]
            )
            expect(
                is.na(symbols$symbol[2]),
                "Expected NA as 2nd symbol", "Got: ", symbols$symbol[2]
            )
            expect(
                symbols$symbol[3] == "WASH7P",
                "Expected WASH7P as 3rd symbol", "Got: ", symbols$symbol[3]
            )
            expect(
                symbols$symbol[4] == "WASH7P",
                "Expected WASH7P as 4th symbol", "Got: ", symbols$symbol[4]
            )
        """

    @r_test
    def test_gene_name_conversion_keep_dup(self):
        return """
            ensgs = c(
                "ENSG00000230373",
                "ENSG00000236269",
                "ENSG00000227232.5"
            )
            symbols = gene_name_conversion(
                ensgs, "ensg", "symbol", ";", suppress_messages = TRUE
            )
            if (!"GOLGA6L3P;GOLGA6L17P" %in% symbols$symbol) {
                stop("GOLGA6L3P;GOLGA6L17P not found")
            }
        """

    @r_test
    def test_gene_name_conversion_skip(self):
        return """
            ensgs = c(
                "ENSG00000230373",
                "ENSG00000236269",
                "ENSG00000227232.5"
            )
            symbols = gene_name_conversion(
                ensgs, "ensg", "symbol", notfound="skip", suppress_messages = TRUE
            )
            if (anyNA(symbols$symbol)) {
                stop("Missing queries were not skipped")
            }
        """


if __name__ == "__main__":
    main(verbosity=2)
