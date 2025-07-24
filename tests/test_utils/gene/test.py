from unittest import TestCase, main
from biopipen.utils.gene import gene_name_conversion

# from pathlib import Path

# from biopipen.core.filters import r
# from biopipen.core.testing import r_test


class TestGenePy(TestCase):

    def test_gene_name_conversion(self):
        # gene_name_conversion(genes, species, infmt, outfmt, notfound)
        out = gene_name_conversion(
            [
                # "ENSG00000230373",
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
            # ['GOLGA6L17P', 'NA', 'WASH7P', 'WASH7P'],
            ['NA', 'WASH7P', 'WASH7P'],
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
                # "ENSG00000230373",
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
        self.assertIn("WASH7P", out["symbol"].tolist())

    def test_gene_name_conversion_skip(self):
        out = gene_name_conversion(
            [
                # "ENSG00000230373",
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


if __name__ == "__main__":
    main(verbosity=2)
