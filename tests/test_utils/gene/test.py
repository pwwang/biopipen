from unittest import TestCase, main
from biopipen.utils.gene import gene_name_conversion


class TestGene(TestCase):

    def test_gene_name_conversion1(self):
        # gene_name_conversion(genes, species, infmt, outfmt, notfound)
        out = gene_name_conversion(
            ['1255_g_at', '1316_at', '1320_at'],
            "human",
            "reporter",
            "symbol",
            "ignore"
        )
        self.assertListEqual(out["symbol"].tolist(), ["GUCA1A", "THRA", "PTPN21"])

    def test_gene_name_conversion2(self):
        out = gene_name_conversion(
            ['NM_003466', 'CDK2', "695", '1320_at', 'Q08345'],
            "human",
            "refseq,symbol,entrezgene,reporter,uniprot",
            "symbol",
            "ignore"
        )
        self.assertListEqual(
            out["symbol"].tolist(),
            ["PAX8", "CDK2", "BTK", "PTPN21", "DDR1"]
        )


if __name__ == "__main__":
    main(verbosity=2)
