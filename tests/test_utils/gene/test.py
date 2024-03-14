# TODO: use unittest
from biopipen.utils.gene import gene_name_conversion


def run_except(genes, infmt, outfmt, notfound, expect, species):
    try:
        gene_name_conversion(genes, species, infmt, outfmt, notfound)
    except expect:
        pass
    else:
        raise ValueError("Test Failed!")


def run_expect(genes, infmt, outfmt, notfound, expect, species):
    out = gene_name_conversion(genes, species, infmt, outfmt, notfound)
    assert out[outfmt].tolist() == expect


def run(genes, infmt, outfmt, notfound, expect, species="human"):
    print(f">>> TESTING {genes} - {infmt} - {outfmt} - {notfound} ")
    if isinstance(expect, Exception):
        run_except(genes, infmt, outfmt, notfound, expect, species)

    else:
        run_expect(genes, infmt, outfmt, notfound, expect, species)
    print(">>> PASSED")
    print(">>> ")


if __name__ == "__main__":
    run(
        ['1255_g_at', '1316_at', '1320_at'],
        infmt="reporter",
        outfmt="symbol",
        notfound="ignore",
        expect=["GUCA1A", "THRA", "PTPN21"]
    )
    run(
        ['NM_003466', 'CDK2', "695", '1320_at', 'Q08345'],
        infmt="refseq,symbol,entrezgene,reporter,uniprot",
        outfmt="symbol",
        notfound="ignore",
        expect=["PAX8", "CDK2", "BTK", "PTPN21", "DDR1"]
    )
