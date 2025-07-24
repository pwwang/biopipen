"""Convert Seurat objects to AnnData format back and forth.

Need R and R packages Seurat, SeuratDisk and biopipen.utils.R installed.
"""


def convert_seurat_to_anndata(
    input_file,
    output_file,
    assay=None,
    subset=None,
    rscript="Rscript",
):
    """Convert Seurat object to AnnData format.

    Args:
        input_file (str): Path to the input Seurat RDS or qs/qs2 file.
        output_file (str): Path to the output AnnData H5AD file.
        assay (str): Name of the assay to use in the Seurat object.
        subset (str): An R expression to subset the Seurat object to convert.
        rscript (RScript): R script executor.
    """
    from biopipen.utils.misc import run_command

    script = f"""
        library(biopipen.utils)

        assay <- {repr(assay) if assay else 'NULL'}
        subset <- {repr(subset) if subset else 'NULL'}

        ConvertSeuratToAnnData(
            "{input_file}", "{output_file}", assay = assay, subset = subset
        )
    """

    # Save the script to a temporary file
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(suffix=".R", delete=False) as temp_script:
        temp_script.write(script.encode('utf-8'))
        temp_script_path = temp_script.name

    # Run the R script using the provided Rscript command
    cmd = [rscript, temp_script_path]
    run_command(cmd, fg=True)


def convert_anndata_to_seurat(
    input_file,
    output_file,
    assay=None,
    rscript="Rscript",
):
    """Convert AnnData object to Seurat format.

    Args:
        input_file (str): Path to the input AnnData H5AD file.
        output_file (str): Path to the output Seurat RDS or qs/qs2 file.
        assay (str): Name of the assay to use in the Seurat object.
        rscript (RScript): R script executor.
    """
    from biopipen.utils.misc import run_command

    script = f"""
        library(biopipen.utils)

        assay <- {repr(assay) if assay else 'NULL'}

        ConvertAnnDataToSeurat(
            "{input_file}", "{output_file}", assay = assay
        )
    """

    # Save the script to a temporary file
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(suffix=".R", delete=False) as temp_script:
        temp_script.write(script.encode('utf-8'))
        temp_script_path = temp_script.name

    # Run the R script using the provided Rscript command
    cmd = [rscript, temp_script_path]
    run_command(cmd, fg=True)
