from biopipen.core.proc import Proc
from biopipen.core.config import config
from biopipen.ns.scrna import Slingshot as Slingshot_
from biopipen.core.testing import get_pipeline


class PrepareData(Proc):
    lang = config.lang.rscript
    input = "var"
    input_data = ["slingshot_data"]
    output = "outfile:file:{{in.var}}.qs"
    script = """
        library(Seurat)
        library(biopipen.utils)

        means <- rbind(
            # non-DE genes
            matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
                ncol = 300, byrow = TRUE),
            # early deactivation
            matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
            # late deactivation
            matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
            # early activation
            matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
            # late activation
            matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
            # transient
            matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50),
                ncol = 300, byrow = TRUE)
        )
        counts <- apply(means,2,function(cell_means){
            total <- rnbinom(1, mu = 7500, size = 4)
            rmultinom(1, total, cell_means)
        })
        rownames(counts) <- paste0('G',1:750)
        colnames(counts) <- paste0('c',1:300)

        sobj <- CreateSeuratObject(counts = counts, project = "Slingshot")
        sobj <- subset(sobj, subset = nFeature_RNA >= 3 & nFeature_RNA >= 10)
        sobj <- NormalizeData(sobj)
        sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 500)
        sobj <- ScaleData(sobj, features = rownames(sobj))
        sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
        sobj <- RunUMAP(sobj, dims = 1:10)
        sobj <- FindNeighbors(sobj, dims = 1:10)
        sobj <- FindClusters(sobj, resolution = 0.5)

        save_obj(sobj, "{{out.outfile}}")
    """


class Slingshot(Slingshot_):
    requires = PrepareData


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareData)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
