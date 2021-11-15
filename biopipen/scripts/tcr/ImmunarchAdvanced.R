# Advanced analysis
library(immunarch)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)

theme_set(theme_prism())

immfile = {{ in.immdata | quote }}
outdir = {{ out.outdir | quote }}
gu_by = {{ envs.gu_by | r }}
gu_top = {{ envs.gu_top | r }}
gua_methods = {{ envs.gua_methods | r }}
spect = {{ envs.spect | enumerate | dict | r }}
div_methods = {{ envs.div_methods | r }}
div_by = {{ envs.div_by | r }}
raref_by = {{ envs.raref_by | r }}
tracking_target = {{ envs.tracking_target | r }}
tracking_samples = {{ envs.tracking_samples | r }}
kmers_args = {{ envs.kmers | r: ignoreintkey=False }}

vec_to_list = function(vec) {
    key = paste(vec, sep = "_")
    vec = list(vec)
    names(vec) = key
    vec
}

if (!is.list(gu_by)) { gu_by = vec_to_list(gu_by) }
if (!is.list(div_by)) { gu_by = vec_to_list(div_by) }
if (!is.list(raref_by)) { gu_by = vec_to_list(raref_by) }

immdata = readRDS(immfile)
n_samples = length(immdata$data)

# Gene usage
# https://immunarch.com/articles/web_only/v5_gene_usage.html

gu_dir = file.path(outdir, "gene_usage")
dir.create(gu_dir, showWarnings = FALSE)
imm_gu = geneUsage(immdata$single, "hs.trbv") %>%
    separate_rows(Names, sep=";") %>%
    group_by(Names) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>%
    arrange(desc(rowSums(select(., -"Names")), na.rm = TRUE))
imm_gu_top = imm_gu %>% head(gu_top)

class(imm_gu) = append("immunr_gene_usage", class(imm_gu))
class(imm_gu_top) = append("immunr_gene_usage", class(imm_gu_top))
gupng = file.path(gu_dir, "gene_usage.png")
width = 2000 + ceiling(n_samples / 15) * 500
png(gupng, res=300, height=1500, width=width)
print(vis(imm_gu_top))
dev.off()

# png(
#     file.path(gu_dir, paste0("gene_usage-_grid.png")),
#     res = 300,
#     width = 2500,
#     height = 1500
# )
# print(vis(imm_gu, .grid = TRUE))
# dev.off()

for (name in names(gu_by)) {
    png(
        file.path(gu_dir, paste0("gene_usage-", name, ".png")),
        res = 300,
        width = 2500,
        height = 1500
    )
    print(vis(imm_gu_top, .by = gu_by[[name]], .meta = immdata$meta))
    dev.off()
}

# Gene usage analysis
gua_dir = file.path(outdir, "gene_usage_analysis")
dir.create(gua_dir, showWarnings = FALSE)
gua_method_names = c(
    js = "Jensen-Shannon Divergence",
    cor = "Correlation",
    cosine = "Cosine similarity",
    pca = "Principal component analysis",
    mds = "Multi-dimensional scaling",
    tsne = "T-Distributed Stochastic Neighbor Embedding"
)
for (gua_method in gua_methods) {
    imm_gua = geneUsageAnalysis(imm_gu, .method = gua_method, .verbose = FALSE)
    p = vis(
        imm_gua,
        .title = gua_method_names[gua_method],
        .leg.title = toupper(gua_method),
        .text.size = 1.5
    )
    png(
        file.path(gua_dir, paste0("gene_usage_analysis-", gua_method, ".png")),
        res = 300,
        width = 2000,
        height = 2000
    )
    print(p)
    dev.off()
}

# Spectratyping
spect_dir = file.path(outdir, "spectratyping")
dir.create(spect_dir, showWarnings = FALSE)

for (sample in names(immdata$data)) {
    spect_sam_dir = file.path(spect_dir, sample)
    dir.create(spect_sam_dir, showWarnings = FALSE)
    for (idx in seq_along(spect)) {
        spec_obj = spectratype(
            immdata$single[[sample]],
            .quant = spect[[idx]]$quant,
            .col = spect[[idx]]$col
        )
        png(
            file.path(spect_sam_dir, paste0("spectratyping-", idx, ".png")),
            res = 150,
            width = 1000,
            height = 1000
        )
        print(vis(spec_obj))
        dev.off()
    }
}


# Diversity estimation
# https://immunarch.com/articles/web_only/v6_diversity.html
div_dir = file.path(outdir, "diversity")
dir.create(div_dir, showWarnings = FALSE)

plot_div = function(div, method, ...) {
    if (method != "gini") {
        do.call(vis, list(div, ...))
    } else {
        ginidiv = as.data.frame(div) %>%
            rownames_to_column("Sample") %>%
            rename(`Gini-coefficient`=V1)
        ggplot(ginidiv) +
            geom_col(aes(x=Sample, y=`Gini-coefficient`, fill=Sample)) +
            theme(
                legend.position="none",
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
            )
    }
}

for (div_method in div_methods) {
    met_dir = file.path(div_dir, div_method)
    dir.create(met_dir, showWarnings = FALSE)
    div = repDiversity(immdata$data, div_method)
    divpng = file.path(met_dir, paste0("diversity-1-.png"))
    png(divpng, res=300, height=2000, width=2000)
    print(plot_div(div, div_method))
    dev.off()

    # repFilter only supported in immunarch 0.6.7
    div1 = repDiversity(
        repFilter(immdata, .method = "by.clonotype", .query = list(Clones = morethan(1)))$data,
        div_method
    )
    divpng1 = file.path(met_dir, paste0("diversity-2-.png"))
    png(divpng1, res=300, height=2000, width=2000)
    print(plot_div(div1, div_method))
    dev.off()

    div2 = repDiversity(
        repFilter(immdata, .method = "by.clonotype", .query = list(Clones = morethan(2)))$data,
        div_method
    )
    divpng2 = file.path(met_dir, paste0("diversity-3-.png"))
    png(divpng2, res=300, height=2000, width=2000)
    print(plot_div(div2, div_method))
    dev.off()

    for (name in names(div_by)) {
        divpng = file.path(met_dir, paste0("diversity-1-", name, ".png"))
        png(divpng, res=300, height=2000, width=2000)
        print(plot_div(div, div_method, .by=div_by[[name]], .meta=immdata$meta))
        dev.off()

        divpng = file.path(met_dir, paste0("diversity-2-", name, ".png"))
        png(divpng, res=300, height=2000, width=2000)
        print(plot_div(div1, div_method, .by=div_by[[name]], .meta=immdata$meta))
        dev.off()

        divpng = file.path(met_dir, paste0("diversity-3-", name, ".png"))
        png(divpng, res=300, height=2000, width=2000)
        print(plot_div(div2, div_method, .by=div_by[[name]], .meta=immdata$meta))
        dev.off()
    }
}

# Rarefaction
raref_dir = file.path(outdir, "raref")
dir.create(raref_dir, showWarnings = FALSE)
imm_raref = tryCatch({
    repDiversity(immdata$data, "raref", .verbose = F)
}, error=function(e) {
    # https://github.com/immunomind/immunarch/issues/44
    valid_samples = c()
    for (sam in names(immdata$data)) {
        vsam = tryCatch({
            repDiversity(immdata$data[sam], "raref", .verbose = F)
            sam
        }, error=function(e) {c()})
        valid_samples = c(valid_samples, vsam)
    }
    repDiversity(immdata$data[valid_samples], "raref", .verbose = F)
})
rarefpng = file.path(raref_dir, "raref-.png")

width = 1700 + ceiling(n_samples / 15) * 500
png(rarefpng, res=300, width=width, height=2000)
print(vis(imm_raref))
dev.off()

for (name in names(raref_by)) {
    rfpng = file.path(raref_dir, paste0("raref-", name, ".png"))
    png(rfpng, res=300, width=width, height=2000)
    print(vis(imm_raref, .by=raref_by[[name]], .meta=immdata$data))
    dev.off()
}


# Clonotype tracking
tracking_dir = file.path(outdir, "tracking")
dir.create(tracking_dir, showWarnings = FALSE)
for (name in names(tracking_target)) {
    target = tracking_target[[name]]
    samples = tracking_samples[[name]]
    if (is.null(samples)) {
        samples = names(immdata$data)
    }
    if (length(samples) == 1) {
        stop(paste0("Cannot track clonotypes for only one sample: ", samples))
    }
    if (is.list(target)) {
        target = list(names(target), unname(unlist(target)))
    }
    if (target[[1]] == "TOP") {
        top_clones = NULL
        for (sample in samples) {
            if (is.null(top_clones)) {
                top_clones = immdata$data[[sample]]
            } else {
                top_clones = inner_join(top_clones, immdata$data[[sample]], by = "CDR3.aa") %>%
                    rowwise() %>%
                    mutate(Clones=sum(Clones.x, Clones.y)) %>%
                    select(-c(Clones.x, Clones.y))

            }
        }
        target = top_clones %>% slice_max(Clones, n=target[[2]]) %>% pull(CDR3.aa)
    }

    imm_tracking = trackClonotypes(immdata$data, target, .col = "aa")
    tracking_png = file.path(tracking_dir, paste0("tracking_", name, ".png"))
    png(tracking_png, res=300, height=2000, width=3000)
    print(vis(imm_tracking, .order = samples))
    dev.off()
}

# K-mer analysis
kmer_dir = file.path(outdir, "kmer")
dir.create(kmer_dir, showWarnings = FALSE)
for (k in names(kmers_args)) {
    k_dir = file.path(kmer_dir, paste0("kmer_", k))
    dir.create(k_dir, showWarnings = FALSE)
    kmer_args = kmers_args[[k]]
    k = as.integer(k)

    kmers = getKmers(immdata$single, k)
    head = kmer_args$head
    if (is.null(head)) { head = 10 }
    position = kmer_args$position
    if (is.null(position)) { position = "stack" }
    logg = kmer_args$log
    if (is.null(logg)) { logg = FALSE }

    for (h in head) {
        kmerpng = file.path(
            k_dir,
            paste0("head_", h, ".position_", position, ".log_", logg, ".png")
        )
        width = 200 * h + ceiling(n_samples / 15 - 1) * 500
        png(kmerpng, res=300, height=2000, width=width)
        print(vis(kmers, .head=h, .position=position, .log=logg))
        dev.off()
    }

    # motif analysis
    for (mot in kmer_args$motif) {
        mot_dir = file.path(k_dir, paste0("motif_", mot))
        dir.create(mot_dir, showWarnings = FALSE)
        # multiple samples not supported as of 0.6.7
        for (sample in names(immdata$data)) {
            kmers_sample = getKmers(immdata$single[[sample]], k)
            motif = kmer_profile(kmers_sample, mot)
            motpng = file.path(mot_dir, paste0("motif_", sample, ".png"))
            png(motpng, res=300, width=2000, height=2000)
            print(vis(motif) + ggtitle(subtitle = paste("Method:", mot), label = sample))
            dev.off()
            motseqpng = file.path(mot_dir, paste0("motif_", sample, ".seq.png"))
            png(motseqpng, res=300, width=2000, height=2000)
            print(vis(motif, .plot = "seq") + ggtitle(subtitle = paste("Method:", mot), label = sample))
            dev.off()
        }
    }
}
