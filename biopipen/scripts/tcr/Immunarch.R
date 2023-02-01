# Basic analysis and clonality
# TODO: How about TRA chain?
library(immunarch)
library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
library(tibble)

theme_set(theme_prism())

immfile = {{ in.immdata | quote }}
outdir = {{ out.outdir | quote }}
volume_by = {{ envs.volume_by | r }}
len_by = {{ envs.len_by | r }}
count_by = {{ envs.count_by | r }}
top_clone_marks = {{ envs.top_clone_marks | r }}
top_clone_by = {{ envs.top_clone_by | r }}
rare_clone_marks = {{ envs.rare_clone_marks | r }}
rare_clone_by = {{ envs.rare_clone_by | r }}
hom_clone_marks = {{ envs.hom_clone_marks | r }}
hom_clone_by = {{ envs.hom_clone_by | r }}
overlap_methods = {{ envs.overlap_methods | r }}
overlap_redim = {{ envs.overlap_redim | r }}
gu_by = {{ envs.gu_by | r }}
gu_top = {{ envs.gu_top | r }}
gua_methods = {{ envs.gua_methods | r }}
spect = {{ envs.spect | enumerate | dict | r }}
div_methods = {{ envs.div_methods | r }}
div_by = {{ envs.div_by | r }}
raref = {{ envs.raref | r }}
tracking_target = {{ envs.tracking_target | r }}
tracking_samples = {{ envs.tracking_samples | r }}
kmers_args = {{ envs.kmers | r: ignoreintkey=False }}

hom_clone_marks = hom_clone_marks[order(unlist(hom_clone_marks))]

# Convert {0 = "Sex", 1 = ["Status", "Sex"]}, to
# {Sex = "Sex", "Status_Sex": ["Status", "Sex"]}
# and ["Status", "Sex"] to
# {"Status_Sex": ["Status", "Sex"]}
norm_list = function(vec) {
    if (is.null(vec)) {
        return(vec)
    }
    if (!is.list(vec)) {
        # ["Status", "Sex"]
        key = paste(vec, sep = "_")
        vec = list(vec)
        names(vec) = key
        return(vec)
    }
    any_na = anyNA(as.numeric(names(vec)))
    if (any_na) {
        return(vec)
    }
    newnames = sapply(names(vec), function(x) {
        paste(vec[[x]], sep = "_")
    })
    names(vec) = newnames
    return(vec)
}

volume_by = norm_list(volume_by)
len_by = norm_list(len_by)
count_by = norm_list(count_by)
top_clone_by = norm_list(top_clone_by)
rare_clone_by = norm_list(rare_clone_by)
hom_clone_by = norm_list(hom_clone_by)
gu_by = norm_list(gu_by)
div_by = norm_list(div_by)
raref$by = norm_list(raref$by)

immdata = readRDS(immfile)
n_samples = length(immdata$data)

# volume
print("- All samples volume")
volume_dir = file.path(outdir, "volume")
dir.create(volume_dir, showWarnings = FALSE)

exp_vol = repExplore(immdata$data, .method = "volume")
png(file.path(volume_dir, "volume.png"), res = 300, width = 2000, height = 2000)
print(vis(exp_vol))
dev.off()

# volume_by
print("- Volume by")
for (name in names(volume_by)) {
    png(
        file.path(volume_dir, paste0("volume-", name, ".png")),
        res = 300,
        width = 2000,
        height = 2000
    )
    print(vis(exp_vol, .by = volume_by[[name]], .meta = immdata$meta))
    dev.off()
}

for (col in c("aa", "nt")) {
    # len
    print(paste0("- len by ", col))
    len_dir = file.path(outdir, paste0("len-", col))
    dir.create(len_dir, showWarnings = FALSE)

    exp_len = repExplore(immdata$data, .method = "len", .col = col)
    png(file.path(len_dir, "len.png"), res = 300, width = 2500, height = 2000)
    print(vis(exp_len))
    dev.off()

    # len_by
    print(paste0("- len_by by ", col))
    for (name in names(len_by)) {
        print(paste0("  ", name))
        png(
            file.path(len_dir, paste0("len_", col, "-", name, ".png")),
            res = 300,
            width = 2000,
            height = 2000
        )
        print(vis(exp_len, .by = len_by[[name]], .meta = immdata$meta))
        dev.off()
    }

    # length distribution
    # lend_dir = file.path(outdir, paste0("lend-", col))
    # dir.create(lend_dir, showWarnings = FALSE)

    # lendata = list()
    # for (sample in names(immdata$data)) {
    #     lendata[[sample]] = immdata$data[[sample]] %>%
    #         select(CDR3=paste0("CDR3.", col)) %>%
    #         mutate(Length = nchar(CDR3), Sample = sample) %>%
    #         select(Sample, Length)
    # }
    # lendata = do.call(bind_rows, lendata) %>%
    #     left_join(immdata$meta, by = "Sample")

    # p = ggplot(lendata, aes(x = Length))
    # png(file.path(lend_dir, "lend.png"), res = 300, width = 2500, height = 2000)
    # print(p + geom_histogram(bindwidth = .5))
    # dev.off()

    # # length distribution by
    # for (name in names(len_by)) {
    #     png(
    #         file.path(lend_dir, paste0("lend_", col, "-", name, ".png")),
    #         res = 300,
    #         width = 2000,
    #         height = 2000
    #     )
    #     print(p + geom_histogram(aes_string(fill=len_by[[name]]), alpha=.5, position="dodge"))
    #     dev.off()
    # }
}


# count
print("- All samples count")
count_dir = file.path(outdir, "count")
dir.create(count_dir, showWarnings = FALSE)

exp_count = repExplore(immdata$data, .method = "count")
png(file.path(count_dir, "count.png"), res = 300, width = 2000, height = 2000)
print(vis(exp_count))
dev.off()

# count_by
print("- Count by")
for (name in names(count_by)) {
    png(
        file.path(count_dir, paste0("count-", name, ".png")),
        res = 300,
        width = 2000,
        height = 2000
    )
    print(vis(exp_count, .by = count_by[[name]], .meta = immdata$meta))
    dev.off()
}


# top_clone
print("- Top clones")
imm_top = repClonality(immdata$data, .method = "top", .head = top_clone_marks)
top_dir = file.path(outdir, "top_clones")
dir.create(top_dir, showWarnings = FALSE)

png(
    file.path(top_dir, "top_clones.png"),
    res = 300,
    width = 2000,
    height = 2000
)
print(vis(imm_top))
dev.off()

for (name in names(top_clone_by)) {
    print(paste0("  by ", name))
    png(
        file.path(top_dir, paste0("top_clones-", name, ".png")),
        res = 300,
        width = 2000,
        height = 2000
    )
    print(vis(imm_top, .by = top_clone_by[[name]], .meta = immdata$meta))
    dev.off()
}


# rare_clone
print("- Rare clones")
imm_rare = repClonality(immdata$data, .method = "rare", .bound = rare_clone_marks)
rare_dir = file.path(outdir, "rare_clones")
dir.create(rare_dir, showWarnings = FALSE)

png(
    file.path(rare_dir, "rare_clones.png"),
    res = 300,
    width = 2000,
    height = 2000
)
print(vis(imm_rare))
dev.off()

for (name in names(rare_clone_by)) {
    print(paste0("  by ", name))
    png(
        file.path(rare_dir, paste0("rare_clones-", name, ".png")),
        res = 300,
        width = 2000,
        height = 2000
    )
    print(vis(imm_rare, .by = rare_clone_by[[name]], .meta = immdata$meta))
    dev.off()
}


# homeo_clone
print("- Homeo clones")
imm_hom = repClonality(
    immdata$data,
    .method = "homeo",
    .clone.types = hom_clone_marks
)
hom_dir = file.path(outdir, "homeo_clones")
dir.create(hom_dir, showWarnings = FALSE)

png(
    file.path(hom_dir, "hom_clones.png"),
    res = 300,
    width = 2000,
    height = 2000
)
print(vis(imm_hom))
dev.off()

for (name in names(hom_clone_by)) {
    print(paste0("  by ", name))
    png(
        file.path(hom_dir, paste0("hom_clones-", name, ".png")),
        res = 300,
        width = 2000,
        height = 2000
    )
    print(vis(imm_hom, .by = hom_clone_by[[name]], .meta = immdata$meta))
    dev.off()
}

# overlap
ov_dir = file.path(outdir, "overlap")
dir.create(ov_dir, showWarnings = FALSE)

print("- Overlap")
for (method in overlap_methods) {
    print(paste0("  ", method))
    ovpng = file.path(ov_dir, paste0("overlap-", method, ".png"))
    imm_ov = repOverlap(immdata$data, .method=method, .verbose=FALSE)
    png(ovpng, res=300, height = 2000, width = 2000)
    if (method == "public") {
        print(vis(imm_ov))
    } else {
        print(vis(imm_ov, .text.size=2))
    }
    dev.off()

    for (red in overlap_redim) {
        ovapng = file.path(ov_dir, paste0("overlapanalysis-", method, "-", red, ".png"))
        ova = tryCatch({
            repOverlapAnalysis(imm_ov, .method = red)
        }, error = function(e) NULL)
        # in case too few samples
        if (!is.null(ova)) {
            png(ovapng, res=300, height = 2000, width = 2500)
            print(vis(ova))
            dev.off()
        }
    }
}

# Gene usage
# https://immunarch.com/articles/web_only/v5_gene_usage.html
print("- Gene usage")
gu_dir = file.path(outdir, "gene_usage")
dir.create(gu_dir, showWarnings = FALSE)
imm_gu = geneUsage(immdata$data, "hs.trbv") %>%
    separate_rows(Names, sep=";") %>%
    group_by(Names) %>%
    summarise(across(everything(), ~ sum(., na.rm = TRUE)))
imm_gu = imm_gu %>%
    arrange(desc(rowSums(select(imm_gu, -"Names"), na.rm = TRUE)))
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
    print(paste0("  by ", name))
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
print("- Gene usage analysis")
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
    print(paste0("  ", gua_method))
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
print("- Spectratyping")
spect_dir = file.path(outdir, "spectratyping")
dir.create(spect_dir, showWarnings = FALSE)

for (sample in names(immdata$data)) {
    print(paste0("  ", sample))
    spect_sam_dir = file.path(spect_dir, sample)
    dir.create(spect_sam_dir, showWarnings = FALSE)
    for (idx in seq_along(spect)) {
        spec_obj = spectratype(
            immdata$data[[sample]],
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

print("- Diversity estimation")
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
    print(paste0("  ", div_method))
    met_dir = file.path(div_dir, div_method)
    dir.create(met_dir, showWarnings = FALSE)
    div = repDiversity(immdata$data, div_method)
    divpng = file.path(met_dir, paste0("diversity-1-.png"))
    png(divpng, res=300, height=2000, width=2000)
    print(plot_div(div, div_method))
    dev.off()

    for (name in names(div_by)) {
        divpng = file.path(met_dir, paste0("diversity-1-", name, ".png"))
        png(divpng, res=300, height=2000, width=2000)
        print(plot_div(div, div_method, .by=div_by[[name]], .meta=immdata$meta))
        dev.off()
    }

    tryCatch({
        # repFilter only supported in immunarch 0.6.7
        div1 = repDiversity(
            repFilter(immdata, .method = "by.clonotype", .query = list(Clones = morethan(1)))$data,
            div_method
        )
        divpng1 = file.path(met_dir, paste0("diversity-2-.png"))
        png(divpng1, res=300, height=2000, width=2000)
        print(plot_div(div1, div_method))
        dev.off()

        for (name in names(div_by)) {
            divpng = file.path(met_dir, paste0("diversity-2-", name, ".png"))
            png(divpng, res=300, height=2000, width=2000)
            print(plot_div(div1, div_method, .by=div_by[[name]], .meta=immdata$meta))
            dev.off()
        }

        div2 = repDiversity(
            repFilter(immdata, .method = "by.clonotype", .query = list(Clones = morethan(2)))$data,
            div_method
        )
        divpng2 = file.path(met_dir, paste0("diversity-3-.png"))
        png(divpng2, res=300, height=2000, width=2000)
        print(plot_div(div2, div_method))
        dev.off()

        for (name in names(div_by)) {
            divpng = file.path(met_dir, paste0("diversity-3-", name, ".png"))
            png(divpng, res=300, height=2000, width=2000)
            print(plot_div(div2, div_method, .by=div_by[[name]], .meta=immdata$meta))
            dev.off()
        }
    }, error=function(e) warning(as.character(e)))

}

# Rarefaction
print("- Rarefaction")
raref_dir = file.path(outdir, "raref")
dir.create(raref_dir, showWarnings = FALSE)

raref_x = NULL
raref_y = NULL
raref_log_x = NULL
raref_analysis = function(idata, sepname, get_max = FALSE) {
    raref_pms = raref  # copy the parameters
    raref_by = raref_pms$by
    raref_pms$by = NULL
    raref_align_x = raref_pms$align_x
    raref_pms$align_x = NULL
    raref_align_y = raref_pms$align_y
    raref_pms$align_y = NULL
    raref_log = raref_pms$log
    raref_pms$log = NULL
    raref_pms$.method = "raref"

    if (is.null(raref_pms$.verbose)) {
        raref_pms$.verbose = F
    }

    imm_raref = tryCatch({
        raref_pms$.data = idata$data
        do.call(repDiversity, raref_pms)
    }, error=function(e) {
        # https://github.com/immunomind/immunarch/issues/44
        valid_samples = c()
        for (sam in names(idata$data)) {
            raref_pms$.data = idata$data[sam]
            vsam = tryCatch({
                do.call(repDiversity, raref_pms)
                sam
            }, error=function(e) {
                warning(
                    paste("Rarefraction analysis failed for sample", sam, ":", as.character(e))
                )
                c()
            })
            valid_samples = c(valid_samples, vsam)
        }
        raref_pms$.data = idata$data[valid_samples]
        do.call(repDiversity, raref_pms)
    })
    rarefpng = file.path(raref_dir, paste0("raref-", sub("-$", "", sepname), ".png"))

    width = 1800 + ceiling(length(idata$data) / 20) * 500
    png(rarefpng, res=300, width=width, height=2000)
    p = vis(imm_raref) + xlab("Sample size (cells)")
    if (get_max) {
        raref_x <<- layer_scales(p)$x$range$range[2]
        raref_y <<- layer_scales(p)$y$range$range[2]
    } else {
        if (!is.null(raref_x)) {
            p = p + xlim(0, raref_x)
        }
        if (!is.null(raref_y)) {
            p = p + ylim(0, raref_y)
        }
    }
    print(p)
    dev.off()

    if (isTRUE(raref_log)) {
        rarefpng = file.path(raref_dir, paste0("raref-", sub("-$", "", sepname), "(log).png"))
        png(rarefpng, res=300, width=width, height=2000)
        p_log = vis(imm_raref, .log = TRUE) + xlab("Sample size (cells)")
        if (get_max) {
            raref_log_x <<- layer_scales(p_log)$x$range$range[2]
        } else {
            if (!is.null(raref_log_x)) {
                p_log = p_log + scale_x_log10(limits = c(1, 10 ^ raref_log_x))
            }
            if (!is.null(raref_y)) {
                p_log = p_log + ylim(0, raref_y)
            }
        }
        print(p_log)
        dev.off()
    }


    for (name in names(raref_by)) {
        print(paste0("  * by ", name))
        rfpng = file.path(raref_dir, paste0("raref-", sepname, name, ".png"))
        png(rfpng, res=300, width=2200, height=2000)
        p = vis(imm_raref, .by=raref_by[[name]], .meta=idata$meta) + xlab("Sample size (cells)")
        if (!is.null(raref_x)) {
            p = p + xlim(c(0, raref_x))
        }
        if (!is.null(raref_y)) {
            p = p + ylim(c(0, raref_y))
        }
        print(p)
        dev.off()

        if (isTRUE(raref_log)) {
            rfpng = file.path(raref_dir, paste0("raref-", sepname, name, "(log).png"))
            png(rfpng, res=300, width=2200, height=2000)
            p_log = vis(imm_raref, .by=raref_by[[name]], .meta=idata$meta, .log=TRUE) + xlab("Sample size (cells)")
            if (!is.null(raref_log_x)) {
                p_log = p_log + scale_x_log10(limits = c(1, 10 ^ raref_log_x))
            }
            if (!is.null(raref_y)) {
                p_log = p_log + ylim(c(0, raref_y))
            }
            print(p_log)
            dev.off()
        }
    }
}

print("  All samples")
raref_sep_by = raref$separate_by
raref$separate_by = NULL
raref_analysis(immdata, "ALL-", TRUE)
if (!is.null(raref_sep_by)) {
    sepvars = unique(immdata$meta[[raref_sep_by]])

    for (sepvar in sepvars) {
        print(paste0("  ", raref_sep_by, ": ", sepvar))
        q = list(include(sepvar))
        names(q) = raref_sep_by
        sepdata = repFilter(immdata, .method = "by.meta", .query = q)
        raref_analysis(sepdata, paste0(sepvar, "-"))
    }
}


# Clonotype tracking
print("- Clonotype tracking")
tracking_dir = file.path(outdir, "tracking")
dir.create(tracking_dir, showWarnings = FALSE)
for (name in names(tracking_target)) {
    print(paste0("  ", name))
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
                    select(-c(Clones.x, Clones.y)) %>%
                    ungroup()
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
print("- K-mer analysis")
kmer_dir = file.path(outdir, "kmer")
dir.create(kmer_dir, showWarnings = FALSE)
for (k in names(kmers_args)) {
    print(paste0("  k=", k))
    k_dir = file.path(kmer_dir, paste0("kmer_", k))
    dir.create(k_dir, showWarnings = FALSE)
    kmer_args = kmers_args[[k]]
    k = as.integer(k)

    kmers = getKmers(immdata$data, k)
    head = kmer_args$head
    if (is.null(head)) { head = 10 }
    position = kmer_args$position
    if (is.null(position)) { position = "stack" }
    logg = kmer_args$log
    if (is.null(logg)) { logg = FALSE }

    for (h in head) {
        print(paste0("    head: ", h))
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
        print(paste0("    motif: ", mot))
        mot_dir = file.path(k_dir, paste0("motif_", mot))
        dir.create(mot_dir, showWarnings = FALSE)
        # multiple samples not supported as of 0.6.7
        for (sample in names(immdata$data)) {
            kmers_sample = getKmers(immdata$data[[sample]], k)
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
