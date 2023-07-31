source("{{biopipen_dir}}/utils/misc.R")
# Basic analysis and clonality
# TODO: How about TRA chain?
library(rlang)
library(immunarch)  # 0.8.0 or 0.9.0 tested
library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
library(tibble)

theme_set(theme_prism())

immfile = {{ in.immdata | quote }}
outdir = {{ out.outdir | quote }}
mutaters = {{ envs.mutaters | r }}
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

immdata = readRDS(immfile)
if (!is.null(mutaters) && length(mutaters) > 0) {
    expr = list()
    for (key in names(mutaters)) {
        expr[[key]] = parse_expr(mutaters[[key]])
    }
    immdata$meta = mutate(immdata$meta, !!!expr)
}

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
    # lendata = do_call(bind_rows, lendata) %>%
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
    png(ovpng, res=100, height = 400 + 40 * n_samples, width = 800 + 40 * n_samples)
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
            png(ovapng, res=100, height = 400 + 40 * n_samples, width = 800 + 40 * n_samples)
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
{% include biopipen_dir + "/scripts/tcr/Immunarch-diversity.R" %}

# Clonotype tracking
{% include biopipen_dir + "/scripts/tcr/Immunarch-tracking.R" %}

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
