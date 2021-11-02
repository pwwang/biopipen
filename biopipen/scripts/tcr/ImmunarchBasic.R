# Basic analysis and clonality
library(immunarch)
library(ggplot2)
library(ggprism)
library(dplyr)

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

hom_clone_marks = hom_clone_marks[order(unlist(hom_clone_marks))]

vec_to_list = function(vec) {
    key = paste(vec, sep = "_")
    vec = list(vec)
    names(vec) = key
    vec
}

if (!is.list(volume_by)) { volume_by = vec_to_list(volume_by) }
if (!is.list(len_by)) { len_by = vec_to_list(len_by) }
if (!is.list(count_by)) { count_by = vec_to_list(count_by) }
if (!is.list(top_clone_by)) { top_clone_by = vec_to_list(top_clone_by) }
if (!is.list(rare_clone_by)) { rare_clone_by = vec_to_list(rare_clone_by) }
if (!is.list(hom_clone_by)) { hom_clone_by = vec_to_list(hom_clone_by) }

immdata = readRDS(immfile)

# volume
volume_dir = file.path(outdir, "volume")
dir.create(volume_dir, showWarnings = FALSE)

exp_vol = repExplore(immdata$data, .method = "volume")
png(file.path(volume_dir, "volume.png"), res = 300, width = 2000, height = 2000)
print(vis(exp_vol))
dev.off()

# volume_by
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
    len_dir = file.path(outdir, paste0("len-", col))
    dir.create(len_dir, showWarnings = FALSE)

    exp_len = repExplore(immdata$single, .method = "len", .col = col)
    png(file.path(len_dir, "len.png"), res = 300, width = 2500, height = 2000)
    print(vis(exp_len))
    dev.off()

    # len_by
    for (name in names(len_by)) {
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
count_dir = file.path(outdir, "count")
dir.create(count_dir, showWarnings = FALSE)

exp_count = repExplore(immdata$data, .method = "count")
png(file.path(count_dir, "count.png"), res = 300, width = 2000, height = 2000)
print(vis(exp_count))
dev.off()

# count_by
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

for (method in overlap_methods) {
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
