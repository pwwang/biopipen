# biopipen + the Parse Biosciences Pipeline (`split-pipe`)
#
# Build locally with the WSL container CLI (`wslc`). The build context is the biopipen
# checkout with the extracted Parse Biosciences Pipeline package inside it; `PARSE_SRC`
# points at the package directory in the context (default:
# `ParseBiosciences-Pipeline.1.9.0`):
#
#   wslc build -f docker/splitpipe_pipeline/Dockerfile \
#       -t biopipen/splitpipe:1.9.0 \
#       D:\path\to\biopipen
#
#   # or, when the package is extracted somewhere else in the checkout:
#   wslc build -f ... --build-arg PARSE_SRC=path/to/ParseBiosciences-Pipeline.1.9.0 -t ... D:\path\to\biopipen
#
# The image installs the checkout itself (editable, from /biopipen), so it embeds the
# development version of biopipen, including the `splitpipe` namespace.
#
# The image embeds the Parse Biosciences Pipeline (research-use-only software).
#
# The build is multi-stage, like immunopipe's: the `builder` stage installs everything
# (and trims it via cleanup.py), and the final stage copies just the two conda envs
# into a clean slim base — so build-time layers (conda package caches, tarballs,
# compiler toolchain) never end up in the published image.

FROM continuumio/miniconda3:latest AS builder

ARG PARSE_SRC=ParseBiosciences-Pipeline.1.9.0
ARG PYTHON_VERSION=3.12.8

COPY ${PARSE_SRC}/ /opt/parse/
RUN test -f /opt/parse/setup.py && test -f /opt/parse/install_dependencies_conda.sh || { \
        echo "ERROR: PARSE_SRC ('${PARSE_SRC}') does not look like an extracted Parse Biosciences Pipeline package"; \
        exit 1; \
    }

# conda.anaconda.org occasionally drops a package download mid-solve (seen on
# r-scpred from bioconda, HTTP 000). Retry a few times instead of failing the
# whole build; this applies to the base env below and the Parse env above.
RUN conda config --set remote_max_retries 10 && \
    conda config --set remote_connect_timeout_secs 30 && \
    conda config --set remote_read_timeout_secs 120

RUN conda create -y -n spipe -c conda-forge python==${PYTHON_VERSION} && \
    conda clean --all --yes


# The biopipen runtime environment (python 3.11, as in the biopipen/base image).
# A newer python resolves a quart release (> 0.22) whose API the released
# pipen-board cannot import yet; see the pipen-board fix for quart >= 0.23.
RUN conda install -y -n base -c conda-forge -c bioconda -c pwwang \
        python=3.11 nodejs=22 r-base=4.5 r-biopipen.utils && \
    conda clean --all --yes

# Accept Anaconda's channel terms of service non-interactively: Parse's own
# install script installs packages from the default channels (pkgs/main, pkgs/r),
# which conda refuses to use until the ToS is accepted.
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

WORKDIR /opt/parse
RUN conda run -n spipe bash ./install_dependencies_conda.sh -i -y && \
    conda run -n spipe pip install --no-cache-dir ./ && \
    conda clean --all --yes

# split-pipe is linked into the base env's bin dir; the Parse tools (STAR,
# samtools, pigz, bgzip, tabix, ...) come from the `spipe` env via PATH. They are
# deliberately NOT symlinked: conda's STAR is a wrapper script that execs its
# sibling STAR-plain, which only resolves in its real directory. The profile.d
# hook restores the PATH entry in login shells, whose PATH is rebuilt by conda's
# own profile script.
RUN ln -sf /opt/conda/envs/spipe/bin/split-pipe /opt/conda/bin/split-pipe && \
    printf '%s\n' 'export PATH="$PATH:/opt/conda/envs/spipe/bin"' \
        > /etc/profile.d/zz-spipe-path.sh && \
    chmod 0644 /etc/profile.d/zz-spipe-path.sh

# Parse's helper tools (STAR, samtools, pigz, bgzip, tabix, ...) live in the
# `spipe` env; append it to PATH so `split-pipe` can find them at run time,
# while the base env keeps python, pipen and R first in the search order.
ENV PATH=${PATH}:/opt/conda/envs/spipe/bin

# Trim the installed envs before the final stage copies them over: the conda
# package cache alone is ~4.7 GB. cleanup.py removes what the manifest lists.
COPY docker/splitpipe_pipeline/cleanup.py docker/splitpipe_pipeline/cleanup.list /tmp/cleanup/
RUN conda clean --all --yes && python /tmp/cleanup/cleanup.py

# The biopipen checkout, installed editable (and copied LAST, so that source
# edits rebuild in about a minute instead of re-running the Parse install).
COPY pyproject.toml README.md /biopipen/
COPY biopipen/ /biopipen/biopipen/
COPY docker/ /biopipen/docker/
# --upgrade keeps this recipe pin-free while still pulling newer releases of
# the stack (e.g. xqute's SpecPath._mounted fix) on a rebuild: without it pip
# reports 'requirement already satisfied' against the cached builder layers.
RUN cd /biopipen && pip install --no-cache-dir --upgrade -e ".[runinfo,log2file]" && \
    pipen report update

# TCR/BCR support: the immune references. Measured 121 MB in total: the IgBlast
# 1.22.0 binaries (116 MB) plus the human and mouse TCR/BCR germline databases
# (1.9 MB, 200 files, families human_tcr, human_bcr, mouse_tcr, mouse_bcr,
# transgenic_mouse). `-i` installs the default germline databases; `-I` would also
# fetch IMGT (larger, and its site is frequently down), so it stays out.
# The installer writes immune_dir into the source tree, while split-pipe reads it
# from the installed package (`os.path.dirname(splitpipe.__file__)`), so the IgBlast
# tree is copied into site-packages. A second `pip install ./` would do the same,
# but it rebuilds the Cython extension and that link fails in this image
# (`x86_64-conda-linux-gnu-ld: cannot find -lgcc`). Do not probe with
# `import splitpipe` here: the WORKDIR is /opt/parse, so that imports the source tree,
# whose Cython extension is not built (`cy_utils C-extension not found`). The shell
# check fails the build if no BLAST index landed in site-packages. The site-packages
# path comes from the env's own interpreter: globbing `lib/python3.*/site-packages`
# also matches a stray `python3.1` directory, which yields two paths.
# The builder stage lacks wget (tarball fetch) and unzip (germline archives);
# without unzip the installer prints "Problem installing IgBlast" and STILL exits 0,
# which is why the log check below is required.
RUN apt-get update && \
    apt-get install -y --no-install-recommends wget unzip && \
    rm -rf /var/lib/apt/lists/*

RUN conda run -n spipe bash ./install_immune_dbs.sh -i -y && \
    echo "STEP 1 installer exited 0" && \
    grep -q "Finished successfully" splitpipe/immune_dir/install_IgBlast.log && \
    echo "STEP 2 installer log confirms success" && \
    SP=$(/opt/conda/envs/spipe/bin/python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])") && \
    echo "STEP 3 site-packages: $SP" && \
    mkdir -p "$SP/splitpipe/immune_dir" && \
    echo "STEP 4 installed package has immune_dir" && \
    cp -a /opt/parse/splitpipe/immune_dir/IgBlast "$SP/splitpipe/immune_dir/" && \
    echo "STEP 5 IgBlast copied" && \
    echo "STEP 6 BLAST indexes: $(ls "$SP/splitpipe/immune_dir/IgBlast/tr_database"/*.ndb 2>/dev/null | wc -l)" && \
    test "$(ls "$SP/splitpipe/immune_dir/IgBlast/tr_database"/*.ndb 2>/dev/null | wc -l)" -gt 0

RUN conda clean --all --yes || true


# ---------------------------------------------------------------------------
# Final stage: a clean slim base carrying only the two conda envs, the Parse
# install, and the biopipen checkout. Everything the builder accumulated on the
# way (caches, tarballs, build tools, intermediate layers) stays behind.
# ---------------------------------------------------------------------------
FROM debian:bookworm-slim

# bash (entrypoint), ca-certificates (https), procps (`ps`, used by split-pipe's
# profiling loop), gawk (Parse's shell helpers), unzip (fastqc archives).
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bash ca-certificates procps gawk unzip && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /opt/parse /opt/parse

# Login shells rebuild PATH from /etc/profile (Debian's default), which would hide
# both conda envs: the continuumio base image ships its own conda profile hook that
# this slim base does not. Put the base env and the Parse tools first explicitly.
RUN printf '%s\n' 'export PATH="/opt/conda/bin:/opt/conda/envs/spipe/bin:$PATH"' \
        > /etc/profile.d/zz-spipe-path.sh && \
    chmod 0644 /etc/profile.d/zz-spipe-path.sh

# Python writes no .pyc: keeps the image read-only-friendly and avoids
# duplicating bytecode at run time.
ENV PYTHONDONTWRITEBYTECODE=true

# The base env (python 3.11, pipen, R) first, then the Parse tools.
ENV PATH=/opt/conda/bin:/opt/conda/envs/spipe/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

COPY pyproject.toml README.md /biopipen/
COPY biopipen/ /biopipen/biopipen/
COPY docker/ /biopipen/docker/
RUN chmod +x /biopipen/docker/splitpipe_pipeline/entry.sh

WORKDIR /workdir

ENTRYPOINT [ "/bin/bash", "/biopipen/docker/splitpipe_pipeline/entry.sh" ]
