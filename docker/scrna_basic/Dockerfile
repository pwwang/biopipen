FROM mambaorg/micromamba:1.4.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

# Install dependencies
RUN micromamba install -y -n base -f /biopipen/docker/scrna_basic/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v && \
    mkdir /biopipen/.pipen-board && \
    patch $(python -c 'from pipen_board import defaults; print(defaults.__file__)') /biopipen/docker/scrna_basic/pipen-board.patch
    # Patch pipen-board to use a non-home directory to avoid failure in singularity

WORKDIR /example
RUN bash /biopipen/docker/scrna_basic/make-examples.sh

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/biopipen/docker/scrna_basic/entry.sh" ]
