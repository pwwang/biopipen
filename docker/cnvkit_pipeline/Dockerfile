FROM mambaorg/micromamba:1.4.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

# Install dependencies
RUN micromamba install -y -n base -f /biopipen/docker/cnvkit_pipeline/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v && \
    mkdir /biopipen/.pipen-board && \
    pipen report update && \
    patch $(python -c 'from pipen_board import defaults; print(defaults.__file__)') \
        /biopipen/docker/cnvkit_pipeline/pipen-board.patch
    # Patch pipen-board to use a non-home directory to avoid failure in singularity

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/biopipen/docker/cnvkit_pipeline/entry.sh" ]
