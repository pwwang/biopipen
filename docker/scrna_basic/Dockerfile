FROM mambaorg/micromamba:1.4.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

# Install dependencies
RUN micromamba install -y -n base -f /biopipen/docker/scrna_basic/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v

WORKDIR /example
RUN bash /biopipen/docker/scrna_basic/make-examples.sh

WORKDIR /workdir
COPY docker/scrna_basic/board.toml /workdir/board.toml

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/biopipen/docker/scrna_basic/entry.sh" ]
