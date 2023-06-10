FROM mambaorg/micromamba:1.4.3

COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp

# Install dependencies
RUN micromamba install -y -n base -f /tmp/docker/scrna_basic/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /tmp
RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v && \
    bash /tmp/docker/scrna_basic/make-examples.sh

WORKDIR /workdir
COPY docker/scrna_basic/board.toml /workdir/board.toml

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/tmp/docker/scrna_basic/entry.sh" ]
