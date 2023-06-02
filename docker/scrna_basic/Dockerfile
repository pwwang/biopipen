FROM mambaorg/micromamba:1.4.3

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/scrna_basic/env.yaml /tmp/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/scrna_basic/entry.sh /tmp/entry.sh

# Install dependencies
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v

WORKDIR /workdir
COPY board.toml /workdir/board.toml

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/tmp/entry.sh" ]
