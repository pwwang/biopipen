FROM mambaorg/micromamba:2.3.0

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

# Install dependencies
RUN micromamba install -y -n base -f /biopipen/tests/conda/env_biopipen.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v && \
    pipen report update

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
