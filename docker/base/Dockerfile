FROM mambaorg/micromamba:2.3.0 AS builder

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

WORKDIR /biopipen

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Install dependencies
# nodejs is for report building (pipen-report)
# procps-ng is for `free` command for memory info by pipen-runinfo
RUN micromamba install -y -n base python=3.10 nodejs=20 git poetry procps-ng && \
    micromamba clean --all --yes && \
    poetry config virtualenvs.create false && \
    poetry install -v && \
    pipen report update

FROM mambaorg/micromamba:2.3.0

COPY --chown=$MAMBA_USER:$MAMBA_USER --from=builder /biopipen /biopipen
COPY --chown=$MAMBA_USER:$MAMBA_USER --from=builder /opt/conda /opt/conda

WORKDIR /workdir
# to test it:
# docker run -it --rm biopipen/base:dev pipen run bed --help
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
