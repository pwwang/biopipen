FROM mambaorg/micromamba:2.3.0

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

# Install dependencies
RUN micromamba install -y -n base python=3.10 nodejs=20 git poetry && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
RUN poetry config virtualenvs.create false && \
    poetry install -v && \
    pipen report update

WORKDIR /workdir

# to test it:
# docker run -it --rm biopipen/base:dev pipen run bed --help
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
