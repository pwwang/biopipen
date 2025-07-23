FROM biopipen/base:dev

# Install dependencies
RUN micromamba install -y -n base -f /biopipen/tests/conda/env_biopipen.yml && \
    micromamba clean --all --yes

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
