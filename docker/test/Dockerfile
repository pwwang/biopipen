FROM mambaorg/micromamba:2.3.0

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Install dependencies
RUN micromamba env install -y -n base -f /biopipen/tests/conda/env_biopipen.yml && \
    micromamba clean --all --yes && \
    poetry config virtualenvs.create false && \
    poetry install -v --all-extras && \
    pipen report update && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /home/$MAMBA_USER/.cache/*

WORKDIR /workdir
# Test it by:
# docker run -it --rm biopipen/test:dev \
#   pipen run vcf Vcf2Bed --in.infile /biopipen/tests/test_vcf/Vcf2Bed/data/tobed.vcf
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
