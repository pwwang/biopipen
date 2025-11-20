FROM biopipen/base:latest AS builder

# in case we have update dependencies
COPY tests/conda/env_biopipen.yml /biopipen/tests/conda/env_biopipen.yml

# Install dependencies
RUN micromamba env update -n base -f /biopipen/tests/conda/env_biopipen.yml && \
    micromamba clean --all --yes && \
    rm -rf /opt/conda/pkgs/*

FROM biopipen/base:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER --from=builder /opt/conda /opt/conda

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /workdir
# Test it by:
# docker run -it --rm biopipen/test:dev \
#   pipen run vcf Vcf2Bed --in.infile /biopipen/tests/test_vcf/Vcf2Bed/data/tobed.vcf
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
