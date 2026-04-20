FROM biopipen/test-deps:latest AS builder

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
# Install biopipen
RUN uv pip install --system -e .[runinfo,log2file] && \
    uv cache clean && \
    pipen report update && \
    python /biopipen/docker/cleanup.py

FROM mambaorg/micromamba:2.5.0

COPY --chown=$MAMBA_USER:$MAMBA_USER --from=builder /biopipen /biopipen
COPY --chown=$MAMBA_USER:$MAMBA_USER --from=builder /opt/conda /opt/conda

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN fc-cache -f -v && \
    echo "cache=/tmp/npm-cache" > /home/$MAMBA_USER/.npmrc

WORKDIR /workdir
# Test it by:
# docker run -it --rm biopipen/test:dev \
#   pipen run vcf Vcf2Bed --in.infile /biopipen/tests/test_vcf/Vcf2Bed/data/tobed.vcf
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh" ]
