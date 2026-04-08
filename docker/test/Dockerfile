FROM biopipen/test-deps:latest AS builder

COPY --chown=$MAMBA_USER:$MAMBA_USER . /biopipen

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /biopipen
# Install biopipen
RUN python -m pip install --no-cache-dir -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install --no-cache -v --all-extras && \
    python -m pip cache purge && \
    python -m poetry cache clear --all pypi && \
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
