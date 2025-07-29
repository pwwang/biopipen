ARG REF_NAME=dev

FROM biopipen/base:${REF_NAME}

# in case we have updates about the pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/cellranger_pipeline /biopipen/docker/cellranger_pipeline

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install dependencies
RUN micromamba install -y -n base -f /biopipen/docker/cnvkit_pipeline/env.yml && \
    micromamba clean --all --yes && \
    chmod +x /biopipen/docker/cnvkit_pipeline/entry.sh && \
    ln -s /biopipen/docker/cnvkit_pipeline/entry.sh /opt/conda/bin/entry.sh

WORKDIR /example
RUN /bin/bash /biopipen/docker/cnvkit_pipeline/make-examples.sh && \
    cp /biopipen/docker/cnvkit_pipeline/board.toml /example/board.toml && \
    cp /biopipen/docker/cnvkit_pipeline/CNVkitPipeline.config.toml /example/CNVkitPipeline.config.toml

WORKDIR /workdir

ENTRYPOINT [ "/bin/bash", "/biopipen/docker/cnvkit_pipeline/entry.sh" ]
