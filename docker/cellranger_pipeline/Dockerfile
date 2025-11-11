ARG REF_NAME=dev

FROM biopipen/cellranger:9.0.1 AS cellranger
FROM biopipen/base:${REF_NAME}

COPY --from=cellranger /cellranger /cellranger

# in case we have updates about the pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/cellranger_pipeline /biopipen/docker/cellranger_pipeline

RUN micromamba env update -n base -f /biopipen/docker/cellranger_pipeline/env.yml && \
    micromamba clean --all --yes && \
    ln -s /cellranger/cellranger-*/bin/cellranger /opt/conda/bin/cellranger && \
    chmod +x /biopipen/docker/cellranger_pipeline/entry.sh && \
    ln -s /biopipen/docker/cellranger_pipeline/entry.sh /opt/conda/bin/entry.sh && \
    /opt/conda/bin/cellranger telemetry disable all

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /example
RUN /bin/bash /biopipen/docker/cellranger_pipeline/make-examples.sh && \
    cp /biopipen/docker/cellranger_pipeline/board.toml /example/board.toml && \
    cp /biopipen/docker/cellranger_pipeline/CellrangerCountPipeline.config.toml /example/CellrangerCountPipeline.config.toml  && \
    cp /biopipen/docker/cellranger_pipeline/CellrangerVdjPipeline.config.toml /example/CellrangerVdjPipeline.config.toml

WORKDIR /workdir

ENTRYPOINT [ "/bin/bash", "/biopipen/docker/cellranger_pipeline/entry.sh" ]
