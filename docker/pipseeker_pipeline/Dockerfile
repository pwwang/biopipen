FROM public.ecr.aws/w3e1n2j6/fluent-pipseeker:3.3.0 AS pipseeker
FROM biopipen/base:latest

COPY --from=pipseeker /usr/local/bin/pipseeker /opt/conda/bin/pipseeker

# in case we have updates about the pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/pipseeker_pipeline /biopipen/docker/pipseeker_pipeline
RUN micromamba env update -n base -f /biopipen/docker/pipseeker_pipeline/env.yml && \
    micromamba clean --all --yes && \
    chmod +x /biopipen/docker/pipseeker_pipeline/entry.sh && \
    ln -s /biopipen/docker/pipseeker_pipeline/entry.sh /opt/conda/bin/entry.sh && \
    rm -rf /opt/conda/pkgs/*

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /example
RUN /bin/bash /biopipen/docker/pipseeker_pipeline/make-examples.sh && \
    cp /biopipen/docker/pipseeker_pipeline/board.toml /example/board.toml && \
    cp /biopipen/docker/pipseeker_pipeline/PipseekerPipeline.config.toml /example/PipseekerPipeline.config.toml

WORKDIR /workdir

ENTRYPOINT [ "/bin/bash", "/biopipen/docker/pipseeker_pipeline/entry.sh" ]