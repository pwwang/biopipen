FROM biopipen/base:dev
FROM biopipen/cellranger:9.0.1

RUN RUN micromamba update -y -n base -f /biopipen/docker/cellranger_pipeline/env.yml && \
    micromamba clean --all --yes

WORKDIR /example
RUN bash /biopipen/docker/cellranger_pipeline/make-examples.sh

WORKDIR /workdir

ENTRYPOINT [ "/bin/bash", "/biopipen/docker/cellranger_pipeline/entry.sh" ]
