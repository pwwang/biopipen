FROM nfcore/cellranger:7.2.0

COPY . /biopipen

WORKDIR /biopipen
RUN python -m pip install -U poetry && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v -E runinfo && \
    pipen report update

WORKDIR /example
RUN bash /biopipen/docker/cellranger_pipeline/make-examples.sh

WORKDIR /workdir

ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh", "bash", "/biopipen/docker/cellranger_pipeline/entry.sh" ]
