FROM ubuntu:rolling

COPY . /biopipen

WORKDIR /biopipen
RUN apt-get update && \
    apt-get install -y python3 python3-pip npm cargo wget && \
    python3 -m pip config set global.break_system_packages on && \
    python3 -m pip install -U poetry setuptools_rust && \
    python3 -m poetry config virtualenvs.create false && \
    python3 -m poetry install -v -E runinfo && \
    pipen report update && \
    apt-get remove -y cargo

# Add a user to avoid running as root
RUN useradd -m -s /bin/bash biopipen
USER biopipen

WORKDIR /example
RUN bash /biopipen/docker/cellranger_pipeline/make-examples.sh

WORKDIR /workdir

ENTRYPOINT [ "bash", "/biopipen/docker/cellranger_pipeline/entry.sh" ]
