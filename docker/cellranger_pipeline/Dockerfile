FROM biopipen/cellranger:9.0.1

# docker build /path/to/biopipen \
#   -f /path/to/biopipen/docker/cellranger_pipeline/Dockerfile
COPY . /biopipen

WORKDIR /biopipen
RUN apt-get update && \
    apt-get install -y python3.10 python3-pip npm cargo r-base && \
    Rscript -e 'install.packages(c("logger", "rlang", "jsonlite", "dplyr", "ggplot2", "ggprism"))' && \
    # alias python3 to python
    ln -s /usr/bin/python3 /usr/bin/python && \
    python -m pip config set global.break_system_packages on && \
    python -m pip install -U poetry setuptools_rust && \
    python -m poetry config virtualenvs.create false && \
    python -m poetry install -v -E runinfo && \
    pipen report update && \
    apt-get clean

# Add a user to avoid running as root
RUN useradd -m -s /bin/bash biopipen
USER biopipen

WORKDIR /example
RUN bash /biopipen/docker/cellranger_pipeline/make-examples.sh

WORKDIR /workdir

ENTRYPOINT [ "bash", "/biopipen/docker/cellranger_pipeline/entry.sh" ]
