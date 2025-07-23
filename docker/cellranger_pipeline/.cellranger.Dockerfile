FROM ubuntu:rolling

ARG CELLRANGER_URL

WORKDIR /cellranger

RUN apt-get update && \
    apt-get install -y wget && \
    apt-get clean && \
    wget -O cellranger.tar.gz "$CELLRANGER_URL" && \
    tar -xzvf cellranger.tar.gz && \
    ln -s /cellranger/cellranger-*/bin/cellranger /usr/bin/cellranger && \
    rm -f /cellranger/cellranger-*.tar.gz

ENTRYPOINT ["cellranger"]
