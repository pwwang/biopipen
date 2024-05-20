FROM ubuntu:rolling

WORKDIR /cellranger
RUN apt-get update && \
    apt-get install -y wget && \
    apt-get clean && \
    wget -O cellranger-8.0.1.tar.gz "cellranger-URL" && \
    tar -xzvf cellranger-8.0.1.tar.gz && \
    ln -s /cellranger/cellranger-8.0.1/bin/cellranger /usr/bin/cellranger && \
    rm -f /cellranger/cellranger-8.0.1.tar.gz

ENTRYPOINT ["cellranger"]
