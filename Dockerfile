FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="LinkSeq-Demux" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

RUN apt update -yqq && \
    apt install -yqq \
    unzip

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/linkseq-demux/bin:$PATH
