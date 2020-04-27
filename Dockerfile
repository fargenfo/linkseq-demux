FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="LinkSeq-Demux -- basecall/demultiplex and trim linked-reads [WIP]" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

RUN apt_get update -yqq && \
    apt_get install -yqq \
    unzip

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/linkseq-demux/bin:$PATH
