# syntax=docker/dockerfile:1
FROM continuumio/miniconda3
RUN apt-get install libgl1-mesa-glx --yes
RUN conda create -n sphire eman-deps-dev=24 -c cryoem -c defaults -c conda-forge --yes
RUN conda init bash
RUN echo "conda activate sphire" > ~/.bashrc
