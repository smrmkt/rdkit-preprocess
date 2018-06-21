FROM continuumio/anaconda3:latest
MAINTAINER smrmkt <smkrmkt@gmail.com>

RUN conda install -y -c rdkit rdkit && pip install boto3 && mkdir -p /data

COPY mp_calculate_descriptors.py /data/

