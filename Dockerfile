FROM continuumio/miniconda3
# RUN apt-get update && apt-get install -y \
#     vim \
#     git
# WORKDIR /
# RUN pip install --extra-index-url https://testpypi.python.org/pypi vcf-annotation-tools
RUN pip install vcfpy pysam
COPY somatic_llr_filter.py /usr/bin/somatic_llr_filter.py
