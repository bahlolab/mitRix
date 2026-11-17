FROM mambaorg/micromamba:1.5.8

USER root

COPY docker/env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes

WORKDIR /app

COPY mitRix .
COPY refs/chrMC_rCRS.fa     refs/
COPY refs/chrMC_rCRS.fa.amb refs/
COPY refs/chrMC_rCRS.fa.ann refs/
COPY refs/chrMC_rCRS.fa.bwt refs/
COPY refs/chrMC_rCRS.fa.pac refs/
COPY refs/chrMC_rCRS.fa.sa  refs/

RUN chmod +x /app/mitRix
