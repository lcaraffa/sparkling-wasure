FROM mambaorg/micromamba:latest as mamba_pdal
COPY environment.yml /environment.yml
RUN micromamba env create -n pdaltools -f /environment.yml


FROM debian:bullseye-slim

# install PDAL
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/bin/pdal /opt/conda/envs/pdaltools/bin/pdal
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/bin/python /opt/conda/envs/pdaltools/bin/python
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/bin/las2las /opt/conda/envs/pdaltools/bin/las2las
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/bin/lasinfo /opt/conda/envs/pdaltools/bin/lasinfo
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/lib/ /opt/conda/envs/pdaltools/lib/
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/ssl /opt/conda/envs/pdaltools/ssl
COPY --from=mamba_pdal /opt/conda/envs/pdaltools/share/proj/proj.db /opt/conda/envs/pdaltools/share/proj/proj.db

ENV PATH=$PATH:/opt/conda/envs/pdaltools/bin/
ENV PROJ_LIB=/opt/conda/envs/pdaltools/share/proj/

WORKDIR /pdaltools
RUN mkdir tmp
COPY pdaltools pdaltools
COPY test test
