# Specifying the base image
FROM condaforge/mambaforge:23.3.1-1

RUN mamba create -n simbal

# May be optional
RUN mamba init bash

# Equivalent command to mamba activate simbal
SHELL ["mamba", "run", "-n", "simbal", "/bin/bash", "-c"]

RUN mamba install -y  -c conda-forge python==3.9.15

# Specifying the install packages in specific order instead of .yml file
RUN mamba install -y -c openalea3 -c conda-forge openalea.plantgl openalea.mtg

RUN mamba install -y -c openalea3 -c conda-forge openalea.caribu

#RUN mamba install -y -c conda-forge xarray==2023.3.0 dask==2023.3.2 bottleneck==1.3.7

#RUN python -m pip install netcdf4==1.6.3

#RUN python -m pip install --force-reinstall charset-normalizer==3.1.0

RUN python -m pip install pandas==1.5.3

RUN python -m pip install matplotlib==3.7.0

RUN python -m pip install scikit-learn==1.2.2

RUN python -m pip install numpy==1.22.4

#RUN mamba install -y -c conda-forge umap-learn==0.5.3

#RUN mamba install -y -c conda-forge hdbscan==0.8.29

RUN python -m pip install tensorflow==2.12.0

#RUN python -m pip install pyncclient

#RUN apt-get update && apt-get -y install libgl1

#python project
RUN mkdir pp 

ADD . ./pp/simbal

WORKDIR ./pp/simbal

#RUN python setup.py develop

VOLUME ./pp/simbal/src/simbal

# Mandatory for good packages indexing
USER root

CMD python src/simbal/main.py
