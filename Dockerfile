# use multi-stage build to install snekmer and avoid package dependency issues
# environment.yml needs to be updated manually if the snekmer version is updated here
FROM condaforge/mambaforge:latest AS mambasetup
COPY ./environment.yml .
RUN mamba env create -f ./environment.yml

# combine mambasetup build into kbase image
FROM kbase/sdkbase2:python
MAINTAINER KBase Developer

# share packages from the mambasetup stage build
# the path must also be updated in the entrypoint.sh
COPY --from=mambasetup /opt/conda/envs/snekmer/. /opt/conda/envs/snekmer/
ENV PATH /opt/conda/envs/snekmer/bin:$PATH

# kbase sdk code
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module
WORKDIR /kb/module
RUN make all
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]
