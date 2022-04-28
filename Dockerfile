# create multi-stage build starting with mambaforge image that has all needed snekmer packages
# might want to upgrade pip (if not in environment.yml)
FROM condaforge/mambaforge:latest AS mambasetup
RUN git clone --branch v0.1.2-beta https://github.com/PNNL-CompBio/Snekmer.git && \
    mamba env update -n base -f ./Snekmer/environment.yml && \
    pip install opencv-python  && \
    apt update && apt install -y libsm6 libxext6 && \
    apt-get install -y libxrender-dev && \
    pip install Snekmer/.

# combine mambaforge build into kbase image
# later might be able to use kbase/sdkbase2:latest since its smaller than :python,
# but I initialized module as python so the tests don't work with :latest
FROM kbase/sdkbase2:python
MAINTAINER KBase Developer

WORKDIR /kb/module
# mamba env section
COPY --from=mambasetup /opt/conda/. /opt/conda/
ENV PATH /opt/conda/bin:$PATH

# kbase sdk code for wrapper
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module
WORKDIR /kb/module
RUN make all
ENTRYPOINT [ "./scripts/entrypoint.sh" ]
CMD [ ]



