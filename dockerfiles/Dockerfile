FROM continuumio/miniconda3:latest

####################
##### megadock #####
####################

RUN apt-get update -y
RUN apt-get install ffmpeg libsm6 libxext6  -y
RUN apt-get install -y --no-install-recommends \
                wget \
                make \
                g++ \
    && rm -rf /var/lib/apt/lists/*


# FFTW
ENV FFTW_TARGET fftw-3.3.8

RUN wget -P /tmp http://www.fftw.org/${FFTW_TARGET}.tar.gz && \
    tar xzf /tmp/${FFTW_TARGET}.tar.gz -C /tmp && \
    rm -rf /tmp/${FFTW_TARGET}.tar.gz && \
    cd /tmp/${FFTW_TARGET} && \
    ./configure --enable-float --enable-sse2 --prefix=/usr/local/${FFTW_TARGET} && \
    make -j$(nproc) && \
    make install


# MEGADOCK
ENV MEGADOCK_WORK_DIR /opt/MEGADOCK

RUN git clone https://github.com/akiyamalab/MEGADOCK.git /opt/MEGADOCK
RUN cp /opt/MEGADOCK/Dockerfiles/cpu/Makefile ${MEGADOCK_WORK_DIR}

RUN cd ${MEGADOCK_WORK_DIR} && \
    make -j$(nproc)

WORKDIR ${MEGADOCK_WORK_DIR}

RUN ln -s ${MEGADOCK_WORK_DIR}/megadock /usr/local/bin/megadock && \
    ln -s ${MEGADOCK_WORK_DIR}/decoygen /usr/local/bin/decoygen && \
    ln -s ${MEGADOCK_WORK_DIR}/block /usr/local/bin/block && \
    ln -s ${MEGADOCK_WORK_DIR}/ppiscore /usr/local/bin/ppiscore


###################
###### conda ######
###################

COPY ./dockerfiles/conda_requirements.txt /home/conda_requirements.txt
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda create -n p4ward python=3.11 --file /home/conda_requirements.txt

###################
##### protacs #####
###################

COPY ./p4ward /home/python/p4ward
ENV PYTHONPATH /home/python:$PYTHONPATH

SHELL ["conda", "run", "-n", "p4ward", "/bin/bash", "-c"]

RUN mkdir /home/data
WORKDIR /home/data

ENTRYPOINT [ "conda", "run", "--no-capture-output", "-n","p4ward", "python", "-m", "p4ward" ]
