# specific version for reproducibilty
FROM phusion/passenger-customizable:0.9.19

LABEL maintainer Austin Zadoks <austin.zadoks@epfl.ch>

ENV HOME /root

############
### init ###
############
CMD ["/sbin/my_init"]

########################
### passeneger setup ###
########################
# buildessentials and git
RUN /pd_build/utilities.sh \
    # python 2.7 and 3.5.x
    && /pd_build/python.sh \
    # enable Nginx and passenger
    && rm -f /etc/service/ngnix/down

########################
### Root environment ###
########################
RUN cp /etc/skel/.bash* /etc/skel/.profile /root/
RUN echo "PS1='docker:\W$ '" >> /root/.bashrc

################
### programs ###
################
RUN apt-get update && \
    apt-get -y install \
    python-pip \
    gfortran \
    liblapack-dev \
    libblas-dev \
    libnetcdf-dev \
    netcdf-bin \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean all

##############
### python ###
##############
RUN pip install --upgrade pip
RUN pip install --no-cache-dir \
    numpy \
    scipy \
    pandas \
    pymatgen \
    ase \
    spglib \
    flask_restful \
    tqdm \
    psutil

#####################
### QUIP + quippy ###
#####################
ENV QUIP_ROOT /opt/QUIP_ROOT
# get QUIP
RUN git clone https://github.com/libAtoms/QUIP.git ${QUIP_ROOT}
WORKDIR ${QUIP_ROOT}
# for reproducibility 05-07-2018
RUN git checkout 531330f
# don't install GAP
ENV BUILD NOGAP
# QUIP for general use is the OpenMP version.
# Installs with no suffix, e.g. quip
# Also installs quippy
ENV QUIP_ARCH linux_x86_64_gfortran_openmp
RUN cd ${QUIP_ROOT} \
    && mkdir -p build/${QUIP_ARCH} \
    && cp docker/arch/${BUILD}_Makefile.${QUIP_ARCH}.inc build/${QUIP_ARCH}/Makefile.inc \
    && make > /dev/null \
    && QUIP_INSTALLDIR=${QUIP_ROOT}/bin make install \
    && make install-quippy > /dev/null \
    && find build/${QUIP_ARCH} -type f ! \( -name '*.a' -o -name 'Makefile.inc' \) -delete

################
### app user ###
################
USER app
WORKDIR /home/app
ENV HOME /home/app

#######################
### SOAP REST Setup ###
#######################
RUN git clone https://github.com/cosmo-epfl/glosim2
RUN git clone https://github.com/zooks97/aiida-fingerprint
RUN mkdir webservice
RUN cp aiida-fingerprint/aiida_fingerprint/soaprest.py webservice/soaprest.py

#############
### Nginx ###
#############
RUN mkdir webservice/public
ENV SP_WSGI_FILE=webservice/soaprest_nginx.py
RUN echo "import sys" > $SP_WSGI_FILE \
    && echo "sys.path.insert(0, '/home/app/webservice')" >> $SP_WSGI_FILE \ 
    && echo "from soaprest import app as application" >> $SP_WSGI_FILE 

#################
### root user ###
#################
USER root
WORKDIR /root

#############
### Nginx ###
#############
RUN rm /etc/nginx/sites-enabled/default
ADD conf/soaprest-nginx.conf /etc/nginx/sites-enabled/soaprest.conf

############
### init ###
############
RUN mkdir -p /etc/my_init.d

###############
### startup ###
###############
# webserver
EXPOSE 80