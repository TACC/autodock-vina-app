FROM tacc/tacc-ubuntu18-mvapich2.3-ib:latest

# Install system updates
RUN apt-get update && apt-get upgrade -y

# Install pip, pip3, and use them to install others
RUN apt-get install -y python3-pip \
    && apt-get install wget -y \ 
    && curl -O https://bootstrap.pypa.io/pip/3.6/get-pip.py \
    && python3 get-pip.py \
    && pip3 install mpi4py==3.1.4 numpy==1.19.5 blosc==1.10.6 \
    && pip install -U numpy vina

# Install ADFRSuite
RUN wget -O ADFR.tar  https://ccsb.scripps.edu/adfr/download/1038/ \
    && tar -xf ADFR.tar \
    && cd ADFRsuite_x86_64Linux_1.0 \ 
    && ./install.sh -c 0 \
    && cd 

# Clean up
RUN rm -f /ADFR.tar /get-pip.py

# Set path to include ADFR scripts
ENV PATH="$PATH:/ADFRsuite_x86_64Linux_1.0/bin"

# Copy scripts into container
COPY ./autodock-src /autodock-src
ENV PATH="$PATH:/autodock-src"
RUN chmod +x /autodock-src/autodock.py

# COPY assets/runner.sh /tapis/assets/runner.sh
# RUN chmod +x /tapis/assets/runner.sh

ENTRYPOINT ["python3", "/autodock-src/autodock.py"]
