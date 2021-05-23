FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential gfortran libblas-dev liblapack-dev liblapacke-dev \
    git python3-dev mpich python3-pip libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone CellConstructor
WORKDIR /home

RUN pip install mpi4py 
RUN git clone https://github.com/SSCHAcode/CellConstructor.git
RUN pip install -r CellConstructor/requirements.txt
RUN cd CellConstructor && python3 setup.py install && cd ..
RUN git clone https://github.com/SSCHAcode/python-sscha.git
RUN cd python-sscha && python3 setup.py install && cd ..

ENTRYPOINT ["/bin/bash"]

