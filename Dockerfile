FROM ubuntu:18.04

WORKDIR /opt
COPY . /opt

USER root
ARG DEBIAN_FRONTEND=noninteractive
ARG PYTHON_VERSION=3.7.1

# Update and install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    libspatialindex-dev \
    build-essential \
    software-properties-common \
    apt-utils \
    libgl1-mesa-glx \
    ffmpeg \
    libsm6 \
    libxext6 \
    libffi-dev \
    libbz2-dev \
    zlib1g-dev \
    libreadline-gplv2-dev \
    libncursesw5-dev \
    libssl-dev \
    libsqlite3-dev \
    sqlite3 \
    tk-dev \
    libgdbm-dev \
    libc6-dev \
    liblzma-dev \
    libxrender-dev \
    libgl1-mesa-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install PROJ 6.3.1
RUN wget https://download.osgeo.org/proj/proj-6.3.1.tar.gz && \
    tar -xvzf proj-6.3.1.tar.gz && \
    cd proj-6.3.1 && \
    ./configure && make -j$(nproc) && make install && ldconfig && \
    cd .. && rm -rf proj-6.3.1 proj-6.3.1.tar.gz

# Install PROJ datumgrid data files
RUN wget https://download.osgeo.org/proj/proj-datumgrid-1.8.tar.gz && \
    tar -xvzf proj-datumgrid-1.8.tar.gz && \
    mkdir -p /usr/local/share/proj && \
    find proj-datumgrid-1.8 -type f \( -name "*.gsb" -o -name "*.gtx" -o -name "*.dat" \) | xargs -I {} cp {} /usr/local/share/proj/ && \
    rm -rf proj-datumgrid-1.8 proj-datumgrid-1.8.tar.gz

# Build Python from source
RUN cd /opt && \
    wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz && \
    tar xzf Python-${PYTHON_VERSION}.tgz && \
    cd Python-${PYTHON_VERSION} && \
    ./configure --with-ensurepip=install && \
    make install && \
    cd .. && rm -rf Python-${PYTHON_VERSION} Python-${PYTHON_VERSION}.tgz

# Update PATH so custom Python is default
ENV PATH="/usr/local/bin:$PATH"

# Ensure pip is installed and upgraded
RUN /usr/local/bin/python3.7 -m ensurepip && \
    /usr/local/bin/python3.7 -m pip install --upgrade pip

# Install GDAL 3.0.4 from source
RUN wget http://download.osgeo.org/gdal/3.0.4/gdal-3.0.4.tar.gz && \
    tar -xvzf gdal-3.0.4.tar.gz && \
    cd gdal-3.0.4 && \
    ./configure && make -j$(nproc) && make install && ldconfig && \
    cd .. && rm -rf gdal-3.0.4 gdal-3.0.4.tar.gz

# Set environment variables for GDAL headers
ENV CPLUS_INCLUDE_PATH=/usr/local/include
ENV C_INCLUDE_PATH=/usr/local/include

# Install Python packages using correct interpreter
RUN /usr/local/bin/python3.7 -m pip install --no-cache-dir -r /opt/requirements.txt

# Install iRODS
RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - && \
    echo "deb [arch=amd64] https://packages.irods.org/apt/ $(lsb_release -sc) main" | tee /etc/apt/sources.list.d/renci-irods.list && \
    apt-get update -y && apt-get upgrade -y && \
    wget -c http://security.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2.24_amd64.deb && \
    apt-get install -y ./libssl1.1_1.1.1f-1ubuntu2.24_amd64.deb && \
    rm -rf ./libssl1.1_1.1.1f-1ubuntu2.24_amd64.deb && \
    apt install -y irods-icommands && \
    mkdir -p /root/.irods && \
    echo "{ \"irods_zone_name\": \"iplant\", \"irods_host\": \"data.cyverse.org\", \"irods_port\": 1247, \"irods_user_name\": \"$IRODS_USER\" }" > /root/.irods/irods_environment.json && \
    apt-get autoremove -y && apt-get clean

# Dynamically set PYTHONPATH
RUN SITE_PACKAGES=$(/usr/local/bin/python3.7 -c "import site; print(site.getsitepackages()[0])") && \
    echo "PYTHONPATH=$SITE_PACKAGES:\$PYTHONPATH" >> /etc/environment

ENV PYTHONPATH=/usr/local/lib/python3.7/site-packages:/usr/lib/python3/dist-packages:$PYTHONPATH

ENTRYPOINT ["/usr/local/bin/python3.7", "/opt/environmental_association.py"]