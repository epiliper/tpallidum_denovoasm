FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
        python3 python3-pip \
        git \
        wget \
        ncbi-blast+ \
        miniasm \
        racon \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install SPAdes 4.1.0 from prebuilt binary
RUN wget -q https://github.com/ablab/spades/releases/download/v4.1.0/SPAdes-4.1.0-Linux.tar.gz \
    && tar -xzf SPAdes-4.1.0-Linux.tar.gz -C /opt \
    && ln -s /opt/SPAdes-4.1.0-Linux/bin/spades.py /usr/local/bin/spades.py \
    && rm SPAdes-4.1.0-Linux.tar.gz

# Install Unicycler from GitHub
RUN pip3 install git+https://github.com/rrwick/Unicycler.git@v0.5.1 psutil

RUN unicycler --version \
    && spades.py --version \
    && python3 -c "import psutil; print('psutil', psutil.__version__)"

CMD ["unicycler", "--help"]
