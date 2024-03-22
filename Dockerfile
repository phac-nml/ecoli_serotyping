FROM ubuntu:22.04
ENV DEBIAN_FRONTEND="noninteractive" TZ="America/New_York"
RUN apt update && apt install git python3-pip -y
RUN apt install libcurl4-openssl-dev libssl-dev -y
RUN pip3 install Cython numpy
RUN apt install mash ncbi-blast+  bowtie2 seqtk samtools bcftools -y
RUN git clone https://github.com/phac-nml/ecoli_serotyping.git
# install the tool and initilize its species ID MASH database
RUN cd ecoli_serotyping && git checkout v2.0.0 && pip3 install .
RUN ectyper_init

#build image:  docker build --tag  ectyper:2.0.0 .
#type a sample: docker run -it --rm -v $PWD:/mnt ectyper:2.0.0 ectyper -i /mnt/assembly.fasta -o temp/ --pathotype