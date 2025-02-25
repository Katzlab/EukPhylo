FROM python:3.8.20-bullseye
RUN apt-get update -y && apt-get upgrade -y && \
    apt-get install git wget ncbi-blast+-legacy mafft build-essential vim -y

# PhyloToL installation
# Tested on: RUN git clone -b patch-1 https://github.com/DLBPointon/EukPhylo.git
# Will work on below container after PR merge

RUN git clone -b Docker https://github.com/Katzlab/EukPhylo.git

# IQ-TREE installation
RUN wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/iqtree-1.6.12-Linux.tar.gz && \
    tar xzf iqtree-1.6.12-Linux.tar.gz && \
    mv iqtree-1.6.12-Linux/ /usr/local/iqtree && \
    ln -s /usr/local/iqtree/bin/iqtree /usr/local/bin/iqtree

# Install python dependencies
RUN pip install biopython==1.75 ete3==3.1.2 tqdm==4.66.4 six==1.17.0

# VSEARCH installation
RUN wget https://github.com/torognes/vsearch/releases/download/v2.21.0/vsearch-2.21.0-linux-aarch64.tar.gz && \
    tar xzf vsearch-2.21.0-linux-aarch64.tar.gz && \
    mv vsearch-2.21.0-linux-aarch64 /usr/local/vsearch && \
    ln -s /usr/local/vsearch/bin/vsearch /usr/local/bin/vsearch

# Guidance installation
RUN mkdir guidance && \
    wget https://taux.evolseq.net/guidance/static/download/guidance.v2.02.tar.gz && \
    tar -xzvf guidance.v2.02.tar.gz -C guidance --no-same-owner && \
    cd guidance/guidance.v2.02 && make

# Install pre-compiled version of Trimal 1.5.0
RUN wget https://github.com/inab/trimal/releases/download/v1.5.0/trimAl_Linux_x86-64.zip && \
    unzip trimAl_Linux_x86-64.zip && \
    mv trimAl_Linux_x86-64 /usr/local/trimal

# Clean the container
RUN rm -rf iqtree-1.6.12-Linux.tar.gz vsearch-2.21.0-linux-aarch64.tar.gz guidance.v2.02.tar.gz trimAl_Linux_x86-64.zip

# Add executables to path
ENV PATH="$PATH:/iqtree/bin/iqtree:/usr/local/vsearch/bin/vsearch:/usr/local/trimal/trimal:/usr/local/bin"

# Change executable permissions on wrapper script
RUN ["chmod", "+x", "/EukPhylo/PTL2/run_eukphylo.sh"]

# Change entrypoint on execution of the container
ENTRYPOINT ["bash","/EukPhylo/PTL2/run_eukphylo.sh"]
