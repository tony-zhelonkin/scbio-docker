FROM ubuntu:20.04 

# Environment
LABEL maintainers="Anton Zhelonkin (anton.bioinf.md@gmail.com)"
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
ARG DEBIAN_FRONTEND=noninteractive

# User definition 
# Add build arguments for UID/GID
ARG USER_ID=1001
ARG GROUP_ID=1001
ARG USER=mogilenko_lab
ARG GROUP=mogilenko_lab

# Create group and user matching host's UID/GID
RUN groupadd -g ${GROUP_ID} ${GROUP} && \
    useradd -m -u ${USER_ID} -g ${GROUP_ID} ${USER}

# Set user as owner of necessary directories
RUN mkdir -p /home/${USER}/reference /home/${USER}/data && \
    chown -R ${USER}:${GROUP} /home/${USER}
RUN echo "alias ll='ls -la -G'" >> /home/${USER}/.profile

# Previous user handling 
#RUN groupadd -g 2000 training && useradd -m -u 2000 -g 2000 training
#RUN echo 'training:training' | chpasswd
#RUN mkdir /home/training/reference
#RUN mkdir /home/training/data
#RUN chsh -s /bin/bash training
#RUN echo "alias ll='ls -la -G'" >> /home/training/.profile

#########
### Aptitude packages (updated openJDK from 8 to 17)
#########
RUN apt update && apt install -y --fix-missing \
    wget libghc-zlib-dev libfreetype6-dev libpng-dev libxft-dev \
    libncurses-dev git unzip ftp libzmq3-dev libboost-all-dev \
    vim nano \
    python3-dev python3-pip ftp \
    apache2 openssl libcurl4-openssl-dev curl \
    supervisor openssh-server libbz2-dev \
    liblzma-dev automake parallel cmake \ 
    openjdk-17-jre-headless shellinabox
    #usermod -G training,www-data training
# After creating the build directory
RUN mkdir build && \
    chown -R ${USER}:${GROUP} /build && \
    chmod 755 /build
WORKDIR /build


#########
### HTSlib tools: Samtools, Bcftools
#########
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 \
    https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 \
    https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -xf htslib-1.21.tar.bz2

WORKDIR /build/htslib-1.21
RUN ./configure && make -j 4 && make install
WORKDIR /build
RUN tar -xf bcftools-1.21.tar.bz2
WORKDIR /build/bcftools-1.21
RUN make -j 4 && make install
WORKDIR /build
RUN tar -xf samtools-1.21.tar.bz2
WORKDIR /build/samtools-1.21
RUN make -j 4 && make install
WORKDIR /build

#########
### STAR (2.7.11b update)
#########
#RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz && \
RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz && \
    tar -xf 2.7.11b.tar.gz
WORKDIR /build/STAR-2.7.11b/bin/Linux_x86_64_static
RUN cp STAR* /usr/local/bin
WORKDIR /build

#########
### kallisto
#########
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_linux-v0.51.1.tar.gz && \
    tar -xf kallisto_linux-v0.51.1.tar.gz
WORKDIR /build/kallisto
RUN cp kallisto /usr/local/bin
WORKDIR /build

#########
### Salmon
#########
WORKDIR /build
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar -xf salmon-1.10.0_linux_x86_64.tar.gz
WORKDIR /build/salmon-latest_linux_x86_64/
RUN cp bin/salmon /usr/local/bin && \
    cp lib/liblzma.so.5* /usr/local/lib && \
    cp lib/libtbb* /usr/local/lib
## building from source won`t work because of boost problems
#RUN wget https://github.com/COMBINE-lab/salmon/archive/refs/tags/v1.10.1.tar.gz  && \
#    tar -xf v1.10.1.tar.gz
#WORKDIR /build/salmon-1.10.1/
#RUN mkdir build
#WORKDIR /build/salmon-1.10.1/build
#RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/bin -DFETCH_BOOST=TRUE /build/salmon-1.10.1/ && \
#    make -j 4 && make install
#RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/bin /build/salmon-1.10.1/ && \
#    make -j 4 && make install
WORKDIR /build

#########
### BEDtools (updated version)
#########
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xf bedtools-2.31.1.tar.gz && cd bedtools2 && \
    make -j 4 && make install

#########
### Java tools: Trimmomatic, FastQC (check, update versions)
#########
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc && \
    cp -r FastQC /usr/share/ && \
    ln -s /usr/share/FastQC/fastqc /usr/bin/ && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    cp -r Trimmomatic-0.39 /usr/share/ && \
    echo '#!/bin/bash\njava -jar /usr/share/Trimmomatic-0.39/trimmomatic-0.39.jar $@' > /usr/bin/trimmomatic && \
    chmod +x /usr/bin/trimmomatic

#########
### Picard
########
# Clone the Picard repository
RUN git clone https://github.com/broadinstitute/picard.git
# Set the working directory to the cloned Picard directory
WORKDIR /build/picard
# Build Picard using Gradle and clean up after
RUN ./gradlew shadowJar
WORKDIR /build
RUN cp -r picard /usr/share/ && \
    echo '#!/bin/bash\njava -jar /usr/share/picard/build/libs/picard.jar $@' > /usr/bin/picard && \
    chmod +x /usr/bin/picard
# Return to the previous working directory
WORKDIR /build

#########
### Pip installs: HTSeq, PySam, MACS2, RSeQC, MultiQC, cutadapt
#########
RUN pip3 install numpy htseq PySam matplotlib scipy anndata loompy SWIG Cython==0.29.32 MACS3 RSeQC multiqc==1.25.2 cutadapt
WORKDIR /build

#########
### Trim Galore
#########
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz && \
    tar xvzf trim_galore.tar.gz
WORKDIR /build/TrimGalore-0.6.10
RUN cp trim_galore /usr/local/bin
WORKDIR /build

########
### quant3p
########
WORKDIR /build
RUN wget https://github.com/ctlab/quant3p/archive/refs/heads/master.zip && \ 
    unzip master.zip && \
    mv quant3p-master quant3p && \
    ln -s /usr/bin/python3 /usr/bin/python
WORKDIR /build/quant3p
RUN python setup.py install
WORKDIR /build

########
### featureCounts
#######
WORKDIR /build
RUN wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-Linux-x86_64.tar.gz && \
    tar -xvf subread-2.0.2-Linux-x86_64.tar.gz
WORKDIR /build/subread-2.0.2-Linux-x86_64/bin
RUN cp featureCounts /usr/bin/
WORKDIR /build
        
#########
### SortmeRNA
#########
#WORKDIR /build
#COPY sortmerna-3.0.4-Linux.sh /build/sortmerna-3.0.4-Linux.sh
#RUN sh sortmerna-3.0.4-Linux.sh --skip-license --prefix=/usr --exclude-subdir
#RUN mkdir /usr/share/rRNA_databases
#WORKDIR /usr/share/rRNA_databases
#RUN wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/rfam-5.8s-database-id98.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/rfam-5s-database-id98.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-arc-16s-id95.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-arc-23s-id98.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-euk-28s-id98.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-euk-18s-id95.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-bac-16s-id90.fasta && \
#    wget https://raw.githubusercontent.com/biocore/sortmerna/master/rRNA_databases/silva-bac-23s-id98.fasta
#COPY scripts/indexdb.sh /usr/share/rRNA_databases/
#COPY scripts/makeDbList.sh /usr/share/rRNA_databases/
#RUN ./indexdb.sh
#ENV SORTMERNA_DB="/usr/share/rRNA_databases/rfam-5s-database-id98.fasta,/usr/share/rRNA_databases/index/rfam-5s-database-id98:/usr/share/rRNA_databases/silva-euk-28s-id98.fasta,/usr/share/rRNA_databases/index/silva-euk-28s-id98:/usr/share/rRNA_databases/silva-euk-18s-id95.fasta,/usr/share/rRNA_databases/index/silva-euk-18s-id95:/usr/share/rRNA_databases/silva-arc-16s-id95.fasta,/usr/share/rRNA_databases/index/silva-arc-16s-id95:/usr/share/rRNA_databases/silva-bac-16s-id90.fasta,/usr/share/rRNA_databases/index/silva-bac-16s-id90:/usr/share/rRNA_databases/silva-arc-23s-id98.fasta,/usr/share/rRNA_databases/index/silva-arc-23s-id98:/usr/share/rRNA_databases/silva-bac-23s-id98.fasta,/usr/share/rRNA_databases/index/silva-bac-23s-id98:/usr/share/rRNA_databases/rfam-5.8s-database-id98.fasta,/usr/share/rRNA_databases/index/rfam-5.8s-database-id98"
#ADD scripts/merge-paired-reads.sh /usr/local/bin/
#ADD scripts/unmerge-paired-reads.sh /usr/local/bin/
#WORKDIR /build

#########
### BWA (updated version)
#########
RUN wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.18.tar.gz && \
    tar -xf v0.7.18.tar.gz
WORKDIR /build/bwa-0.7.18
RUN make -j 4 && cp bwa /usr/local/bin

######### CHIP-seq tools
### USeq and Sissrs
#########
RUN mkdir /data && \
    wget https://github.com/HuntsmanCancerInstitute/USeq/releases/download/USeq_9.3.5/USeq_9.3.5.zip && \
    unzip USeq_9.3.5.zip && mv USeq_9.3.5 /usr/share/USeq && \
    wget http://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/sissrs_v1.4.tar.gz && \
    tar -xf sissrs_v1.4.tar.gz && \
    cp sissrs.pl /usr/local/bin

############
### Install scripts
#ADD scripts/script.sh /usr/local/bin/


######### (updated version)
### JBrowse to easily view biological data formats
#########
#RUN rm /var/www/html/index.html && \
#    ln -sf /home/training /var/www/html/home
#WORKDIR /var/www/html
#RUN wget https://github.com/GMOD/jbrowse/releases/download/1.16.11-release/JBrowse-1.16.11.zip && \
#    unzip JBrowse-1.16.11.zip
#WORKDIR /var/www/html/JBrowse-1.16.11
#RUN ./setup.sh && \
#    mkdir -p data/bam && \
#    chown -R www-data:www-data /var/www/html && \
#    chmod a+rwx /var/www/html/JBrowse-1.16.11/data/bam && \
#    touch data/tracks.conf && \
#    chmod a+rw    /var/www/html/JBrowse-1.16.11/data/tracks.conf
#WORKDIR /
#RUN rm /var/www/html/JBrowse-1.16.11.zip
#ADD scripts/add_JBrowse_tracks.sh /usr/local/bin/add_JBrowse_tracks.sh

#########
### R installs
#########
#ENV R_LIBS="/home/training/.r-library"
#ADD scripts/R_installs.R /build/R_installs.R
#ADD .Renviron /home/training/.Renviron
#RUN mkdir /home/training/.r-library && \
#    chown -R training:training /home/training/.Renviron && \
#    chown -R training:training /home/training/.r-library && \
#    chmod 755 /home/training/.r-library && \
#    sudo -u training Rscript /build/R_installs.R

#ADD supervisord.conf /etc/supervisor/conf.d/supervisord.conf
#VOLUME /home/training/share
#EXPOSE 80 8888 8787 443
#CMD ["/usr/bin/supervisord","-c","/etc/supervisor/conf.d/supervisord.conf"]

#########
### Cleaning up the build 
#########
WORKDIR /build
RUN rm -rf /build/*

# Switch to non-root user for remaining operations
USER ${USER}
RUN mkdir -p /home/${USER}/analysis \
    /home/${USER}/data \
    /home/${USER}/scripts
#RUN rm *tar* && \
#    rm *zip
