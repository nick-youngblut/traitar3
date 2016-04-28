############################################################
# Dockerfile to run Traitar - the microbial trait analyzer 
# Based on Ubuntu Image
############################################################

# Set the base image to use to Ubuntu
FROM ubuntu:trusty 

RUN apt-get update
MAINTAINER Aaron Weimann (weimann@hhu.de)
RUN apt-get install -y python-scipy python-matplotlib python-pip python-pandas
RUN echo "deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse ">> /etc/apt/sources.list 
RUN apt-get update
RUN apt-get install -y hmmer prodigal
RUN apt-get install -y wget 
RUN apt-get install -y git
RUN mkdir /home/traitar
WORKDIR  /home/traitar
RUN git clone https://github.com/aweimann/traitar
WORKDIR  /home/traitar/traitar
RUN python setup.py sdist
RUN pip install traitar  --find-links file:///home/traitar/traitar/dist
RUN wget -O /home/traitar/Pfam-A.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz
RUN gunzip /home/traitar/Pfam-A.hmm.gz 
RUN traitar pfam --local /home/traitar
