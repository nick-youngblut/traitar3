############################################################
# Dockerfile to run traitar - the microbial trait analyzer 
# Based on Ubuntu Image
############################################################

# Set the base image to use to Ubuntu
FROM ubuntu:14.04

RUN apt-get update
MAINTAINER Aaron Weimann (weimann@hhu.de)
RUN apt-get install -y python-scipy python-matplotlib python-pip python-pandas
RUN echo "deb http://archive.ubuntu.com/ubuntu trusty-backports main restricted universe multiverse ">> /etc/apt/sources.list 
RUN sudo apt-get update
RUN apt-get install -y hmmer prodigal
RUN apt-get install -y wget 
COPY dist/traitar-1.0.0.tar.gz /tmp
COPY traitar/data/sample_data /tmp/sample_data
RUN pip install /tmp/traitar-1.0.0.tar.gz
COPY Pfam-A.hmm /tmp
RUN traitar pfam --local /tmp
