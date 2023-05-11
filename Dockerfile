# latch base image + dependencies for latch SDK --- removing these will break the workflow
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9c8f-main
run pip install latch==2.19.7
run mkdir /opt/latch

run apt-get update -y &&\
    apt-get install -y autoconf samtools bwa bcftools

# copy all code from package (use .dockerignore to skip files)
copy . /root/

# latch internal tagging system + expected root directory --- changing these lines will break the workflow
arg tag
env FLYTE_INTERNAL_IMAGE $tag
workdir /root
