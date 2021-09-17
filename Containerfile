FROM ubi8:latest

RUN yum update \
  && yum install -y python3-pip python3-devel pigz \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 --no-cache-dir install --upgrade pip \
  && yum clean all \
  && echo "system packages installed"

RUN python -m pip install pysam

WORKDIR /opt/
