FROM andrewosh/binder-base

# for use with mybinder.org

MAINTAINER Hanno Rein  <hanno@hanno-rein.de>

USER root
COPY . $HOME/
RUN pip install -v -e .
RUN $HOME/anaconda2/envs/python3/bin/pip install -v -e .
