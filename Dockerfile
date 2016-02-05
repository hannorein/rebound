FROM andrewosh/binder-base

# for use with mybinder.org

MAINTAINER Hanno Rein  <hanno@hanno-rein.de>

USER root
COPY . $HOME/
RUN ls *
RUN pip install -v -e .
