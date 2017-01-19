FROM andrewosh/binder-base

# for use with mybinder.org

MAINTAINER Hanno Rein  <hanno@hanno-rein.de>

USER root
COPY . $HOME/
RUN $HOME/anaconda2/envs/python3/bin/conda upgrade -y notebook 
RUN $HOME/anaconda2/envs/python3/bin/conda install -y -c conda-forge ipywidgets  
RUN conda upgrade -y notebook 
RUN conda install -y -c conda-forge ipywidgets  
RUN pip install -v -e .
#RUN pip install ipywidgets
RUN $HOME/anaconda2/envs/python3/bin/pip install -v -e .
#RUN $HOME/anaconda2/envs/python3/bin/pip install ipywidgets
#RUN jupyter nbextension enable --py --sys-prefix widgetsnbextension
