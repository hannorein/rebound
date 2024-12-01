#!/bin/zsh
export CC=gcc-mp-14 && make clean && make -j14 -B && pip install -e .
