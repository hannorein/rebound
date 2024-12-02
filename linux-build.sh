#!/bin/bash
make clean
make -j24 -B
pip install -e .
