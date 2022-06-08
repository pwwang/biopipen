#!/bin/bash

cd /workdir

pip install -U pip
pip install -U poetry
poetry config virtualenvs.create false

# Use wheel to skip building from rust, which takes too long
pip install -U rtoml
poetry install

python $1
