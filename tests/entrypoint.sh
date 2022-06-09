#!/bin/bash

cd /workdir

pip install -U pip
pip install -U poetry
poetry config virtualenvs.create false

poetry install

python $1
