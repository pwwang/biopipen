#!/bin/bash

cd /workdir

pip install -U pip
pip install -U poetry
# poetry config virtualenvs.create false
poetry export -f requirements.txt --output requirements.txt

pip install -r requirements.txt

python $1
