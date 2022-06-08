#!/bin/bash

cd /biopipen

pip install -U poetry
poetry config virtualenvs.create false

poetry install

python $1
