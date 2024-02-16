#!/bin/bash

yes | conda create --name wlcstat python=3.9.12
eval "$(conda shell.bash hook)"
conda activate wlcstat
yes | pip install -r requirements.txt
pip install -e .