#!/bin/bash

export CUDA_VISIBLE_DEVICES=1

source /nfs/polizzi/shared/miniconda3/etc/profile.d/conda.sh &&
conda activate /nfs/polizzi/shared/miniconda3/envs/SE3nv &&
python rfaa.py
