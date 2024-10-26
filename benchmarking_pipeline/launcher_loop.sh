#!/bin/bash


for i in {21..100}; do
python launch_pipeline.py --seed $i --config_template config/template.yaml
done

