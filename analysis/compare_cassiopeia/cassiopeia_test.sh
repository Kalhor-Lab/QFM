#!/bin/bash
ml python/3.6
ml gurobi/9.0.2
source ../../../Cassiopeia/venv/bin/activate
i=$SLURM_ARRAY_TASK_ID
#i=2
python3 -u cassiopeia_testing.py -t $i -hg 50
