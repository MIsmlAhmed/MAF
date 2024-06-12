#!/bin/bash
# load the needed modules
module restore scimods
# load python environments (that contains all needed packages)
source $HOME/virtual-envs/scienv/bin/activate
#run the model agnostic framework
./model-agnostic.sh model-agnostic.json