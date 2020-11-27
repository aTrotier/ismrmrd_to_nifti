#!/bin/bash

if [ $# = 2 ]
then
    docker run --rm -v $1:/data ismrmrd2nifti -i /data/$2 -o /data/out.nii.gz
else
    echo "Usage: ./run.sh /path/to/datadir data.h5"
fi
