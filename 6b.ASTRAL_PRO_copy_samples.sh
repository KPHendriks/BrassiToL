#!/bin/bash

sample=$1

#Copy sample specific mapped HybPiper results.
cp mapped_Nikolov2019/${sample}_mapped_Nikolov2019.tar.gz /tmp/paralog_retrieval/

#Move to the directory.
cd /tmp/paralog_retrieval/

#Unpack
tar -xf ${sample}_mapped_Nikolov2019.tar.gz
