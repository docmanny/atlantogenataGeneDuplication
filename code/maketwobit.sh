#!/bin/bash

makeTwoBitIfNotExist(){
	if [[ ! -e $2 ]]; then
	    faToTwoBit $1 $2
	fi
}

export -f makeTwoBitIfNotExist

ls /usr/db/genomes/*.fa | parallel makeTwoBitIfNotExist {} /usr/db/BLAT/{/.}.2bit
