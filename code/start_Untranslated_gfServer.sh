#!/usr/bin/env bash

spc=$1
outfile=$2
twoBitDir=$3
port_table="data/portTable.csv"
gfServ=$(which gfServer)

if [ -z "$gfServ" ]
then
      wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/blat/gfServer -P "$CONDA_PREFIX""/bin/" && \
      chmod +x "$CONDA_PREFIX""/bin/gfServer"
fi

if [ ! -f "$port_table" ]
then
    (>&2 echo "No port table file found in $PWD/$port_table!" ); exit 1
fi

gfinfo=$(cat $port_table | \
    sed "s/\"//g; /,,/d; /----/d; /^#.*$/d; s/^[0-9]* *$//g; /^$/d" | \
    grep "$spc")
utPort=$(echo $gfinfo | cut -d"," -f 2)
twoBitFile="$(echo $gfinfo | cut -d"," -f 3)".2bit

if $($gfServ status localhost $utPort &> /dev/null); then
    (>&2 echo "gfServer for $gfinfo already running" )
else
    (>&2 echo "Starting gfServer at localhost:$utPort for file $twoBitDir/$twoBitFile" )
    ($gfServ start localhost "$utPort" -canStop -stepSize=5 -repMatch=2253 "$twoBitDir/$twoBitFile" &> /dev/null) &
    (>&2 echo "Waiting 1 minute for gfServer to activate..." )
    sleep 1m
    secs=$((SECONDS+3600))
    while (( SECONDS < secs )); do
        (>&2 echo "Checking to see if gfServer is live" )
        $($gfServ status localhost $utPort &> /dev/null) && break || \
            ( (>&2 echo "Waiting 1 minutes for gfServer to activate..." ) && sleep 1m)
    done
fi
$($gfServ status localhost $utPort &> /dev/null) && \
	echo "gfServer is live at localhost $utPort" || \
        { >&2 echo "[T: 60m] Something's wrong with the gfServer!"; exit 1; }

touch "$outfile" && exit 0
