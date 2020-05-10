#!/usr/bin/env bash

spc=$1
outfile=$2
twoBitDir=$3
gfServ=$(which gfServer)

if [ -z "$gfServ" ]
then
      wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/blat/gfServer -P "$CONDA_PREFIX""/bin/" && \
      chmod +x "$CONDA_PREFIX""/bin/gfServer"
fi

port_table="./data/portTable.csv"

if [ ! -f "$port_table" ]
then
    (>&2 echo "No port table file found in $port_table!" ); exit 1
fi

gfinfo=$(cat $port_table | \
    sed "s/\"//g; /,,/d; /----/d; /^#.*$/d; s/^[0-9]* *$//g; /^$/d" | \
    grep "$spc")
tPort=$(echo $gfinfo | cut -d"," -f 1)
twoBitFile="$(echo $gfinfo | cut -d"," -f 3)".2bit

if $($gfServ status localhost $tPort &> /dev/null); then
    (>&2 echo "gfServer for $gfinfo already running" )
else
    (>&2 echo "Starting gfServer at localhost:$tPort for file $twoBitDir/$twoBitFile" )
    ($gfServ start localhost "$tPort" -canStop -trans -mask "$twoBitDir"/"$twoBitFile" &> /dev/null) &
    (>&2 echo "Waiting 1 minute for gfServer to activate..." )
    sleep 1m
    secs=$((SECONDS+3600))
    while (( SECONDS < secs )); do
        (>&2 echo "Checking to see if gfServer is live" )
        $($gfServ status localhost $tPort &> /dev/null) && break || \
            ( (>&2 echo "Waiting 1 minute for gfServer to activate..." ) && sleep 1m)
    done
fi
$($gfServ status localhost $tPort &> /dev/null) && \
	echo "gfServer is live at localhost $tPort! " || \
        { >&2 echo "[T: 60m] Something's wrong with the gfServer!"; exit 1; }

touch "$outfile" && exit 0
