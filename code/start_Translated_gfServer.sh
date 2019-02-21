#!/usr/bin/env bash

spc=$1
gfServ=$(which gfServer)

if [ -z "$gfServ" ]
then
      wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/blat/gfServer -P "$CONDA_PREFIX""/bin/" && \
      chmod +x "$CONDA_PREFIX""/bin/gfServer"
fi

if [ -z "$BLATDB" ]
then
    BLATDB="$PWD/data/2bit"
fi

if [ ! -d "$BLATDB" ]
then
    (>&2 echo "No BLATDB variable was found in environment, and no '$BLATDB'  folder was identified. Please set the directory for all 2bit files and try again!" ); exit 1
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
    (>&2 echo "Starting gfServer at localhost:$tPort for file $BLATDB/$twoBitFile" )
    ($gfServ start localhost "$tPort" -canStop -trans -mask "$BLATDB"/"$twoBitFile" &> /dev/null) &
    (>&2 echo "Waiting 5 minutes for gfServer to activate..." )
    sleep 5m
    secs=$((SECONDS+3600))
    while (( SECONDS < secs )); do
        (>&2 echo "Checking to see if gfServer is live" )
        $($gfServ status localhost $tPort &> /dev/null) && break || \
            ( (>&2 echo "Waiting 5 minutes for gfServer to activate..." ) && sleep 5m)
    done
fi
$($gfServ status localhost $tPort &> /dev/null) || \
        { >&2 echo "[T: 60m] Something's wrong with the gfServer!"; exit 1; }

touch "flags/translated-$spc-$tPort"; exit 0