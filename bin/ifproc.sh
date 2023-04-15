#! /usr/bin/env bash
#

obsnum=$1

if [ -z "$obsnum" ]; then
    echo usage: $0 obsnum
    echo Reports the ifproc Header variables
    exit 0
fi

tmp=/tmp/tmp$$.rc

ifproc="obsnum $obsnum"
lmtinfo.py $obsnum > $tmp.rc
source $tmp.rc
rm $tmp.rc

if [ ! -z "$rawnc" ]; then
    echo "$rawnc"
    ncdump -l 256 $rawnc | grep '^ Header' | column -s= -t 
else
    echo $ifproc does not exist
fi
