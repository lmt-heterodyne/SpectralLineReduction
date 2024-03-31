#! /usr/bin/env bash
#

obsnum=$1

if [ -z "$obsnum" ]; then
    echo usage: $0 obsnum
    echo "Reports the ifproc Header variables for SEQ."
    echo "For RSR it actually also reports the netCDF variables of the RedshiftChassis0 file"
    exit 0
fi

tmp=/tmp/tmp$$.rc

ifproc="obsnum $obsnum"
lmtinfo.py $obsnum > $tmp.rc
if [ $? -ne 0 ]; then
    echo "ifproc.sh:   $obsnum not a valid obsnum"
    exit 0
fi
source $tmp.rc
rm $tmp.rc

if [ ! -z "$rawnc" ]; then
    echo "$rawnc"
    ncdump -l 256 $rawnc | grep '^ Header' | column -s= -t 
else
    echo $ifproc does not exist
fi
