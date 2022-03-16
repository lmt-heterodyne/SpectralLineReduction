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
lmtinfo.py $obsnum | grep '# ifproc' | cut -c3- > $tmp.rc
source $tmp.rc
rm $tmp.rc

if [ -e "$ifproc" ]; then
    echo "$ifproc"
    ncdump -l 256 $ifproc | grep '^ Header' | column -s= -t 
else
    echo $ifproc does not exist
fi
