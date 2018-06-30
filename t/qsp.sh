#!/bin/bash
cmd=$(echo /run/shm/lf.*/A)
if [[ -z "$2" ]]; then sec="1"; else sec="$2"; fi
if [[ -z "$3" ]]; then rpt=$(expr 10 / "$sec"); else rpt="$3"; fi
if [[ -z "$4" ]]; then slp=3; else slp="$4"; fi
if [[ -z "$5" ]]; then pfx="perf"; else pfx="$5"; fi
tot=$(expr 44100 '*' $sec)
for blk in 441 294 210 150 100 70 50 36 25 18; do
	nb=$(expr $tot / $blk)
	echo "X$1\$%$pfx $1 $rpt*$nb*$blk :$rpt $nb $blk" >> $cmd
	sleep $slp
done

