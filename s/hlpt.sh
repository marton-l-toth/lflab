#!/bin/bash
txf=$1
lnksed="s/^==> 0./==> /"
lnum=1; ec=0; tc=0; lc=0; cnt=1;
for a in $(head -1 $txf); do
	ds="$ds, \"$a\""
	lnksed="s/^==> $cnt./==> $(echo $a|sed 's/^./&./')./;$lnksed"
	cnt=$(expr $cnt + 1)
done
#echo "lnksed=\"$lnksed\"" >&2
grep -n ^ $1 | sed "s/\([0-9]*\):{{{.*\$/{{{?._\1}}}/;s/^[0-9]*://;$lnksed"
