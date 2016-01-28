#!/bin/bash
SDIR="$HOME/g/lflab/s"
PDIR="$HOME/lf-prod/s"
mkdir -p "$PDIR/c"
OKCMD=""; FAILCMD="";
if [[ "$1" == "-okcmd"   ]]; then OKCMD="$2";   shift 2; fi
if [[ "$1" == "-failcmd" ]]; then FAILCMD="$2"; shift 2; fi
cd "$SDIR"
for a in * c/*; do
	b="$PDIR/$a"
	if [[ -d "$a" ]]; then continue; fi
	if [[ ! -f "$b" ]] || [[ "$a" -nt "$b" ]]; then cp $a $b; fi
done
cd "$PDIR"
CP=$(echo /run/shm/lf.*/A)
if [[ -z "$OKCMD$FAILCMD" ]]; then exec make $*; fi
if [[ ! -p "$CP" ]]; then echo "no pipe found ('$CP')"; exec make $*; fi
if make $*; then
	[[ -n "$OKCMD"   ]] && (echo "$OKCMD"   > "$CP")
else
	[[ -n "$FAILCMD" ]] && (echo "$FAILCMD" > "$CP")
fi
