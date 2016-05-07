#!/bin/bash
SDIR="$HOME/g/lflab/s"
PDIR="$HOME/lf-prod/s"
mkdir -p "$PDIR/c"
OCMD=""; FCMD=""; JARG=""; cond="y"; SSE="i"; UCFLG="";
OPTARG="-O2 -fpredictive-commoning -fgcse-after-reload"
while [[ -n "$cond" ]]; do
        case "$1" in
                "")   cond="" ;;
                "--") cond=""; shift ;;
                -j*)  JARG=$1; shift ;;
                "-S") SSE="";  shift ;;
                "-s") SSE="-msse2";  shift ;;
		"-0") OPTARG="";    shift ;;
		"-cf")  UCFLG="$2"; shift 2 ;;
		"-ocmd") OCMD="$2"; shift 2 ;;
		"-fcmd") FCMD="$2"; shift 2 ;;
		"-rsm")  OCMD="^q2"; FCMD="^C47\$E>00ffffec"; shift ;;
		"-rsg")  OCMD="^q";  FCMD="^C47\$E>00ffffeb"; shift ;;
		*) cond="" ;;
	esac
done
[[ -z "$JARG" ]] && JARG=-j$(expr $(grep '^processor[[:space:]]*:' /proc/cpuinfo | wc -l) + 1)
if [[ "$SSE" == "i" ]]; then
	SSE=""; grep -q 'flag.*sse2' /proc/cpuinfo && SSE="-msse2"
fi

cd "$SDIR"
for a in * c/*; do
	b="$PDIR/$a"
	if [[ -d "$a" ]]; then continue; fi
	if [[ ! -f "$b" ]] || [[ "$a" -nt "$b" ]]; then cp $a $b; fi
done
cd "$PDIR"
CP=$(echo /run/shm/lf.*/A)
CFLG="CFLAGS=$SSE $OPTARG $UCFLG"
echo make $JARG "\"$CFLG\"" $*
if [[ -z "$OCMD$FCMD" ]]; then exec make $JARG "$CFLG" $*; fi
if [[ ! -p "$CP" ]]; then echo "no pipe found ('$CP')"; exec make $JARG "$CFLG" $*; fi
if make $JARG "$CFLG" $*; then
	[[ -n "$OCMD" ]] && (echo "$OCMD" > "$CP")
else
	[[ -n "$FCMD" ]] && (echo "$FCMD" > "$CP")
fi
