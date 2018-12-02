#!/bin/bash
SDIR="$PWD"
PDIR="$HOME/tmp/lf-prod/s"
OCMD=""; FCMD=""; JARG=""; cond="y"; SSE="i"; UCFLG=""; DEFS=""; STFLG="";
CGARG="-fno-exceptions -fno-rtti"
WNARG="-Wall -Wno-unused-function -Wno-misleading-indentation"
OPTARG="-O2 -fpredictive-commoning -fgcse-after-reload"
while [[ -n "$cond" ]]; do
        case "$1" in
                "")   cond="" ;;
                "--") cond=""; shift ;;
                -j*)  JARG=$1; shift ;;
		"-BB") rm -r "$PDIR"; shift;;
                "-S") SSE="-mno-sse2";  shift ;;
                "-s") SSE="-msse2";  shift ;;
                "-st") STFLG="y"; shift ;;
		"-md") DEFS="$DEFS -DLF_C_MEMDEBUG"; shift ;;
		"-0") OPTARG="";    shift ;;
		"-0g") OPTARG="-Og";    shift ;;
		"-cf")  UCFLG="$2"; shift 2 ;;
		"-ocmd") OCMD="$2"; shift 2 ;;
		"-fcmd") FCMD="$2"; shift 2 ;;
		"-rsm")  OCMD="^q2"; FCMD="^C47\$E>00ffffec"; shift ;;
		"-rsg")  OCMD="^q";  FCMD="^C47\$E>00ffffeb"; shift ;;
		"-sd") SDIR="$2"; shift 2 ;;
		"-pd") PDIR="$2"; shift 2 ;;
		*) cond="" ;;
	esac
done
mkdir -p "$PDIR/c" "$PDIR/co" || exit 1
[[ -z "$JARG" ]] && JARG=-j$(expr $(nproc) + 1)
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
if [[ -n "$STFLG" ]]; then
	LF_TMPROOT="$(readlink -f /run/shm)"
	[[ -d "$LF_TMPROOT" ]] || LF_TMPROOT="$(readlink -f /dev/shm)"
	[[ -d "$LF_TMPROOT" ]] || LF_TMPROOT="$(readlink -f /tmp)"
	stfil=$LF_TMPROOT/mkstat-$$
	make stat | awk  '/^ *[0-9]/{print $3 " print \"" $0 "\\t(\"; print eXpR1 " $3 " eXpR2; print \"%)\\n\"" }' | sort -n > $stfil
	total=$(tail -1 $stfil | grep -o '^[0-9]*')
	(echo scale=3; sed "s/^[0-9]* //;s/eXpR1/(100.0*/;s/eXpR2/+$total\/2000)\/$total/" $stfil) | bc
	rm $stfil; exit 0
fi
CP=$(echo /run/shm/lf.*/A)
CFLG="CFLAGS=$SSE -pipe $CGARG $WNARG $OPTARG $UCFLG$DEFS"
echo make $JARG "\"$CFLG\"" $*
if [[ -z "$OCMD$FCMD" ]]; then exec make $JARG "$CFLG" $*; fi
if [[ ! -p "$CP" ]]; then echo "no pipe found ('$CP')"; exec make $JARG "$CFLG" $*; fi
if make $JARG "$CFLG" $*; then
	[[ -n "$OCMD" ]] && (echo "$OCMD" > "$CP")
else
	[[ -n "$FCMD" ]] && (echo "$FCMD" > "$CP")
	exit 1
fi
