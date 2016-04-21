#!/bin/bash
function he() {
	ans="";
	while [[ "$ans" != "y" ]]; do echo "$1 (type \"y\")"; read ans; done
}

function oops() {
	echo "FATAL: $1"; exit 1;
}

function run1() {
	"$TOOL" >>$LF_TMPROOT/lf.*/A <$TCDIR/$1
}

TDIR="$(dirname $(readlink -f "$0"))"
TCDIR="$TDIR/tc"
WDIR="$HOME/.lflab"
TOOL="$HOME/lf.cmd.playback"
cc -g "-o$TOOL" "$TDIR/cmdplay.c" || oops "cmdplay compile failed"
if [[ "$1" == "-k" ]]; then shift; else rm -r "$WDIR.old"; mv "$WDIR" "$WDIR.old"; fi
if [[ -z "$1" ]]; then LAUNCH="/usr/bin/lflab"; else LAUNCH="$1"; fi
LF_TMPROOT="$(readlink -f /run/shm)"
[[ -d "$LF_TMPROOT" ]] || LF_TMPROOT="$(readlink -f /tmp)"
echo "tdir=\"$TDIR\" launch=\"$LAUNCH\""
$LAUNCH &
sleep 1; he "configure audio"
run1 t_0
he "write wav(2,2) and flac(all)"
ls -ltr "$HOME" | tail -3
echo '^V57$F_01c1' >> $LF_TMPROOT/lf.*/A
he "change aodir"
echo '^V57$F>c1' >> $LF_TMPROOT/lf.*/A
sleep 0.5; AODIR="$(grep LF_WAV_DIR $HOME/.lflab/lf.ini | cut '-d"' -f2)"
[[ -z "$AODIR" ]] && oops "LF_WAV_DIR undefined/empty, did you set it?"
run1 t_rec
sleep .5
echo "found AODIR: \"$AODIR\""; ls -l "$AODIR"
he "check wav/flac files"
run1 t_ro
he "check results"
run1 t_err
he "check results"
echo '^V7$.mc1' >> $LF_TMPROOT/lf.*/A
sleep .3; echo '^m:exit(autosv)' >> $LF_TMPROOT/lf.*/A
he "check autosaves & log"
$LAUNCH &
sleep .5; run1 t_start
he "recover autosave"
run1 t_trkdep0
he "test trk editor (drag while playing)"
run1 t_common
$TOOL -a-1$TCDIR >> $LF_TMPROOT/lf.*/A
echo '^V7$.sc1' >> $LF_TMPROOT/lf.*/A
he "save"
echo '^V7$.mc1' >> $LF_TMPROOT/lf.*/A
sleep .3; echo '^m:exit(autosv)' >> $LF_TMPROOT/lf.*/A
sleep .5; $LAUNCH &
sleep 1; run1 t_start
he "load saved file"
echo '^V7$.2c9' >> $LF_TMPROOT/lf.*/A
$TOOL -a-a$TCDIR >> $LF_TMPROOT/lf.*/A
he "check (no err, t[A-Z] removed) -- done!"
