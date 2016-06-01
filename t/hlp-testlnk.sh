#!/bin/bash
LF_TMPROOT="$(readlink -f /run/shm)"
[[ -d "$LF_TMPROOT" ]] || LF_TMPROOT="$(readlink -f /tmp)"
qw=$1
[[ -z "$qw" ]] && qw="/usr/share/lflab"
nm="$qw/help.txt"
echo $nm
sed -n 's/^==> //p' $nm | sed 's/^\([^.]\)/.!b.?.\1/;s/^/^^I/;s/ -- .*$//;s/$/\$0/' > $LF_TMPROOT/lf.*/A
