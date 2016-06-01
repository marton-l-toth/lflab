#!/bin/bash
htx=$1
ds=""; ln=0
for a in $(head -1 $htx); do ds="$ds, \"$a\""; done
echo "const char * hlpn_dir[] = {\"\"$ds, 0};"
echo "const char * hlpn_box[] = {"
grep -n '^{{{' $htx | sed 's/\([0-9]*\):{{{\(.\)\([^}]*\)}}}/"\2_\1\3",/'
echo "0};"
