#!/bin/bash
if [[ -z "$1" ]] || [[ -z "$2" ]]; then echo "usage: $0 <cfgtab> <sh-src>"; exit 1; fi
gl=$(grep -n '#=#=#=#' "$2" | cut -d: -f1)
head -$(expr $gl - 1) "$2"
echo "#####begin generated part ($0 at $(date))"
awk '-F|' '{gsub(/ *$/,"",$3); gsub(/ *$/,"",$1); print "                     echo \"    -" $6 ($2=="s"?"s":"x") "   " $7 " [" $3 "/\x24LF_" $1 "]"  "\""}' < "$1"
echo "                    exit 0;;"
sed 's/ *|/|/g' "$1" | awk '-F|' '{print "                 \"-" $6 "\") LF_" $1 "=\"$2\"; shift 2;; "}'
echo '                 -*) echo "unknown option $1"; shift;;'
echo '                 *) cond="";; '
echo -e "        esac\ndone\n"
vls=$(cut '-d|' -f1 $1 | sed 's/^/LF_/g')
echo "export" $vls
echo "#####end generated part ($0 at $(date))"
tail -$(expr $(wc -l < "$2") - $gl) "$2"
