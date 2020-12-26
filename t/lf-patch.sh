#!/bin/bash

GIT_CMD=git
if [[ -n "$1" ]]; then GIT_CMD="$1"; fi
BIN_DIR="$HOME/bin"
LFDIR="$HOME/g/lflab"
CURR_DATE="$(date +%s)"
PROD_BASE="$HOME/lf-patch"

PROD_DIR="$PROD_BASE/$CURR_DATE"
BIN_BASE="$BIN_DIR/lflab"

if [[ -e "$PROD_DIR" ]]; then
	echo "$PROD_DIR already exists -- you are either very quick, or you found a bug";
	echo "nothing changed"; exit 1
fi
mkdir -p "$PROD_DIR"

cd "$LFDIR/s"
if $GIT_CMD pull ; then
	echo "Download OK"
else
	echo "Download failed, nothing changed"; exit 1
fi

PATCH_LS=$(echo $(git log -n5 | grep '[a-zA-Z0-9]' | egrep -v '^(Author:|commit [0-9a-f]|Date:)' | sed 's/^ */|/;s/ *$//'))
PATCH_STR="patch: $(date -d @$CURR_DATE +%F,%T)$PATCH_LS"

if ../t/prod.sh -j$(nproc) -pstr "$PATCH_STR" -pd "$PROD_DIR"; then
	echo "compile OK"
else
	echo "compile failed, nothing changed"; exit 1
fi
	
if [[ -x "$BIN_BASE-5" ]]; then
	lnk=$(readlink -f $BIN_BASE-5)
	if echo $lnk | grep "$PROD_BASE/[0-9]*/lf.bin" | grep -qv $CURR_DATE; then
		olddir=$(dirname $lnk)
		echo "removing old version '$olddir'"
		rm -rf $olddir
	else
		echo "$BIN_BASE-5 points to '$lnk', not removing"
	fi
	rm "$BIN_BASE-5"
fi

OLD_LS="";
for j in $(seq 5 -1 1); do
	k=$(expr $j - 1)
	if [[ -x "$BIN_BASE-$k" ]]; then mv "$BIN_BASE-$k" "$BIN_BASE-$j"; OLD_LS="lflab-$j $OLD_LS"; fi
done
if ln -sf "$PROD_DIR/lf.bin" $BIN_BASE-0 && ln -sf "$PROD_DIR/lf.bin" $BIN_BASE; then
	echo "patch done (lflab), old versions: $OLD_LS"
else
	echo "hmmm... some strange error happened; good luck anyway"
fi
