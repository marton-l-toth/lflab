#!/bin/bash
. $(dirname $(readlink -f $0))/tdir.sh
for d in $(ls "$LF_TMPROOT" | grep 'lf\.[0-9]*$'); do
	dir="$LF_TMPROOT/$d"
	echo "cleaning up $dir"
	fuser -k "$dir/killer-file"
	rm -rf "$dir"
done
