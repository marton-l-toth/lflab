#!/bin/bash
if [[ -z "$1" ]]; then echo "usage: $0 targetdir"; exit 1; fi
mkdir -p $1
SDIR=$HOME/lf-prod/s
cp $SDIR/*rpm $SDIR/*deb $1/
tar -czvf $1/debugbin.tgz $SDIR/lf.bin $SDIR/lf.gui $SDIR/lf.bb $SDIR/lf.pump
