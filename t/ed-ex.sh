#!/bin/bash
TDIR="$(dirname $(readlink -f "$0"))"
cp "$TDIR/../s/ex.lf" "$HOME/ex0.lf"
$HOME/lf-prod/s/lf -n "$HOME/ex0.lf" &
