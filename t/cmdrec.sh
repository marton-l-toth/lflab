#!/bin/bash
if [[ -z "$1" ]]; then logf="$HOME/.lflab/lf.log"; else logf="$1"; fi
grep '\. cmd_rec:' "$logf" | sed 's/^.*\. cmd_rec://'
