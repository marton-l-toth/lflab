#!/bin/bash
grep 'save clk:' | cut -d: -f2 | awk 'BEGIN{a=b=c=d=0}END{print a " " b " " c " " d " " a+b+c+d}{a+=$1;b+=$2;c+=$3;d+=$4}'
