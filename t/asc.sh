#!/bin/bash
if [[ -z "$1" ]]; then exec xterm -font -misc-fixed-medium-r-normal--18-120-100-100-c-90-iso8859-1 -geometry 23x20 -e "$0 x"; fi
while true; do
	for i in $(seq 32 50); do for j in $(seq $i 19 126); do
		printf " %0$(expr 2 + '(' $j '>' 88 ')')d" $j; echo -en "\e[1;32m\x$(printf %02x $j)\e[0m"
	done; echo; done
	read qw; clear
	for i in $(seq 32 50); do for j in $(seq $i 19 126); do
		printf " % $(expr 2 + '(' $j '>' 88 ')')x" $j; echo -en "\e[1;32m\x$(printf %02x $j)\e[0m"
	done; echo; done
	read qw; clear
done

