#!/bin/bash

c0='\e[0m'
c1='\e[31;1m'
c2='\e[33;1m'
c3='\e[32;1m'
c4='\e[35;1m'
if [[ -z "$1" ]]; then
	a1='========'
	a2='[Mc][0-9()]*:.* [Mc][0-9()]*:'
	a3='M[0-9()]*:'
	a4='c[0-9()]*:'
else
	a1="$1"
	a2="$2"
	a3="$3"
	a4="$4"
fi

while read qw; do
	if   [[ -n "$a1" ]] && [[ "$qw" =~ $a1 ]]; then echo -en "$c1"
	elif [[ -n "$a2" ]] && [[ "$qw" =~ $a2 ]]; then echo -en "$c2"
	elif [[ -n "$a3" ]] && [[ "$qw" =~ $a3 ]]; then echo -en "$c3"
	elif [[ -n "$a4" ]] && [[ "$qw" =~ $a4 ]]; then echo -en "$c4"
	else echo -en "$c0"; fi
	echo $qw
done
echo -en $c0
