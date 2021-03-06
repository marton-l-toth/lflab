#!/bin/bash
if [[ ! -f "c/upd_contrib.sh" ]]; then 
	echo "$0: please run from the source directory"; exit 1; fi
cd c
if [[ -n "$1" ]]; then
	for cc in [a-zA-Z]*.cc; do
		out="$out $2${cc:0:(-3)}.$1"
	done; echo $out; exit 0;
fi
stamp="generated by $0 at $(date)"
head1="// $stamp\n#include \"contrib.h\""
head2="contrib_ent contrib_list[] = {"
stl=""
dcl=""
for cc in [a-zA-Z]*.cc; do
	bn=${cc:0:(-3)}
	fn="ci_${bn}"
	dcl="$dcl,$fn"
	stl="$stl {\"$bn\", &$fn}, "
done
echo -e "$head1\nc_inifun ${dcl:1};\n$head2\n$stl {0,0}};" > ../contrib.cc
