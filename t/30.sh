#!/bin/bash
cp $1 $2
for a in $(seq 29); do sed "/^[CENX:]/{s/\!e/\!e$a/g}" $1 >> $2; done
