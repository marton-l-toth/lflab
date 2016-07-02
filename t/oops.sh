#!/bin/bash
git stash -a; ../t/prod.sh -rsm && echo "--- prod. OK" ; git stash pop
