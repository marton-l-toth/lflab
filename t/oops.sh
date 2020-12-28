#!/bin/bash
git stash -a; ../t/prod.sh -rsm && echo -e '\e[1;32m--- prod. OK\e[0m'; git stash pop
