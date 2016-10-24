#!/bin/sh
export LF_ZZZARGS="$*"
set | grep "LF_"
exec /bin/bash
