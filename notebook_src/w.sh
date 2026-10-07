#!/bin/bash
# usage: w.sh target.wl [args]
B="$(cd "$(dirname "$0")" && pwd)"
"/c/Program Files/Wolfram Research/Wolfram/15.0/wolfram.exe" -noinit -noprompt -script "$(cygpath -w "$B/run.wl")" "$@" 2>&1
