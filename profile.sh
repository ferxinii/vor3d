#!/bin/bash
CPUPROFILE=$1.prof DYLD_INSERT_LIBRARIES=/opt/homebrew/Cellar/gperftools/2.16/lib/libprofiler.dylib ./$1
pprof --pdf $1 $1.prof > $1.pdf
echo "Profiling results: $1.pdf"
