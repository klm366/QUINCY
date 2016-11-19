#!/bin/ksh 

set -eu

path=$1
file=$2

pftp << EOF
cd ${path}
get ${file} 
EOF
