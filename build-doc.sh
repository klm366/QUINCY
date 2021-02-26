#!/bin/bash

# Build the documentation using FORD

FORD_MD='project.md'

no_search=""
quiet=""
if [[ $# > 0 ]]; then
  case $1 in
    -f)
      no_search="--no-search"
      ;;
    -q)
      quiet="-q"
      ;;
  esac
fi

output_dir=$(grep output_dir project.md | cut -d' ' -f2)

if hash ford 2>/dev/null; then
  echo "Building documentation ..."
  mkdir pp_temp$$
  python3 scripts/dsl4jsb/dsl4jsb.py -n -d src -t pp_temp$$ -k
  ford $no_search $quiet -d pp_temp$$ ${FORD_MD}
  find $output_dir -type f -name '*.html' -exec sed -i '/github.com\/cmacmackin/s/cmacmackin/Fortran-FOSS-Programmers/' '{}' \;
  rm -rf pp_temp$$
else
  echo "FORD not found! Cannot build documentation"
fi

