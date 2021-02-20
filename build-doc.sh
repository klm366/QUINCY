#!/bin/bash

# Build the documentation using FORD

FORD_MD='jsbach-ford.md'

if hash ford 2>/dev/null; then
  echo "Building documentation ..."
  mkdir pp_temp$$
  python scripts/dsl4jsb/dsl4jsb.py -d src -t pp_temp$$ -k
  ford -d pp_temp$$ ${FORD_MD}
  rm -rf pp_temp$$
else
  echo "FORD not found! Cannot build documentation"
fi

