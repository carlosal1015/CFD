#!/usr/bin/env bash

DIR=$(pwd)/build
if [ -d "$DIR" ]; then
  printf '%s\n' "Removing Lock ($DIR)"
  rm -rf "$DIR"
fi

cmake -S . -B build
cmake --build build

pushd build/
./FTCS
./BTCS
popd

./readu.py
