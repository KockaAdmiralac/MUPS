#!/bin/bash

git pull
hash="$(git rev-parse --short HEAD)"
cd "$(dirname "$0")"
cp -r src ../sl190368d
cp run.py ../sl190368d
cp Makefile ../sl190368d
cd ../sl190368d
svn add --force *
svn ci -m "Transferring Git commit $hash."
