#!/bin/bash

mkdir -p $HOME/bin

curl -sSL https://github.com/rogeryoungh/hugo-patch/releases/download/latest/hugo-extended-linux-amd64.tar.gz | tar -xz --directory=$HOME/bin

export PATH="$HOME/bin:$PATH"

cd docs && hugo
