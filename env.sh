#!/bin/bash

module load gcc
module load cmake
path=$(which gcc)
export CC=$path
module load R
module load git
module load valgrind
ssh -T git@github.com