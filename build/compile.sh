#!/usr/bin/bash

clear
mkdir build
cd build
cmake ~/01_util/02_FPCSL/02_projects/1d2psolver/build
cmake --build .
mv 1D2Psolver.x ..
rm ~/01_util/02_FPCSL/02_projects/1d2psolver/src/compile_commands.json
ln -s ~/01_util/02_FPCSL/02_projects/cases/test/build/compile_commands.json ~/01_util/02_FPCSL/02_projects/1d2psolver/src
cd ..
