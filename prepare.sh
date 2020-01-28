#!/bin/bash

git clone https://github.com/quarkslab/NFLlib.git external/NFLlib
cd external/NFLlib || return
mkdir build
cd build || return
cmake .. -DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=ON
# cmake .. -DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=OFF
make -j
sudo make install
cd ../..

git clone https://github.com/fionser/YELL.git external/YELL
cd external/YELL || return
cmake . -DCMAKE_BUILD_TYPE=Release
make -j
cd ../..
