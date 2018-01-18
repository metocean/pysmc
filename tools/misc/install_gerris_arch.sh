#!/bin/bash

echo "-- Installing libs"
#sudo pacman -S darcs ftgl glib2 glibc m4 ode fftw3 automake libtool
#yaourt hypre

echo "-- Downloading source code"
darcs get http://gerris.dalembert.upmc.fr/darcs/gts-stable
darcs get http://gerris.dalembert.upmc.fr/darcs/gerris-stable
darcs get http://gerris.dalembert.upmc.fr/darcs/gfsview-stable

cd gts-stable
./configure
make -j4
sudo make install

cd ../gerris-stable
./configure "LIBS=-lm"
make -j4
sudo make install

cd ../gfsview-stable
./configure
make -j4
sudo make install

echo "-- Finished"
