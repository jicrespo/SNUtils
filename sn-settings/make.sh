#!/bin/sh
echo -e "Compiling suggest_supernova_settings.cc\n"
g++ suggest_supernova_settings.cc -Wall -o suggest_supernova_settings.exe `root-config --cflags  --glibs`
