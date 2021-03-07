#!/bin/bash
set -e

toolchain_path="C:/msys64/mingw64/bin/"

if [[ ! -d "build" ]]; then
    mkdir build
fi

cpp_files=()
cpp_includes=(
    "-Ideps/matplotlib-cpp"
    "-IC:\msys64\mingw64\include\python3.8"
)
library_paths=(
    "-LC:\msys64\mingw64\lib"
    "-LC:\msys64\mingw64\lib\python3.8"
)

${toolchain_path}clang++ --std=c++2a -Wno-ignored-attributes -fopenmp ${cpp_includes[@]} -g ${cpp_files[@]} example.cpp -lpython3.8 -o build/bin
echo "Compilation done."
echo
echo "Running program..."
./build/bin.exe