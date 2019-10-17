#!/bin/sh
set -eu

BASE_DIR=$(realpath "$(dirname "$0")/..")
BUILD_DIR=${BUILD_DIR:-"${BASE_DIR}/../build-gatbl"}
BUILD_DIR=$(realpath "$BUILD_DIR")

COMPILERS=${COMPILERS:-"guic gcc clang"}
CXX_STANDARDS=${CXX_STANDARDS:-"11 14 17"}
CMAKE_GENERATOR=$(command -v ninja > /dev/null && echo "Ninja" || echo "Unix Makefiles")
CFLAGS=${CFLAGS:-"-O2 -g"}
CXXFLAGS=${CXXFLAGS:-$CFLAGS}

mkdir -p "$BUILD_DIR"


for comp in $COMPILERS; do
    CC=$(command -v "$comp") || continue
    CXX="${CC%cc}++"

    for standard in $CXX_STANDARDS; do
        CURRENT_BUILD_DIR="$comp-$standard"
        echo "Testing for ${CURRENT_BUILD_DIR}..."
        cd "$BUILD_DIR"
        mkdir -p "$CURRENT_BUILD_DIR"
        cd "$CURRENT_BUILD_DIR"

        cmake "$BASE_DIR" -G "$CMAKE_GENERATOR" \
            -DCMAKE_C_FLAGS="$CFLAGS" -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
            -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_C_COMPILER="$CC" \
            -DCMAKE_CXX_STANDARD="$standard" || continue
        cmake --build . -j --target tests || continue
        cmake --build . -j --target test
        echo
    done
done

