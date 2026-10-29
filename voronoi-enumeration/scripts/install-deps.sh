#!/usr/bin/env bash
# Build pinned dependencies locally; no sudo or system installation.
set -euo pipefail
project_dir=$(cd "$(dirname "$0")/.." && pwd)
deps_dir="$project_dir/.deps"
if [ "$(uname -s)" = Darwin ]; then
    task_sdk=$(xcrun --show-sdk-path)
    if [ -d "$task_sdk/usr/include/c++/v1" ]; then
        export CXXFLAGS="${CXXFLAGS:-} -isystem $task_sdk/usr/include/c++/v1"
    fi
fi
mkdir -p "$deps_dir/downloads" "$deps_dir/src" "$deps_dir/build" "$deps_dir/prefix"
fetch() {
    local url="$1" archive="$2"
    if [ ! -f "$deps_dir/downloads/$archive" ]; then
        curl --fail --location --http1.1 --retry 5 --retry-all-errors --continue-at - \
            --max-time 600 "$url" -o "$deps_dir/downloads/$archive.part"
        mv "$deps_dir/downloads/$archive.part" "$deps_dir/downloads/$archive"
    fi
}
fetch https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz gmp-6.3.0.tar.xz
fetch https://archives.boost.io/release/1.89.0/source/boost_1_89_0.tar.bz2 boost_1_89_0.tar.bz2
fetch https://github.com/scipopt/soplex/archive/refs/tags/v8.0.3.tar.gz soplex-8.0.3.tar.gz
(cd "$deps_dir/downloads" && shasum -a 256 -c "$project_dir/dependencies.sha256")
if [ ! -d "$deps_dir/src/gmp-6.3.0" ]; then
    tar -xf "$deps_dir/downloads/gmp-6.3.0.tar.xz" -C "$deps_dir/src"
fi
if [ ! -d "$deps_dir/src/boost_1_89_0" ]; then
    tar -xf "$deps_dir/downloads/boost_1_89_0.tar.bz2" -C "$deps_dir/src" boost_1_89_0/boost boost_1_89_0/LICENSE_1_0.txt
fi
if [ ! -d "$deps_dir/src/soplex-8.0.3" ]; then
    tar -xf "$deps_dir/downloads/soplex-8.0.3.tar.gz" -C "$deps_dir/src"
fi
if [ ! -f "$deps_dir/prefix/lib/libgmpxx.a" ]; then
    mkdir -p "$deps_dir/build/gmp"
    cd "$deps_dir/build/gmp"
    "$deps_dir/src/gmp-6.3.0/configure" --prefix="$deps_dir/prefix" --enable-cxx --disable-shared --with-pic
    make -j4
    make check -j4
    make install
fi
cmake -S "$deps_dir/src/soplex-8.0.3" -B "$deps_dir/build/soplex" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$deps_dir/prefix" \
    -DCMAKE_MODULE_PATH="$project_dir/cmake" \
    -DCMAKE_PREFIX_PATH="$deps_dir/prefix" -DGMP=on -DMPFR=off -DBOOST=on \
    -DBoost_INCLUDE_DIR="$deps_dir/src/boost_1_89_0" -DBUILD_SHARED_LIBS=off
cmake --build "$deps_dir/build/soplex" -j4
cmake --install "$deps_dir/build/soplex"
if ! grep -q '^#define SOPLEX_WITH_BOOST' "$deps_dir/prefix/include/soplex/config.h"; then
    echo "ERROR: SoPlex was built without the required exact rational support" >&2
    exit 1
fi
cd "$deps_dir/downloads"
shasum -a 256 gmp-6.3.0.tar.xz boost_1_89_0.tar.bz2 soplex-8.0.3.tar.gz > "$deps_dir/SHA256SUMS"
echo "Dependencies installed in $deps_dir/prefix"
