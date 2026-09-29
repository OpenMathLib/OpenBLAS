#!/usr/bin/env bash
# Deep WASM numerical suite: scalar-oracle CBLAS checks under Node.
# Builds TARGET=WASM128_GENERIC twice (IEEE and relaxed SIMD) and runs the suite.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DIR="$ROOT/test/wasm"
OUT="${OUT:-$DIR/out}"
JOBS="${JOBS:-20}"
COMMON_OPT="${COMMON_OPT:--O2}"

activate_emscripten() {
  if command -v emcc >/dev/null 2>&1; then
    return 0
  fi
  local prefix="${OPENBLAS_EM_PREFIX:-}"
  # Optional in-tree symlink/dir to an emscripten-forge prefix (not committed).
  if [[ -z "$prefix" && -f "$ROOT/.em-prefix/bin/activate_emscripten.sh" ]]; then
    prefix="$ROOT/.em-prefix"
  fi
  if [[ -z "$prefix" || ! -f "$prefix/bin/activate_emscripten.sh" ]]; then
    echo "emcc not found. Set OPENBLAS_EM_PREFIX to an emscripten-forge env," >&2
    echo "or put activate_emscripten.sh under \$ROOT/.em-prefix/bin/." >&2
    exit 1
  fi
  export CONDA_PREFIX="$prefix"
  export PREFIX="${PREFIX:-$prefix}"
  export PATH="$prefix/bin:$prefix/opt/emsdk/upstream/emscripten:$prefix/opt/emsdk/upstream/bin:$PATH"
  export LDFLAGS="${LDFLAGS:-}"
  export CFLAGS="${CFLAGS:-}"
  # shellcheck disable=SC1091
  source "$prefix/bin/activate_emscripten.sh"
}

activate_emscripten
mkdir -p "$OUT"
cd "$ROOT"

SRCS=(
  "$DIR/main.c"
  "$DIR/ref_l1.c"
  "$DIR/ref_l2.c"
  "$DIR/ref_l3.c"
  "$DIR/check_l1.c"
  "$DIR/check_l2.c"
  "$DIR/check_l3.c"
)

LINKFLAGS=(
  -O2
  -msimd128
  -sALLOW_MEMORY_GROWTH=1
  -sSTACK_SIZE=8MB
  -sINITIAL_MEMORY=512MB
  -sMAXIMUM_MEMORY=4GB
  -sEXIT_RUNTIME=1
)

build_and_run() {
  local relaxed="$1"
  local tag
  if [[ "$relaxed" == "1" ]]; then
    tag="relaxed"
  else
    tag="ieee"
  fi

  echo "==> make clean ($tag)"
  # Host leftovers may make `make clean` fail under emcc (-march=native); still wipe objects.
  emmake make clean COLORCODE=0 >/dev/null 2>&1 || true
  find "$ROOT" -name '*.o' -not -path '*/.git/*' -not -path '*/test/wasm/out/*' -delete 2>/dev/null || true
  rm -f "$ROOT"/libopenblas_wasm128*.a "$ROOT"/libopenblas.a

  echo "==> OpenBLAS libs TARGET=WASM128_GENERIC WASM_RELAXED_SIMD=$relaxed"
  emmake make libs \
    -j"$JOBS" \
    HOSTCC=gcc \
    CC=emcc \
    AR=emar \
    RANLIB=emranlib \
    TARGET=WASM128_GENERIC \
    USE_THREAD=0 \
    WASM_RELAXED_SIMD="$relaxed" \
    NOFORTRAN=1 \
    NO_LAPACK=1 \
    NO_LAPACKE=1 \
    COMMON_OPT="$COMMON_OPT" \
    COLORCODE=0

  local lib
  lib=$(ls -1 "$ROOT"/libopenblas_wasm128*.a | head -n1)
  if [[ -z "$lib" ]]; then
    echo "libopenblas_wasm128*.a not found" >&2
    exit 1
  fi

  local defs=()
  local lf=("${LINKFLAGS[@]}")
  if [[ "$relaxed" == "1" ]]; then
    defs+=(-DTEST_WASM_RELAXED)
    lf+=(-mrelaxed-simd)
  fi

  local js="$OUT/numerical_${tag}.js"
  echo "==> linking $js"
  emcc "${lf[@]}" "${defs[@]}" -I"$ROOT" -I"$DIR" -o "$js" "${SRCS[@]}" "$lib"

  echo "==> node $js"
  node "$js"
}

echo "test/wasm run.sh JOBS=$JOBS"
# If WASM_RELAXED_SIMD is set to 0 or 1, run only that mode (used by CI matrix).
# Otherwise run both IEEE and relaxed builds.
case "${WASM_RELAXED_SIMD-}" in
  0)
    build_and_run 0
    echo "test/wasm: IEEE suite passed (WASM_RELAXED_SIMD=0)"
    ;;
  1)
    build_and_run 1
    echo "test/wasm: relaxed suite passed (WASM_RELAXED_SIMD=1)"
    ;;
  "")
    build_and_run 0
    build_and_run 1
    echo "test/wasm: IEEE + relaxed suites passed"
    ;;
  *)
    echo "WASM_RELAXED_SIMD must be unset, 0, or 1 (got: ${WASM_RELAXED_SIMD})" >&2
    exit 1
    ;;
esac
