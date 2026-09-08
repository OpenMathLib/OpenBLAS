#!/usr/bin/env bash
# Verify WASM relaxed SIMD is opt-in:
#   * the default compile does not pass -mrelaxed-simd and emits no relaxed_madd
#   * WASM_RELAXED_SIMD=1 passes the flag and kernels use relaxed_madd
#
# Requires a WASM-configured tree (ARCH=wasm in Makefile.conf), emcc, and
# wasm-dis (Binaryen, shipped with Emscripten). Node is used to instantiate
# the default module when available.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
FAILS=0

activate_emscripten() {
  if command -v emcc >/dev/null 2>&1; then
    return 0
  fi
  local prefix="${OPENBLAS_EM_PREFIX:-}"
  if [[ -z "$prefix" && -n "${CONDA_PREFIX:-}" && -f "${CONDA_PREFIX}/bin/activate_emscripten.sh" ]]; then
    prefix="$CONDA_PREFIX"
  fi
  if [[ -z "$prefix" && -n "${EMSDK:-}" && -x "${EMSDK}/upstream/emscripten/emcc" ]]; then
    export PATH="${EMSDK}/upstream/emscripten:${EMSDK}/upstream/bin:${PATH}"
    return 0
  fi
  if [[ -z "$prefix" && -f "$ROOT/.em-prefix/bin/activate_emscripten.sh" ]]; then
    prefix="$ROOT/.em-prefix"
  fi
  if [[ -z "$prefix" || ! -f "$prefix/bin/activate_emscripten.sh" ]]; then
    echo "emcc not found. Set OPENBLAS_EM_PREFIX to an Emscripten prefix." >&2
    exit 1
  fi
  export CONDA_PREFIX="$prefix"
  export PREFIX="${PREFIX:-$prefix}"
  export PATH="$prefix/bin:$prefix/opt/emsdk/upstream/emscripten:$prefix/opt/emsdk/upstream/bin:$PATH"
  # shellcheck disable=SC1091
  source "$prefix/bin/activate_emscripten.sh"
}

count_op() {
  # wasm-dis prints a warning on relocatable objects; keep stdout only.
  local n
  n="$(wasm-dis "$1" 2>/dev/null | grep -c "$2" || true)"
  echo "${n:-0}"
}

expect() {
  local name="$1" got="$2" pred="$3" want="$4"
  if [[ "$pred" == "eq" && "$got" == "$want" ]]; then
    echo "PASS  $name ($got)"
    return 0
  fi
  if [[ "$pred" == "gt" && "$got" -gt "$want" ]]; then
    echo "PASS  $name ($got > $want)"
    return 0
  fi
  echo "FAIL  $name (got $got, expected $pred $want)" >&2
  FAILS=$((FAILS + 1))
}

activate_emscripten
cd "$ROOT"

if ! command -v wasm-dis >/dev/null 2>&1; then
  echo "wasm-dis not found (Binaryen). It is usually on PATH after activating Emscripten." >&2
  exit 1
fi

if [[ ! -f Makefile.conf ]] || ! grep -q '^ARCH=wasm$' Makefile.conf; then
  echo "Makefile.conf is not a WASM config (need ARCH=wasm)." >&2
  echo "Configure first, for example:" >&2
  echo "  make TARGET=WASM128_GENERIC HOSTCC=gcc CC=emcc AR=emar RANLIB=emranlib USE_THREAD=0 NOFORTRAN=1 NO_LAPACK=1" >&2
  exit 1
fi

MAKE_COMMON=(
  HOSTCC=gcc
  CC=emcc
  AR=emar
  RANLIB=emranlib
  TARGET=WASM128_GENERIC
  USE_THREAD=0
  NOFORTRAN=1
  NO_LAPACK=1
  COLORCODE=0
)

# Kernels that emit madd when __wasm_relaxed_simd__ is set.
OBJS=(
  sgemm_kernel.o
  dgemm_kernel.o
  cgemm_kernel_n.o
  zgemm_kernel_n.o
  sdot_k.o
  ddot_k.o
  dsdot_k.o
  srot_k.o
  dtrmm_kernel_LN.o
)

EXPORTS=_sgemm_kernel,_dgemm_kernel,_cgemm_kernel_n,_zgemm_kernel_n,_sdot_k,_ddot_k,_dsdot_k,_srot_k,_dtrmm_kernel_LN

MAKE_BIN="${MAKE:-make}"

sgemm_line() {
  local extra=("$@")
  rm -f kernel/sgemm_kernel.o
  # GNU make -n prints the compile line only when the target is out of date.
  "$MAKE_BIN" -C kernel -n sgemm_kernel.o "${MAKE_COMMON[@]}" "${extra[@]}" 2>/dev/null \
    | grep -E '(^|[[:space:]])emcc[[:space:]]' | tail -1 || true
}

echo "==> compile-line flags"
DEFAULT_LINE="$(sgemm_line)"
RELAXED_LINE="$(sgemm_line WASM_RELAXED_SIMD=1)"
echo "  default: $DEFAULT_LINE"
echo "  relaxed: $RELAXED_LINE"
if [[ -z "$DEFAULT_LINE" || -z "$RELAXED_LINE" ]]; then
  echo "FAIL  could not extract an emcc compile line from make -n" >&2
  exit 1
fi

case "$DEFAULT_LINE" in
  *-msimd128*) expect "default passes -msimd128" 1 eq 1 ;;
  *) expect "default passes -msimd128" 0 eq 1 ;;
esac
case "$DEFAULT_LINE" in
  *-mrelaxed-simd*) expect "default omits -mrelaxed-simd" 1 eq 0 ;;
  *) expect "default omits -mrelaxed-simd" 0 eq 0 ;;
esac
case "$RELAXED_LINE" in
  *-mrelaxed-simd*) expect "WASM_RELAXED_SIMD=1 passes -mrelaxed-simd" 1 eq 1 ;;
  *) expect "WASM_RELAXED_SIMD=1 passes -mrelaxed-simd" 0 eq 1 ;;
esac

OUT="$(mktemp -d "${TMPDIR:-/tmp}/openblas-relaxed-simd.XXXXXX")"
trap 'rm -rf "$OUT"' EXIT

build_variant() {
  local dest="$1"
  shift
  rm -f "${OBJS[@]/#/kernel/}"
  "$MAKE_BIN" -C kernel -j"$JOBS" "${OBJS[@]}" "${MAKE_COMMON[@]}" "$@"
  mkdir -p "$dest"
  local o
  for o in "${OBJS[@]}"; do
    cp "kernel/$o" "$dest/"
  done
}

echo "==> compile default kernels"
build_variant "$OUT/default"

echo "==> compile WASM_RELAXED_SIMD=1 kernels"
build_variant "$OUT/relaxed" WASM_RELAXED_SIMD=1

sum_op() {
  local dir="$1" op="$2" total=0 n
  local o
  for o in "$dir"/*.o; do
    n="$(count_op "$o" "$op")"
    total=$((total + n))
  done
  echo "$total"
}

DEFAULT_MADD="$(sum_op "$OUT/default" 'relaxed_madd')"
RELAXED_MADD="$(sum_op "$OUT/relaxed" 'relaxed_madd')"
DEFAULT_MUL="$(sum_op "$OUT/default" 'f32x4.mul')"
RELAXED_MUL="$(sum_op "$OUT/relaxed" 'f32x4.mul')"

echo "==> opcode counts (wasm-dis)"
echo "  default  relaxed_madd=$DEFAULT_MADD  f32x4.mul=$DEFAULT_MUL"
echo "  relaxed  relaxed_madd=$RELAXED_MADD  f32x4.mul=$RELAXED_MUL"

expect "default objects have no relaxed_madd" "$DEFAULT_MADD" eq 0
expect "WASM_RELAXED_SIMD=1 objects use relaxed_madd" "$RELAXED_MADD" gt 0
expect "default objects still use f32x4.mul" "$DEFAULT_MUL" gt 0

link_wasm() {
  local dir="$1"
  local extra=()
  if [[ "${2:-}" == "relaxed" ]]; then
    extra+=(-mrelaxed-simd)
  fi
  emcc -O2 -msimd128 "${extra[@]}" --no-entry \
    -sEXPORTED_FUNCTIONS="$EXPORTS" \
    -sERROR_ON_UNDEFINED_SYMBOLS=0 \
    -sSTANDALONE_WASM=1 \
    "$dir"/*.o -o "$dir/kernels.wasm"
}

echo "==> link and instantiate"
link_wasm "$OUT/default"
link_wasm "$OUT/relaxed" relaxed

LINKED_DEFAULT_MADD="$(count_op "$OUT/default/kernels.wasm" 'relaxed_madd')"
LINKED_RELAXED_MADD="$(count_op "$OUT/relaxed/kernels.wasm" 'relaxed_madd')"
expect "linked default module has no relaxed_madd" "$LINKED_DEFAULT_MADD" eq 0
expect "linked WASM_RELAXED_SIMD=1 module has relaxed_madd" "$LINKED_RELAXED_MADD" gt 0

if command -v node >/dev/null 2>&1; then
  node -e '
const fs = require("fs");
const defb = fs.readFileSync(process.argv[1]);
const relb = fs.readFileSync(process.argv[2]);
if (!WebAssembly.validate(defb)) {
  console.error("FAIL  default module failed WebAssembly.validate");
  process.exit(1);
}
WebAssembly.instantiate(defb).then(function () {
  console.log("PASS  default module instantiates (SIMD128 only)");
  if (!WebAssembly.validate(relb)) {
    console.log("SKIP  WASM_RELAXED_SIMD=1 module does not validate on this engine (expected without relaxed SIMD)");
    return;
  }
  return WebAssembly.instantiate(relb).then(function () {
    console.log("PASS  WASM_RELAXED_SIMD=1 module instantiates on this engine");
  });
}).catch(function (err) {
  console.error(err);
  process.exit(1);
});
' "$OUT/default/kernels.wasm" "$OUT/relaxed/kernels.wasm"
else
  echo "SKIP  node not found; opcode checks still apply"
fi

# Rebuild default objects so a following 'make libs' is not mixed.
rm -f "${OBJS[@]/#/kernel/}"
"$MAKE_BIN" -C kernel -j"$JOBS" "${OBJS[@]}" "${MAKE_COMMON[@]}"

if [[ "$FAILS" -ne 0 ]]; then
  echo "$FAILS check(s) failed" >&2
  exit 1
fi
echo "relaxed-SIMD opt-in checks passed"
