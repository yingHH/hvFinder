#!/usr/bin/env bash
# =============================================================================
# Build the hvFinder desktop client on macOS or Linux, inside a conda environment.
#
#   bash build/build_unix.sh                  # create/reuse the conda env `hvFinder`
#   bash build/build_unix.sh --refresh-env    # conda env update --prune, then build
#   bash build/build_unix.sh --env-name hvFinder-414
#   PYTHON=/usr/bin/python3 bash build/build_unix.sh    # bypass conda entirely
#
# The environment comes from build/environment-unix.yml (Python 3.9, PySide6,
# paramiko, PyYAML, PyInstaller). PyInstaller cannot cross-compile: this produces a
# binary for the platform it runs on (macOS -> Mach-O, Linux -> ELF). There is no
# .app bundle; the result is a single executable, `dist/hvFinder`.
#
# Output:
#   dist/hvFinder            single-file executable
#   dist/settings.ini        template copied next to it
#   build/BUILD_INFO.txt     size, SHA256 and tool versions
# =============================================================================
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$REPO/build"
DIST_DIR="$REPO/dist"
SPEC="$REPO/hvFinder.spec"
ENV_YML="$BUILD_DIR/environment-unix.yml"
ENTRY="$REPO/client/gui.py"
SETTINGS="$REPO/client/settings.ini"

ENV_NAME="hvFinder"
REFRESH=0
CLEAN=1
PY_OVERRIDE="${PYTHON:-}"

ok()   { printf '  [OK]   %s\n' "$*"; }
warn() { printf '  [WARN] %s\n' "$*"; }
die()  { printf 'ERROR: %b\n' "$*" >&2; exit 1; }
sec()  { printf '\n== %s ==\n' "$*"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --refresh-env) REFRESH=1; shift ;;
        --env-name)    ENV_NAME="${2:-}"; shift 2 ;;
        --no-clean)    CLEAN=0; shift ;;
        -h|--help)     sed -n '3,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)             die "unknown argument: $1 (try --help)" ;;
    esac
done

[[ -n "$ENV_NAME" ]] || die "--env-name must not be empty"

# ---------------------------------------------------------------- 1. platform
sec "platform"
UNAME="$(uname -s)"
case "$UNAME" in
    Darwin) PLATFORM="macOS $(uname -r)" ;;
    Linux)  PLATFORM="Linux $(uname -r)" ;;
    *)      die "unsupported platform: $UNAME (use build/build_windows.ps1 on Windows)" ;;
esac
ok "$PLATFORM"
ok "repository: $REPO"

[[ -f "$SPEC" ]]  || die "missing $SPEC"
[[ -f "$ENTRY" ]] || die "missing $ENTRY (is this the repository root?)"

# ---------------------------------------------------------------- 2. conda
sec "conda environment"
CONDA=""
if [[ -n "$PY_OVERRIDE" ]]; then
    warn "PYTHON=$PY_OVERRIDE given; skipping conda entirely"
    PY="$PY_OVERRIDE"
    command -v "$PY" >/dev/null 2>&1 || die "python not found: $PY"
else
    CONDA="${CONDA_EXE:-}"
    if [[ -z "$CONDA" || ! -x "$CONDA" ]]; then
        CONDA="$(command -v conda 2>/dev/null || true)"
    fi
    [[ -n "$CONDA" ]] || die "conda not found; install Miniconda/Anaconda, or pass PYTHON=... to bypass"
    ok "conda  $CONDA"
    printf '         %s\n' "$("$CONDA" --version 2>&1 | head -1)"
    [[ -f "$ENV_YML" ]] || die "missing $ENV_YML"

    # `conda env list` prints one env per line, name in the first column
    if "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        if [[ "$REFRESH" -eq 1 ]]; then
            printf '  refreshing conda env %s from %s\n' "$ENV_NAME" "$(basename "$ENV_YML")"
            "$CONDA" env update -n "$ENV_NAME" -f "$ENV_YML" --prune \
                || die "conda env update failed"
        else
            ok "reusing conda env $ENV_NAME (use --refresh-env to update it)"
        fi
    else
        printf '  creating conda env %s from %s\n' "$ENV_NAME" "$(basename "$ENV_YML")"
        "$CONDA" env create -n "$ENV_NAME" -f "$ENV_YML" \
            || die "conda env create failed"
        ok "conda env $ENV_NAME created"
    fi
    PY="$CONDA"
fi

# everything below goes through `conda run -n <env>` when conda is in play
run_py() {
    if [[ -n "$CONDA" ]]; then
        "$CONDA" run -n "$ENV_NAME" python "$@"
    else
        "$PY" "$@"
    fi
}

# ---------------------------------------------------------------- 3. python
sec "python"
PY_VER="$(run_py -c 'import sys; print(sys.version.split()[0])')" \
    || die "cannot run python in the environment"
ok "python $PY_VER"
run_py -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)' \
    || die "python >= 3.9 required"

# ---------------------------------------------------------------- 4. deps
sec "dependencies"
MISSING=""
for mod in PySide6 paramiko yaml PyInstaller; do
    if run_py -c "import $mod" >/dev/null 2>&1; then
        ok "$mod"
    else
        warn "$mod missing"
        MISSING="$MISSING $mod"
    fi
done
run_py -c "import cryptography" >/dev/null 2>&1 && ok "cryptography" \
    || warn "cryptography missing (paramiko needs it for key-based auth)"

if [[ -n "$MISSING" ]]; then
    printf '\n'
    if [[ -n "$CONDA" ]]; then
        die "missing modules:$MISSING - refresh the environment with:\n  $CONDA env update -n $ENV_NAME -f $ENV_YML --prune"
    fi
    die "missing modules:$MISSING - install them with:\n  $PY -m pip install -e \"$REPO[client]\" pyinstaller"
fi

# ---------------------------------------------------------------- 5. build
sec "build"
cd "$REPO" || die "cannot cd to $REPO"
# keep intermediates out of build/ (they would sit next to the build scripts)
ARGS=(--noconfirm --workpath "$BUILD_DIR/build" --distpath "$DIST_DIR")
[[ "$CLEAN" -eq 1 ]] && ARGS=(--clean "${ARGS[@]}")
printf '  pyinstaller %s hvFinder.spec\n' "${ARGS[*]}"
run_py -m PyInstaller "${ARGS[@]}" "$SPEC" \
    || die "PyInstaller failed; see the tail of the log above"

# ---------------------------------------------------------------- 6. verify
sec "artifact"
BIN="$DIST_DIR/hvFinder"
[[ -f "$BIN" ]] || die "expected $BIN but it was not produced"
chmod +x "$BIN"
ok "binary: $BIN ($(du -h "$BIN" | cut -f1))"

if [[ -f "$SETTINGS" ]]; then
    cp -f "$SETTINGS" "$DIST_DIR/settings.ini"
    ok "settings.ini copied into dist/"
else
    warn "client/settings.ini not found; a default one is created on first start"
fi

if command -v sha256sum >/dev/null 2>&1; then
    SHA="$(sha256sum "$BIN" | cut -d' ' -f1)"
elif command -v shasum >/dev/null 2>&1; then
    SHA="$(shasum -a 256 "$BIN" | cut -d' ' -f1)"
else
    SHA="unavailable"
fi
SIZE="$(du -h "$BIN" | cut -f1)"
{
    echo "platform   : $PLATFORM"
    echo "built      : $(date -Is)"
    echo "binary     : dist/hvFinder"
    echo "size       : $SIZE"
    echo "sha256     : $SHA"
    if [[ -n "$CONDA" ]]; then
        echo "conda env  : $ENV_NAME"
        echo "conda      : $("$CONDA" --version 2>&1 | head -1)"
    else
        echo "conda env  : (bypassed, PYTHON=$PY_OVERRIDE)"
    fi
    echo "python     : $PY_VER"
    echo "pyinstaller: $(run_py -m PyInstaller --version 2>/dev/null)"
    echo "pyside6    : $(run_py -c 'import PySide6; print(PySide6.__version__)' 2>/dev/null)"
} > "$BUILD_DIR/BUILD_INFO.txt"
ok "build/BUILD_INFO.txt written"

# ---------------------------------------------------------------- 7. next steps
sec "next steps"
printf '  1) run it from the dist directory:  cd dist && ./hvFinder\n'
printf '  2) fill in the SSH connection settings in dist/settings.ini\n'
printf '  3) run one sample through Tab 2 and Tab 3\n'
if [[ "$UNAME" == "Darwin" ]]; then
    printf '\n  macOS: the binary is unsigned, so Gatekeeper may refuse it.\n'
    printf '         xattr -d com.apple.quarantine dist/hvFinder\n'
    printf '         For distribution, codesign + notarise it.\n'
fi
