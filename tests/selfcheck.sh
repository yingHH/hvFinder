#!/usr/bin/env bash
# =============================================================================
# hvFinder server-side self-check (preflight).
#
# Runs **in the foreground** and is **read-only**: it never submits a job, never
# writes into the output directory and never touches reference databases beyond
# `ls` / `blastdbcmd -info`. Use it after installation and before the first real
# sample to find out what is still missing.
#
# Usage:
#   bash tests/selfcheck.sh                          # check everything
#   bash tests/selfcheck.sh --stage preprocess       # only the mNGS stage
#   bash tests/selfcheck.sh --config server/configs/config2.yaml --outdir /data/run
#   bash tests/selfcheck.sh --min-free-gb 200 --quiet
#
# Exit code: 0 = no hard failure (WARN allowed), 1 = at least one [FAIL].
# =============================================================================
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SERVER_DIR="$ROOT/server"
PYTHON="${PYTHON:-python3}"

CONFIG=""
OUTDIR=""
DATA_DIR=""
STAGE="all"
MIN_FREE_GB=20
QUIET=0
DEEP=0

FAILS=0
WARNS=0

ok()   { [[ "$QUIET" -eq 1 ]] || printf '  [OK]   %s\n' "$*"; }
warn() { printf '  [WARN] %s\n' "$*"; WARNS=$((WARNS + 1)); }
bad()  { printf '  [FAIL] %s\n' "$*"; FAILS=$((FAILS + 1)); }
sec()  { [[ "$QUIET" -eq 1 ]] || printf '\n== %s ==\n' "$*"; }
die()  { echo "ERROR: $*" >&2; exit 2; }

usage() { sed -n '3,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)      CONFIG="${2:-}"; shift 2 ;;
        --outdir)      OUTDIR="${2:-}"; shift 2 ;;
        --data-dir)    DATA_DIR="${2:-}"; shift 2 ;;
        --stage)       STAGE="${2:-}"; shift 2 ;;
        --min-free-gb) MIN_FREE_GB="${2:-20}"; shift 2 ;;
        --deep)        DEEP=1; shift ;;
        --quiet)       QUIET=1; shift ;;
        -h|--help)     usage ;;
        *)             die "unknown argument: $1 (try --help)" ;;
    esac
done

case "$STAGE" in
    all|preprocess|detect) ;;
    *) die "--stage must be one of: all, preprocess, detect" ;;
esac

[[ -n "$CONFIG" ]] || CONFIG="$SERVER_DIR/configs/config2.yaml"
[[ -n "$DATA_DIR" ]] || DATA_DIR="$ROOT/data"

# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

# cfg_value <section> [dotted key] -- read a value from config2.yaml.
# Uses PyYAML when available, otherwise a plain grep fallback (a missing PyYAML
# is itself reported by the python section, so the fallback keeps the rest usable).
cfg_value() {
    local section="$1" key="${2:-}"
    [[ -f "$CONFIG" ]] || return 0
    if "$PYTHON" - "$CONFIG" "$section" "$key" <<'PY' 2>/dev/null
import sys
try:
    import yaml
except Exception:
    sys.exit(9)
path, section, key = sys.argv[1], sys.argv[2], sys.argv[3]
try:
    cfg = yaml.safe_load(open(path, encoding="utf-8")) or {}
except Exception:
    sys.exit(9)
node = cfg.get(section) or {}
for part in [p for p in key.split(".") if p]:
    node = (node or {}).get(part) if isinstance(node, dict) else None
print("" if node is None else node)
PY
    then
        return 0
    fi
    # fallback: only handles tools' `bin:` and one-level keys
    if [[ -z "$key" || "$key" == "bin" ]]; then
        awk -v s="$section" '
            $0 ~ "^"s":" {ins=1; next}
            ins && /^[^[:space:]]/ {ins=0}
            ins && $1 == "bin:" {gsub(/"/,"",$2); print $2; exit}
        ' "$CONFIG"
    fi
}

check_bin() {   # <label> <command-or-path>
    local label="$1" target="${2:-}"
    if [[ -z "$target" ]]; then
        bad "$label: not configured (set $label.bin in config2.yaml)"
        return
    fi
    if [[ "$target" == */* ]]; then
        if [[ -x "$target" ]]; then ok "$label  $target"
        else bad "$label not executable: $target"; fi
    elif command -v "$target" >/dev/null 2>&1; then
        ok "$label  $(command -v "$target")"
    else
        bad "$label not found in PATH: $target"
    fi
}

tool_bin() {    # <section> <fallback-name>
    local v; v="$(cfg_value "$1" bin)"
    printf '%s' "${v:-$2}"
}

n_lines() {     # number of non-empty, non-comment lines
    awk 'NF && $1 !~ /^#/' "$1" 2>/dev/null | wc -l | tr -d ' '
}

# --------------------------------------------------------------------------- #
# banner
# --------------------------------------------------------------------------- #
[[ "$QUIET" -eq 1 ]] || {
    echo "hvFinder self-check"
    echo "  time   : $(date -Is)"
    echo "  host   : $(uname -sm)"
    echo "  root   : $ROOT"
    echo "  config : $CONFIG"
    echo "  stage  : $STAGE"
}

# --------------------------------------------------------------------------- #
# 1. python + code
# --------------------------------------------------------------------------- #
sec "python and code"
if command -v "$PYTHON" >/dev/null 2>&1; then
    ok "$("$PYTHON" -V 2>&1)"
    if "$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)' 2>/dev/null; then
        ok "python >= 3.9"
    else
        bad "python >= 3.9 required"
    fi
    if "$PYTHON" -c 'import yaml' >/dev/null 2>&1; then
        ok "pyyaml importable"
    else
        bad "pyyaml missing (pip install pyyaml)"
    fi
    # syntax only: ast.parse never writes __pycache__
    for f in mNGS_pip.py pathogen_detect.py logger.py get_script_base.py; do
        if [[ ! -f "$SERVER_DIR/$f" ]]; then
            bad "missing $SERVER_DIR/$f"
        elif "$PYTHON" -c 'import ast,sys; ast.parse(open(sys.argv[1],encoding="utf-8").read())' \
                "$SERVER_DIR/$f" >/dev/null 2>&1; then
            ok "$f parses"
        else
            bad "$f has a syntax error"
        fi
    done
else
    bad "python not found: $PYTHON"
fi

if [[ -f "$CONFIG" ]]; then
    ok "config2.yaml present"
else
    bad "config not found: $CONFIG (pass --config)"
fi

# --------------------------------------------------------------------------- #
# 2. external tools
# --------------------------------------------------------------------------- #
sec "external tools"
check_bin fastp           "$(tool_bin fastp fastp)"
check_bin bowtie2         "$(tool_bin bowtie2 bowtie2)"
check_bin samtools        "$(tool_bin samtools samtools)"
check_bin megahit         "$(tool_bin megahit megahit)"
check_bin cdhit           "$(tool_bin cdhit cd-hit-est)"

if [[ "$STAGE" != "detect" ]]; then
    rb="$(tool_bin ribodetector ribodetector_cpu)"
    check_bin ribodetector "$rb"
    # ribodetector also needs its model files; only a WARN, they ship with the tool
    check_bin bowtie2-build "$(tool_bin bowtie2 bowtie2)-build"
fi

if [[ "$STAGE" != "preprocess" ]]; then
    check_bin blastn   "$(tool_bin blastn blastn)"
    check_bin bedtools "$(tool_bin bedtools bedtools)"
    check_bin bedGraphToBigWig "$(tool_bin bedgraph_to_bigwig bedGraphToBigWig)"
    bc="$(tool_bin blastdbcmd blastdbcmd)"
    if [[ "$bc" == */* ]]; then
        [[ -x "$bc" ]] && ok "blastdbcmd (optional)  $bc" \
            || warn "blastdbcmd not executable: $bc (only genome length fallback needs it)"
    elif command -v "$bc" >/dev/null 2>&1; then
        ok "blastdbcmd (optional)  $(command -v "$bc")"
    else
        warn "blastdbcmd missing (optional; genome breadth comes from BLASTN slen)"
    fi
fi

check_bin diamond "$(tool_bin diamond diamond)"

# --------------------------------------------------------------------------- #
# 3. reference data
# --------------------------------------------------------------------------- #
sec "reference data"

# ---- target taxids
taxids="$(cfg_value databases virus_taxids)"
[[ -n "$taxids" ]] || taxids="$DATA_DIR/virus_taxids_high_pathogenic.txt"
[[ "$taxids" == /* ]] || taxids="$ROOT/$taxids"
if [[ -f "$taxids" ]]; then
    n=$(n_lines "$taxids")
    if [[ "$n" -gt 0 ]]; then ok "taxid list  $taxids  ($n entries)"
    else bad "taxid list has no entries: $taxids"; fi
else
    bad "taxid list not found: $taxids"
fi

# ---- BLASTN nucleotide database
vindex="$(cfg_value databases blastn_index)"
[[ -n "$vindex" ]] || vindex="$DATA_DIR/hvn/hvn.fa"
[[ "$vindex" == /* ]] || vindex="$ROOT/$vindex"
if [[ -f "$vindex" ]] || ls "${vindex}".n* >/dev/null 2>&1; then
    ok "BLASTN database  $vindex"
    if command -v "$(tool_bin blastdbcmd blastdbcmd)" >/dev/null 2>&1; then
        info="$("$(tool_bin blastdbcmd blastdbcmd)" -db "$vindex" -info 2>/dev/null | head -3)"
        [[ -n "$info" ]] && ok "blastdbcmd -info OK ($(printf '%s' "$info" | tr '\n' ';'))" \
            || warn "blastdbcmd -info returned nothing (is the DB made with makeblastdb?)"
    fi
    [[ "$DEEP" -eq 1 ]] || true
else
    bad "BLASTN database not found: $vindex (run makeblastdb - see data/README.md)"
fi

# ---- bowtie2 host index
hindex="$(cfg_value bowtie2 params.index_path)"
[[ -n "$hindex" ]] || hindex="host_idx/host_idx"
[[ "$hindex" == /* ]] || hindex="$ROOT/$hindex"
if compgen -G "${hindex}*.bt2" >/dev/null || compgen -G "${hindex}.1.bt2" >/dev/null; then
    ok "bowtie2 host index  $hindex"
else
    bad "bowtie2 host index not found: ${hindex}*.bt2 (run bowtie2-build)"
fi

# ---- DIAMOND databases
check_diamond_db() {   # <label> <section> <key> <default>
    local label="$1" section="$2" key="$3" def="$4" p size
    p="$(cfg_value "$section" "$key")"
    [[ -n "$p" ]] || p="$def"
    if [[ "$p" != /* ]]; then
        if [[ -e "$ROOT/$p" || ! -e "$DATA_DIR/${p#data/}" ]]; then p="$ROOT/$p"
        else p="$DATA_DIR/${p#data/}"; fi
    fi
    if [[ ! -f "$p" ]]; then
        bad "$label not found: $p (see data/README.md)"
        return
    fi
    size=$(du -h "$p" 2>/dev/null | cut -f1)
    ok "$label  $p  (${size:-?})"
    if [[ "$DEEP" -eq 1 ]] && command -v "$(tool_bin diamond diamond)" >/dev/null 2>&1; then
        dbinfo="$("$(tool_bin diamond diamond)" dbinfo --db "$p" 2>/dev/null | head -6)"
        [[ -n "$dbinfo" ]] && ok "diamond dbinfo OK" || warn "diamond dbinfo failed for $p"
    fi
}

if [[ "$STAGE" != "detect" ]]; then
    check_diamond_db "hvp protein DB (hvp.dmnd)" prescreen db_path "data/hvp/hvp.dmnd"
fi
if [[ "$STAGE" != "preprocess" ]]; then
    check_diamond_db "nr DB (nr.dmnd)" diamond db_path "data/databases/nr.dmnd"
fi

# placeholder guard: published config ships relative placeholders
if [[ -f "$CONFIG" ]]; then
    # A value counts as a placeholder when it does not start with "/"
    placeholders=$(awk '
        /^[[:space:]]*(index_path|db_path|virus_taxids|blastn_index):/ {
            v = $0
            sub(/^[^:]*:[[:space:]]*/, "", v)
            gsub(/"/, "", v)
            sub(/[[:space:]]*#.*$/, "", v)
            if (v != "" && v !~ /^\//) print v
        }' "$CONFIG" 2>/dev/null)
    if [[ -n "$placeholders" ]]; then
        n_ph=$(printf '%s\n' "$placeholders" | wc -l | tr -d ' ')
        shown=$(printf '%s\n' "$placeholders" | head -3 | paste -sd, - | sed 's/,/, /g')
        [[ "$n_ph" -gt 3 ]] && shown="$shown (+$((n_ph - 3)) more)"
        warn "config still has non-absolute paths (shipped placeholders?): $shown"
    fi
fi

# --------------------------------------------------------------------------- #
# 4. output directory
# --------------------------------------------------------------------------- #
sec "output directory"
if [[ -z "$OUTDIR" ]]; then
    [[ "$QUIET" -eq 1 ]] || printf '  [SKIP] no --outdir given\n'
else
    if [[ -d "$OUTDIR" ]]; then
        if [[ -w "$OUTDIR" ]]; then ok "writable  $OUTDIR"
        else bad "not writable: $OUTDIR"; fi
    else
        parent="$(dirname "$OUTDIR")"
        if [[ -d "$parent" && -w "$parent" ]]; then
            ok "creatable  $OUTDIR (parent is writable)"
        else
            bad "cannot create output directory: $OUTDIR"
        fi
    fi
    target="$OUTDIR"; [[ -d "$target" ]] || target="$(dirname "$OUTDIR")"
    avail=$(df -Pk "$target" 2>/dev/null | awk 'NR==2 {print int($4/1048576)}')
    if [[ -n "${avail:-}" ]]; then
        if [[ "$avail" -ge "$MIN_FREE_GB" ]]; then
            ok "free space  ${avail} GB (>= ${MIN_FREE_GB} GB)"
        else
            warn "only ${avail} GB free on $target (< ${MIN_FREE_GB} GB)"
        fi
    fi
fi

# --------------------------------------------------------------------------- #
# 5. optional: scheduler (only relevant for the GUI's Slurm mode)
# --------------------------------------------------------------------------- #
sec "scheduler (optional)"
if command -v sbatch >/dev/null 2>&1; then
    ok "sbatch present ($(sbatch --version 2>/dev/null | head -1))"
else
    warn "sbatch not found - the GUI must use 'Direct Python' mode on this host"
fi

# --------------------------------------------------------------------------- #
# summary
# --------------------------------------------------------------------------- #
if [[ "$QUIET" -eq 1 ]]; then
    printf 'RESULT: %s (%d fail, %d warn)\n' \
        "$([[ "$FAILS" -eq 0 ]] && echo PASS || echo FAIL)" "$FAILS" "$WARNS"
else
    printf '\n== summary ==\n'
    printf '  fail : %d\n  warn : %d\n' "$FAILS" "$WARNS"
    if [[ "$FAILS" -eq 0 ]]; then
        echo "RESULT: PASS"
        echo "  environment looks ready; WARN items are optional."
    else
        echo "RESULT: FAIL"
        echo "  fix the [FAIL] items above (see README 'Required data' / data/README.md)."
    fi
fi
[[ "$FAILS" -eq 0 ]] && exit 0 || exit 1
