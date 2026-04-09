#!/bin/bash
# Run R and Rust GenomicSEM pipelines and compare outputs.
set -e
cd "$(dirname "$0")"

TRAITS="data/iter1GWAS1.sumstats.gz data/iter1GWAS2.sumstats.gz data/iter1GWAS3.sumstats.gz"
LD="data/eur_w_ld_chr/"
MODEL="F1 =~ NA*V1 + V2 + V3
F1 ~~ 1*F1
V1 ~~ V1
V2 ~~ V2
V3 ~~ V3"

# Build Rust binary
echo "=== Building Rust binary (release) ==="
cargo build --release -p gsem --quiet 2>/dev/null || cargo build --release -p gsem
RUST_BIN="../target/release/genomicsem"

if [ ! -f "$RUST_BIN" ]; then
    echo "ERROR: Rust binary not found at $RUST_BIN"
    exit 1
fi

mkdir -p out_r out_rust

# ── R Pipeline ────────────────────────────────────────────────────────────────
echo ""
echo "=== Running R Pipeline ==="
Rscript r_pipeline.R

# ── Rust Pipeline ─────────────────────────────────────────────────────────────
echo ""
echo "=== Running Rust Pipeline ==="

# LDSC
echo "--- Rust: LDSC ---"
LDSC_START=$(python3 -c "import time; print(time.time())")
RUST_LOG=info $RUST_BIN ldsc \
    --traits $TRAITS \
    --sample-prev "NA,NA,NA" \
    --pop-prev "NA,NA,NA" \
    --ld "$LD" \
    --n-blocks 200 \
    --out out_rust/ldsc.json 2>&1 | head -20
LDSC_END=$(python3 -c "import time; print(time.time())")
LDSC_TIME=$(python3 -c "print(f'{$LDSC_END - $LDSC_START:.2f}')")
echo "  Rust LDSC time: ${LDSC_TIME}s"

# SEM
echo "--- Rust: SEM ---"
SEM_START=$(python3 -c "import time; print(time.time())")
$RUST_BIN sem \
    --covstruc out_rust/ldsc.json \
    --model "$MODEL" \
    --estimation DWLS \
    --out out_rust/sem.tsv 2>&1 | head -20
SEM_END=$(python3 -c "import time; print(time.time())")
SEM_TIME=$(python3 -c "print(f'{$SEM_END - $SEM_START:.2f}')")
echo "  Rust SEM time: ${SEM_TIME}s"

# Save Rust timings
python3 -c "
import json
timings = {'ldsc': float('$LDSC_TIME'), 'sem': float('$SEM_TIME')}
with open('out_rust/timings.json', 'w') as f:
    json.dump(timings, f)
"

# ── Compare ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Comparing Outputs ==="
Rscript compare_outputs.R

echo ""
echo "=== Done ==="
cat results.md
