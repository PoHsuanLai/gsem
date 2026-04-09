#!/bin/bash
# Download reference data for benchmarking.
# Sources: Zenodo mirrors (Broad Institute URLs are dead).
set -e
cd "$(dirname "$0")"

mkdir -p data/eur_w_ld_chr

# ── LD scores ─────────────────────────────────────────────────────────────────
if [ ! -f "data/eur_w_ld_chr/1.l2.ldscore.gz" ]; then
    echo "Downloading LD scores from Zenodo..."
    curl -L -o data/eur_w_ld_chr.tar.gz \
        "https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1"
    cd data && tar xzf eur_w_ld_chr.tar.gz && rm eur_w_ld_chr.tar.gz && cd ..
    echo "LD scores extracted."
else
    echo "LD scores already present."
fi

# ── HapMap3 SNP list ──────────────────────────────────────────────────────────
if [ ! -f "data/w_hm3.snplist" ]; then
    echo "Downloading HapMap3 SNP list from Zenodo..."
    curl -L -o data/w_hm3.snplist.gz \
        "https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1"
    gunzip data/w_hm3.snplist.gz
    echo "HapMap3 SNP list extracted."
else
    echo "HapMap3 SNP list already present."
fi

# ── Simulated GWAS data ──────────────────────────────────────────────────────
if [ ! -f "data/iter1.GWAS1.sumstats.gz" ]; then
    echo "Generating simulated GWAS data with simLDSC..."
    Rscript generate_bench_data.R
else
    echo "Simulated GWAS data already present."
fi

echo ""
echo "Setup complete. Files:"
ls -lh data/eur_w_ld_chr/1.l2.ldscore.gz 2>/dev/null || echo "  (LD scores missing)"
ls -lh data/w_hm3.snplist 2>/dev/null || echo "  (HapMap3 missing)"
ls -lh data/iter1.GWAS*.sumstats.gz 2>/dev/null || echo "  (GWAS data missing)"
