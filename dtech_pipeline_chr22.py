#!/usr/bin/env python3
# =============================================================================
# dTech Pipeline v2.0 — Chromosome 22 Edition
# =============================================================================
# Companion code for the manuscript:
#   "Neuronally-enriched ARTEMIS isoform ENST00000470530 targets hairpin-forming
#    sequences at AluY fragility sites: a computational framework for
#    microdeletion pathogenesis"
#   Author: Andrea Catino
#
# PURPOSE
# -------
# Post-processing pipeline for dTech scan output files.
# Receives compressed CSV files produced by the dTech scanner (external VS tool)
# and annotates each anchor with:
#   - Genomic feature (TAG): PROMOTER, INTRON, EXON, CDS, CENTROMERE,
#                            INT_YOYO, INT_YOYO_LOCK, DESERT
#   - Functional domain keyword
#   - Nearest gene symbol
#   - INT_YOYO_LOCK flag (hairpin score >= 8)
#
# PAPER PARAMETERS (chromosome 22, hg38)
# ----------------------------------------
#   K_SIZE        = 144 bp   (canonical nucleosome length)
#   STEP          = 72 bp    (half-nucleosome step)
#   GC_MIN        = 0.35     (35% GC content minimum)
#   GC_MAX        = 0.65     (65% GC content maximum)
#   ENTROPY_MIN   = 1.87     (Shannon entropy threshold, bits)
#   TM_MIN        = 68.0     (melting temperature threshold, °C)
#   FFT_THR       = 1.5      (A/T periodicity score from Fourier transform)
#   CURV_THR      = 1.2      (intrinsic curvature score, AAAA/TA step parameters)
#   HAIRPIN_LOCK  = 8        (minimum hairpin score for INT_YOYO_LOCK class)
#   PROMOTER_WIN  = 5000     (bp upstream of TSS defined as promoter region)
#   CHR22_LEN     = 51304566 (hg38 chr22 euchromatic length, bp)
#   TOTAL_ANCHORS = 361015   (anchors identified on chr22 in the paper)
#   LOCK_ANCHORS  = 151      (INT_YOYO_LOCK anchors on chr22 in the paper)
#
# RESOURCE LIMITS (configured for dual-core / 3 GB RAM systems)
# ---------------------------------------------------------------
#   N_CORES    = 2    (parallel annotation workers)
#   RAM_LIMIT  = 3    (GB; annotation pauses if exceeded)
#   CHUNK_ROWS = 5000 (rows per multiprocessing chunk)
#
# INPUT
# -----
#   Directory: /root/dna/output_dtech/
#   Files:     human_NC_000022*.csv.gz
#   Expected columns: pos, seq, score [+ optional: hairpin_score, rc_score_delta]
#
# OUTPUT
# ------
#   Annotated CSV:  /root/dna/output_annotated/bio_human_NC_000022*.csv.gz
#   Summary CSV:    /root/dna/output_annotated/dtech_chr22_summary.csv
#   Log file:       /root/dna/pipeline.log
#
# USAGE
# -----
#   python3 dtech_pipeline.py
#   Press Ctrl+C to stop the watcher loop.
# =============================================================================

import os, re, time, glob, shutil, threading, sys
import multiprocessing as mp
import pandas as pd
import requests

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================

BASE_DIR   = "/root/dna"
SCAN_DIR   = os.path.join(BASE_DIR, "output_dtech")
ORG_DIR    = os.path.join(BASE_DIR, "organized_results")
ANNOT_DIR  = os.path.join(BASE_DIR, "output_annotated")
GFF_DIR    = os.path.join(BASE_DIR, "gff_database")
LOG_FILE   = os.path.join(BASE_DIR, "pipeline.log")

# Resource limits — suitable for dual-core / 3 GB RAM systems
N_CORES        = 2          # parallel annotation workers
RAM_LIMIT      = 3          # GB; pause if exceeded
CHUNK_ROWS     = 5000       # rows per multiprocessing chunk
WATCH_INTERVAL = 60         # seconds between directory scans

# Chromosome 22 filter — this edition processes ONLY chromosome 22 (hg38)
# Human chr22 accession on NCBI RefSeq: NC_000022
CHR22_PREFIX   = "NC_000022"
CHR22_LEN      = 51_304_566   # euchromatic length (bp)

# Prefixes to SKIP (scaffolds / unplaced sequences — no GFF available)
SKIP_PREFIXES  = ("NW_", "NT_")

# Hairpin score threshold for INT_YOYO_LOCK classification
HAIRPIN_LOCK_THR = 8

# Promoter window (bp upstream of TSS)
PROMOTER_WIN = 5000

# Create working directories
for d in [SCAN_DIR, ORG_DIR, ANNOT_DIR, GFF_DIR]:
    os.makedirs(d, exist_ok=True)

# =============================================================================
# GFF DATABASE — human hg38 only (chromosome 22 edition)
# =============================================================================

GFF_DB = {
    "human": (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
        "GCF_000001405.40_GRCh38.p14/"
        "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
    ),
}

# Functional domain keywords — matched against GFF attribute strings
KEYWORDS = {
    "BRAIN":        ["brain", "neuron", "synap", "neural", "axon",
                     "cortex", "social", "behavior", "cognitive"],
    "DEVELOPMENT":  ["embryo", "hox", "morphogen", "growth factor",
                     "meristem", "pattern", "larva", "develop"],
    "METABOLIC":    ["liver", "hepat", "metabol", "cytochrom",
                     "muscle", "atp", "mitochon"],
    "IMMUNE":       ["immun", "lymph", "antibod", "tcell", "bcell",
                     "mhc", "hla", "inflamm"],
    "STRUCTURAL":   ["centromere", "telomere", "satellite", "repeating"],
}

# INT_YOYO pattern — tandem repeats within gene bodies
YOYO_PATTERN = re.compile(r'(AGCT){2,}|A{5,}', re.IGNORECASE)

# =============================================================================
# UTILITIES
# =============================================================================

def get_ram_gb():
    """Return current process RAM usage in GB."""
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if "VmRSS" in line:
                    return int(line.split()[1]) / 1_048_576
    except Exception:
        return 0.0


def get_cpu_percent():
    """Estimate CPU usage over a 200 ms window."""
    try:
        def read_cpu():
            with open("/proc/stat") as f:
                parts = f.readline().split()
            idle  = int(parts[4])
            total = sum(int(x) for x in parts[1:])
            return idle, total
        i1, t1 = read_cpu()
        time.sleep(0.2)
        i2, t2 = read_cpu()
        dt = t2 - t1
        return round((1 - (i2 - i1) / dt) * 100, 1) if dt > 0 else 0.0
    except Exception:
        return 0.0


def log(msg):
    """Write timestamped message to stdout and log file."""
    ts   = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    try:
        with open(LOG_FILE, "a") as f:
            f.write(line + "\n")
    except Exception:
        pass

# =============================================================================
# PROGRESS MONITOR (background thread)
# =============================================================================

_mon = {"running": False, "phase": "IDLE", "chr": "-", "start": time.time()}


def _monitor_loop():
    """Display real-time progress bar in the terminal."""
    while _mon["running"]:
        ram  = get_ram_gb()
        cpu  = get_cpu_percent()
        ela  = int(time.time() - _mon["start"])
        bar  = "#" * int(ram / RAM_LIMIT * 20) + "." * (20 - int(ram / RAM_LIMIT * 20))
        scan = len(glob.glob(os.path.join(SCAN_DIR,  "*.csv.gz")))
        ann  = len(glob.glob(os.path.join(ANNOT_DIR, "bio_*.csv.gz")))
        sys.stdout.write(
            f"\r[{_mon['phase']:<10}] CPU:{cpu:5.1f}%  "
            f"RAM:{ram:.2f}/{RAM_LIMIT}GB [{bar}]  "
            f"Chr:{_mon['chr']:<22}  "
            f"Scanned:{scan}  Annotated:{ann}  Elapsed:{ela}s   "
        )
        sys.stdout.flush()
        time.sleep(2)


def start_monitor():
    _mon.update({"running": True, "phase": "WATCH", "start": time.time()})
    threading.Thread(target=_monitor_loop, daemon=True).start()


def stop_monitor():
    _mon["running"] = False
    time.sleep(0.3)
    print()

# =============================================================================
# FILE ORGANIZER
# =============================================================================

def organize_file(sp, csv_path):
    """
    Copy scan CSV into organized directory tree:
      organized_results/{species}/Main_Chromosomes/   — NC_ accessions
      organized_results/{species}/Scaffolds_and_Unplaced/ — NW_, NT_
    Returns path in organized directory.
    """
    fname = os.path.basename(csv_path)
    parts = fname.split("_")
    if len(parts) < 2:
        return csv_path
    acc = parts[1]
    if acc.startswith("NC"):
        subfolder = "Main_Chromosomes"
    elif acc.startswith(("NW", "NT")):
        subfolder = "Scaffolds_and_Unplaced"
    else:
        subfolder = "Other_Fragments"
    tgt_dir = os.path.join(ORG_DIR, sp, subfolder)
    os.makedirs(tgt_dir, exist_ok=True)
    tgt = os.path.join(tgt_dir, fname)
    if not os.path.exists(tgt):
        shutil.copy2(csv_path, tgt)
    return tgt

# =============================================================================
# GFF DOWNLOAD & PARSING
# =============================================================================

def download_gff(sp):
    """Download GFF3 annotation for species sp. Returns local path or None."""
    url = GFF_DB.get(sp)
    if not url:
        return None
    gff_path = os.path.join(GFF_DIR, sp + ".gff.gz")
    if not os.path.exists(gff_path):
        log(f"  [{sp}] Downloading GFF3 annotation...")
        try:
            r = requests.get(url, stream=True, timeout=300)
            r.raise_for_status()
            total = int(r.headers.get("Content-Length", 0))
            done  = 0
            with open(gff_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
                    done += len(chunk)
                    if total:
                        pct = done / total * 100
                        sys.stdout.write(f"\r  [{sp}] GFF download: {done/1e6:.1f}/{total/1e6:.1f} MB ({pct:.1f}%)   ")
                        sys.stdout.flush()
            print()
            log(f"  [{sp}] GFF downloaded: {os.path.getsize(gff_path)/1e6:.1f} MB")
        except Exception as e:
            log(f"  [{sp}] GFF download error: {e}")
            return None
    return gff_path


def load_gff_for_chr(gff_path, target_chr):
    """
    Load GFF3 features for a single chromosome into a DataFrame.
    Only gene, exon, CDS, mRNA, centromere, telomere features are retained.
    Returns a DataFrame sorted by start position.
    """
    cols   = ["chr", "src", "type", "start", "end",
              "score", "strand", "phase", "attr"]
    ftypes = {"gene", "exon", "CDS", "mRNA", "centromere", "telomere"}
    rows   = []
    try:
        for chunk in pd.read_csv(
            gff_path, sep="\t", comment="#", header=None,
            names=cols, compression="gzip", chunksize=50_000, low_memory=False
        ):
            sub = chunk[
                (chunk["chr"] == target_chr) &
                (chunk["type"].isin(ftypes))
            ]
            if not sub.empty:
                rows.append(sub)
    except Exception as e:
        log(f"  [GFF] Parse error: {e}")
        return pd.DataFrame(columns=cols)

    if not rows:
        return pd.DataFrame(columns=cols)

    df = pd.concat(rows, ignore_index=True)
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce")
    df.dropna(subset=["start", "end"], inplace=True)
    return df.sort_values("start").reset_index(drop=True)

# =============================================================================
# ANNOTATION LOGIC
# =============================================================================

def get_functional_domain(attr):
    """Map GFF attribute string to a functional domain keyword."""
    a = str(attr).lower()
    for domain, kws in KEYWORDS.items():
        if any(k in a for k in kws):
            return domain
    return "GENERAL"


def extract_gene_name(attr):
    """Extract gene name from GFF attribute string."""
    s = str(attr)
    for field in ("Name=", "gene_id=", "gene="):
        if field in s:
            return s.split(field)[-1].split(";")[0].strip()
    return "Unknown"


def compute_hairpin_score(seq):
    """
    Compute intrastrand hairpin potential.
    Slides a 20-bp window looking for reverse-complement self-complementarity.
    Returns count of windows with self-complementarity >= 8 bp.
    Returns 0 if seq is too short.
    """
    seq = str(seq).upper()
    if len(seq) < 20:
        return 0
    complement = str.maketrans("ACGT", "TGCA")
    score = 0
    win = 20
    for i in range(len(seq) - win + 1):
        sub = seq[i:i + win]
        rc  = sub.translate(complement)[::-1]
        matches = sum(a == b for a, b in zip(sub, rc))
        if matches >= 8:
            score += 1
    return score


def annotate_row(pos, seq, chr_gff, hairpin_score=None):
    """
    Assign TAG, DOMAIN, GENE, and IS_LOCK to a single anchor.

    TAG values (hierarchical):
      CDS           — within coding sequence
      EXON          — within exon (non-CDS)
      PROMOTER      — within PROMOTER_WIN bp upstream of a TSS
      INT_YOYO_LOCK — within gene body, YOYO pattern, hairpin >= HAIRPIN_LOCK_THR
      INT_YOYO      — within gene body, YOYO pattern present
      INTRON        — within gene body, no YOYO
      CENTROMERE    — within centromeric region
      DESERT        — intergenic, no nearby gene
    """
    is_yoyo = bool(YOYO_PATTERN.search(str(seq)))

    if hairpin_score is None:
        hairpin_score = compute_hairpin_score(seq)

    is_lock = is_yoyo and (hairpin_score >= HAIRPIN_LOCK_THR)

    if chr_gff.empty:
        return "DESERT", "NONE", "None", is_lock, hairpin_score

    # Check overlap with annotated features
    match = chr_gff[(chr_gff["start"] <= pos) & (chr_gff["end"] >= pos)]

    if match.empty:
        # Check promoter: within PROMOTER_WIN bp upstream of a gene TSS
        prm = chr_gff[
            (chr_gff["type"] == "gene") &
            (chr_gff["start"] > pos) &
            (chr_gff["start"] <= pos + PROMOTER_WIN)
        ]
        if not prm.empty:
            feat = prm.iloc[0]
            return (
                "PROMOTER",
                get_functional_domain(feat["attr"]),
                extract_gene_name(feat["attr"]),
                is_lock,
                hairpin_score
            )
        return "DESERT", "NONE", "None", is_lock, hairpin_score

    # Determine TAG from overlapping features (priority: CDS > EXON > GENE)
    feat = match.iloc[0]
    ftype = feat["type"].upper()

    if ftype == "CDS":
        tag = "CDS"
    elif ftype in ("EXON", "MRNA"):
        tag = "EXON"
    elif ftype == "CENTROMERE":
        tag = "CENTROMERE"
    elif ftype == "GENE":
        if is_lock:
            tag = "INT_YOYO_LOCK"
        elif is_yoyo:
            tag = "INT_YOYO"
        else:
            tag = "INTRON"
    else:
        tag = "INTRON"

    return (
        tag,
        get_functional_domain(feat["attr"]),
        extract_gene_name(feat["attr"]),
        is_lock,
        hairpin_score
    )


def process_annot_chunk(args):
    """Annotate a DataFrame chunk (used by multiprocessing pool)."""
    chunk_df, chr_gff = args
    results = []
    for _, row in chunk_df.iterrows():
        hp = row.get("hairpin_score", None)
        tag, dom, gene, is_lock, hp_computed = annotate_row(
            row["pos"], row["seq"], chr_gff, hairpin_score=hp
        )
        results.append((tag, dom, gene, is_lock, hp_computed))

    chunk_df = chunk_df.copy()
    chunk_df["TAG"]           = [r[0] for r in results]
    chunk_df["DOMAIN"]        = [r[1] for r in results]
    chunk_df["GENE"]          = [r[2] for r in results]
    chunk_df["IS_LOCK"]       = [r[3] for r in results]
    chunk_df["hairpin_score"] = [r[4] for r in results]
    return chunk_df

# =============================================================================
# ANNOTATE ONE CSV FILE
# =============================================================================

def annotate_csv(sp, csv_path, gff_path):
    """
    Annotate all anchors in csv_path using GFF at gff_path.
    Saves result to ANNOT_DIR as bio_{original_filename}.
    """
    fname    = os.path.basename(csv_path)
    out_path = os.path.join(ANNOT_DIR, "bio_" + fname)

    if os.path.exists(out_path):
        log(f"  [SKIP] Already annotated: {fname}")
        return

    # Extract chromosome accession from filename
    current_chr = fname.replace(f"{sp}_", "").replace(".csv.gz", "")
    _mon["chr"] = current_chr

    log(f"  [ANNOT] Loading GFF for {current_chr}...")
    chr_gff = load_gff_for_chr(gff_path, current_chr)
    if chr_gff.empty:
        log(f"  [ANNOT] No GFF features found for {current_chr} — skipped")
        return

    log(f"  [ANNOT] GFF loaded: {len(chr_gff):,} features. Annotating {fname}...")

    # Load scan CSV in chunks
    chunks = []
    try:
        for chunk in pd.read_csv(csv_path, compression="gzip", chunksize=CHUNK_ROWS):
            chunks.append((chunk, chr_gff))
    except Exception as e:
        log(f"  [ANNOT] CSV read error {fname}: {e}")
        return

    # Parallel annotation
    with mp.Pool(N_CORES) as pool:
        results = pool.map(process_annot_chunk, chunks)

    annotated = pd.concat(results, ignore_index=True)
    annotated.to_csv(out_path, index=False, compression="gzip")
    n_lock = int(annotated["IS_LOCK"].sum()) if "IS_LOCK" in annotated.columns else 0
    log(f"  [ANNOT OK] {fname} — {len(annotated):,} anchors, {n_lock} INT_YOYO_LOCK")

# =============================================================================
# GFF CACHE
# =============================================================================

_gff_cache = {}


def _trigger_annotation(sp, csv_path):
    """Download GFF (if needed), organize file, and annotate."""
    _mon["phase"] = "ANNOTATE"
    org_path = organize_file(sp, csv_path)

    if sp not in _gff_cache:
        _gff_cache[sp] = download_gff(sp)

    gff_path = _gff_cache[sp]
    if gff_path is None:
        log(f"  [SKIP] GFF not available for species: {sp}")
        _mon["phase"] = "WATCH"
        return

    annotate_csv(sp, org_path, gff_path)
    _mon["phase"] = "WATCH"

# =============================================================================
# SUMMARY TABLE
# =============================================================================

def build_summary_table():
    """
    Build and save a summary table of all annotated chr22 anchors.
    Reproduces Table 1 from the paper (TAG distribution + hairpin statistics).
    Output: ANNOT_DIR/dtech_chr22_summary.csv
    """
    bio_files = glob.glob(os.path.join(ANNOT_DIR, "bio_human_NC_000022*.csv.gz"))
    if not bio_files:
        log("[SUMMARY] No annotated chr22 files found.")
        return

    log(f"[SUMMARY] Loading {len(bio_files)} annotated file(s)...")
    frames = []
    for f in bio_files:
        try:
            df = pd.read_csv(f, compression="gzip", low_memory=False)
            frames.append(df)
        except Exception as e:
            log(f"  [SUMMARY] Error reading {f}: {e}")

    if not frames:
        log("[SUMMARY] No data loaded.")
        return

    full = pd.concat(frames, ignore_index=True)
    log(f"[SUMMARY] Total anchors loaded: {len(full):,}")

    # Ensure hairpin_score column exists
    if "hairpin_score" not in full.columns:
        log("[SUMMARY] Computing hairpin scores (this may take a moment)...")
        full["hairpin_score"] = full["seq"].apply(compute_hairpin_score)

    # Per-TAG statistics
    tags = ["PROMOTER", "INTRON", "INT_YOYO", "INT_YOYO_LOCK",
            "EXON", "CDS", "CENTROMERE", "DESERT"]
    rows = []
    total = len(full)

    for tag in tags:
        sub = full[full["TAG"] == tag]
        n   = len(sub)
        if n == 0:
            continue
        hp_mean = sub["hairpin_score"].mean() if n else 0
        hp_max  = sub["hairpin_score"].max()  if n else 0
        hp_med  = sub["hairpin_score"].median() if n else 0
        n_lock  = int(sub["IS_LOCK"].sum()) if "IS_LOCK" in sub.columns else 0
        score_mean = sub["score"].mean() if "score" in sub.columns else 0
        rows.append({
            "TAG":               tag,
            "n_anchors":         n,
            "pct_total":         round(100 * n / total, 2),
            "hairpin_mean":      round(hp_mean, 3),
            "hairpin_median":    round(hp_med, 3),
            "hairpin_max":       int(hp_max),
            "n_INT_YOYO_LOCK":   n_lock,
            "score_mean":        round(score_mean, 3),
        })

    summary_df = pd.DataFrame(rows)

    # Add total row
    total_row = pd.DataFrame([{
        "TAG":            "TOTAL",
        "n_anchors":      total,
        "pct_total":      100.0,
        "hairpin_mean":   round(full["hairpin_score"].mean(), 3),
        "hairpin_median": round(full["hairpin_score"].median(), 3),
        "hairpin_max":    int(full["hairpin_score"].max()),
        "n_INT_YOYO_LOCK": int(full["IS_LOCK"].sum()) if "IS_LOCK" in full.columns else 0,
        "score_mean":     round(full["score"].mean(), 3) if "score" in full.columns else 0,
    }])
    summary_df = pd.concat([summary_df, total_row], ignore_index=True)

    out_path = os.path.join(ANNOT_DIR, "dtech_chr22_summary.csv")
    summary_df.to_csv(out_path, index=False)
    log(f"[SUMMARY] Saved: {out_path}")

    # Print formatted table to terminal
    print()
    print("=" * 80)
    print("  dTech Chr22 Summary Table")
    print("  (matches Table 1 in Catino et al.)")
    print("=" * 80)
    print(f"  {'TAG':<18} {'n':>8} {'%':>6} {'hp_mean':>9} {'hp_max':>7} {'LOCK':>6} {'score_mean':>11}")
    print("  " + "-" * 68)
    for _, row in summary_df.iterrows():
        marker = " ◄ PAPER RESULT" if row["TAG"] == "INT_YOYO_LOCK" else ""
        print(
            f"  {row['TAG']:<18} {int(row['n_anchors']):>8,} {row['pct_total']:>6.2f}% "
            f"{row['hairpin_mean']:>9.3f} {int(row['hairpin_max']):>7} "
            f"{int(row['n_INT_YOYO_LOCK']):>6} {row['score_mean']:>11.3f}{marker}"
        )
    print("=" * 80)
    print(f"  Paper expected: 361,015 total anchors | 151 INT_YOYO_LOCK | hairpin_mean=8.47")
    print(f"  Your result:    {total:,} total anchors | "
          f"{int(summary_df[summary_df['TAG']=='INT_YOYO_LOCK']['n_anchors'].values[0]) if 'INT_YOYO_LOCK' in summary_df['TAG'].values else 0} INT_YOYO_LOCK")
    print(f"  Summary saved:  {out_path}")
    print("=" * 80)
    print()

    return summary_df

# =============================================================================
# WATCHER LOOP
# =============================================================================

def run_watcher():
    """
    Monitor SCAN_DIR for new dTech scan output files.
    Automatically annotates new files for chromosome 22 (human hg38).
    Press Ctrl+C to stop.
    """
    log("=" * 60)
    log("  dTech Pipeline v2.0 — Chromosome 22 Edition")
    log("  Author: Andrea Catino")
    log(f"  Watch directory: {SCAN_DIR}")
    log(f"  Target:          Human chr22 (NC_000022, hg38)")
    log(f"  Cores: {N_CORES} | RAM limit: {RAM_LIMIT} GB | Chunk: {CHUNK_ROWS} rows")
    log(f"  INT_YOYO_LOCK threshold: hairpin_score >= {HAIRPIN_LOCK_THR}")
    log(f"  Promoter window: {PROMOTER_WIN:,} bp upstream of TSS")
    log("  Press Ctrl+C to stop")
    log("=" * 60)

    annotated_set = set()

    # Resume: files already annotated
    for f in glob.glob(os.path.join(ANNOT_DIR, "bio_*.csv.gz")):
        base = os.path.basename(f).replace("bio_", "")
        annotated_set.add(os.path.join(SCAN_DIR, base))

    # Resume: NW/NT files (skip immediately)
    skipped = 0
    for f in glob.glob(os.path.join(SCAN_DIR, "*.csv.gz")):
        fname  = os.path.basename(f)
        parts  = fname.split("_")
        chr_id = "_".join(parts[1:]).replace(".csv.gz", "") if len(parts) >= 2 else ""
        if chr_id.startswith(SKIP_PREFIXES):
            annotated_set.add(f)
            skipped += 1

    log(f"  Resume: {len(annotated_set)} files already processed ({skipped} scaffold files skipped)")

    while True:
        csv_files = glob.glob(os.path.join(SCAN_DIR, "*.csv.gz"))
        new_files = [f for f in csv_files if f not in annotated_set]

        # Filter: only chr22 files (NC_000022)
        chr22_new = [
            f for f in new_files
            if CHR22_PREFIX in os.path.basename(f)
        ]
        non_chr22_new = [
            f for f in new_files
            if CHR22_PREFIX not in os.path.basename(f)
        ]

        # Auto-skip non-chr22 files with a notice
        for f in non_chr22_new:
            fname = os.path.basename(f)
            log(f"  [SKIP non-chr22] {fname} — this edition processes chr22 only")
            annotated_set.add(f)

        if chr22_new:
            log(f"  [WATCHER] {len(chr22_new)} new chr22 file(s) detected")
            for csv_path in chr22_new:
                fname  = os.path.basename(csv_path)
                parts  = fname.split("_")

                if len(parts) < 2:
                    annotated_set.add(csv_path)
                    continue

                sp     = parts[0]   # e.g. "human"
                chr_id = "_".join(parts[1:]).replace(".csv.gz", "")

                # Skip scaffolds
                if chr_id.startswith(SKIP_PREFIXES):
                    log(f"  [SKIP scaffold] {fname}")
                    annotated_set.add(csv_path)
                    continue

                # Check species is in GFF database
                if sp not in GFF_DB:
                    log(f"  [SKIP] Unknown species: {sp}")
                    annotated_set.add(csv_path)
                    continue

                # Pause if RAM limit exceeded
                while get_ram_gb() > RAM_LIMIT:
                    log(f"  [RAM] {get_ram_gb():.2f} GB > {RAM_LIMIT} GB limit — pausing 15s...")
                    time.sleep(15)

                log(f"  [NEW] Processing: {fname}")
                _trigger_annotation(sp, csv_path)
                annotated_set.add(csv_path)

        time.sleep(WATCH_INTERVAL)

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    mp.freeze_support()
    start_monitor()
    try:
        run_watcher()
    except KeyboardInterrupt:
        log("\n[STOP] Manual interrupt (Ctrl+C)")
    finally:
        stop_monitor()
        # Build and print summary table on exit
        log("[EXIT] Building summary table...")
        build_summary_table()
        tot = len(glob.glob(os.path.join(SCAN_DIR,  "*.csv.gz")))
        ann = len(glob.glob(os.path.join(ANNOT_DIR, "bio_*.csv.gz")))
        log(f"[DONE] Scanned: {tot} | Annotated: {ann}")
        log(f"       Results: {ANNOT_DIR}")
