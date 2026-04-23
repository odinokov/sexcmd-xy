#!/usr/bin/env python3
"""Human XY sex inference from FASTQ using SEXCMD marker FASTA.

Port of SEXCMD.R (XY + human only), with pigz / awk / bwa / sambamba pipeline.
Upstream reference: https://github.com/lovemun/SEXCMD

Requirements (all on PATH):
    bwa, sambamba, pigz, awk, bash, head

Quick start:
    # Fetch the hg38 marker FASTA (one-time)
    wget -O sex_marker.hg38.filtered.final.fasta \\
        https://raw.githubusercontent.com/lovemun/SEXCMD/master/Examples/Human/sex_marker.hg38.filtered.final.fasta

    # Optional: pre-build the bwa index (script auto-builds if missing)
    bwa index -a is sex_marker.hg38.filtered.final.fasta

    # Run (R1 alone is sufficient; R2 only contributes if R1 has
    # fewer reads than --max-reads, e.g. low-coverage cfDNA)
    python sexcmd_xy.py \\
        --marker sex_marker.hg38.filtered.final.fasta \\
        --fastq R1.fq.gz R2.fq.gz \\
        --seq-type wgs \\
        --threads 8

Sequencing-type presets (--seq-type sets --max-reads default):
    wes ->   1,000,000 reads
    rna ->   5,000,000 reads   (stderr warning: tissue-dependent calls)
    wgs -> 150,000,000 reads

Output:
    <first-fastq>.OUTPUT      TSV byte-compatible with legacy SEXCMD.R
    stdout                    one-line summary (call / X reads / Y reads / ratio)
    stderr                    LOG lines: thread split, seq_type, mapped counts
"""
import argparse
import math
import os
import shlex
import shutil
import subprocess
import sys


# SEXCMD read-count targets: WES baseline * multiplier
SEQ_TYPE_MAX_READS = {
    "wes": 1_000_000,
    "rna": 5_000_000,
    "wgs": 150_000_000,
}

FEMALE_THRESHOLD = 0.2  # Y/X < threshold => female
TRIM_LEN = 101          # SEXCMD's hardcoded trim length

CFDNA_NOTE = (
    "cfDNA / liquid biopsy: works on low-pass WGS (0.1-5x); script simply "
    "EOFs early when input is smaller than --max-reads. For NIPT (maternal "
    "plasma with male fetus), expect an 'F' call because fetal Y signal is "
    "typically <20% of maternal X - this reflects the maternal genotype, "
    "not fetal sex."
)


# ----------------------------------------------------------------------------
# Utilities

def die(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def log(msg):
    print(f"LOG: {msg}", file=sys.stderr)


def require_cmd(cmd, hint=""):
    if shutil.which(cmd) is None:
        suffix = f" ({hint})" if hint else ""
        die(f"Required executable not found in PATH: {cmd}{suffix}")


def ensure_file(path):
    if not os.path.exists(path):
        die(f"File not found: {path}")


# ----------------------------------------------------------------------------
# FASTA + bwa index

def parse_fasta_in_order(fasta_path):
    """Return ordered [(name, length), ...] preserving FASTA order."""
    records = []
    name = None
    length = 0
    with open(fasta_path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        records.append((name, length))
    if not records:
        die(f"No FASTA records found in {fasta_path}")
    return records


def ensure_bwa_index(marker_fasta):
    required = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    missing = [marker_fasta + ext for ext in required
               if not os.path.exists(marker_fasta + ext)]
    if missing:
        log("building bwa index")
        rc = subprocess.run(["bwa", "index", "-a", "is", marker_fasta]).returncode
        if rc != 0:
            die(f"bwa index failed (exit {rc})")


# ----------------------------------------------------------------------------
# Thread split

def split_threads(total):
    """Split total thread budget across pigz / bwa / sambamba.

    Returns (decomp, bwa, sambamba). At low totals stages may share cores;
    at T>=4 the split is strict.
    """
    decomp = min(4, max(1, total // 4))
    sambamba = 1
    bwa = max(1, total - decomp - sambamba)
    return decomp, bwa, sambamba


# ----------------------------------------------------------------------------
# Pipeline

# Inline trim + N-filter in awk. Trim length passed via `-v t=N` at call site.
AWK_TRIM_FILTER = r'''
NR%4==1 { h=$0 }
NR%4==2 { s=substr($0,1,t); keep=(length(s)>0 && index(s,"N")==0) }
NR%4==3 { p=$0 }
NR%4==0 { if (keep) { q=substr($0,1,length(s)); print h ORS s ORS p ORS q } }
'''

# Count mapped reads per reference. Skips any SAM header lines defensively.
AWK_COUNT = r'''
/^@/ { next }
{ c[$3]++ }
END { for (k in c) print k "\t" c[k] }
'''


def build_decompress_cmd(fastqs, pigz_threads):
    """Build a bash fragment that streams all FASTQs decompressed."""
    parts = []
    for fq in fastqs:
        q = shlex.quote(fq)
        if fq.endswith(".gz"):
            parts.append(f"pigz -dc -p {pigz_threads} {q}")
        else:
            parts.append(f"cat {q}")
    if len(parts) == 1:
        return parts[0]
    return "( " + "; ".join(parts) + " )"


def run_pipeline(marker_fasta, fastqs, max_reads, trim_len, mapq,
                 decomp_t, bwa_t, sam_t):
    """Run the full decompress|head|trim|bwa|sambamba|count pipeline.

    bwa and sambamba stderr are silenced (2>/dev/null). Other stages keep
    inherited stderr so pigz/awk/bash errors stay visible.
    """
    decomp = build_decompress_cmd(fastqs, decomp_t)

    pipeline = (
        f"{decomp} "
        f"| head -n {max_reads * 4} "
        f"| awk -v t={trim_len} {shlex.quote(AWK_TRIM_FILTER)} "
        f"| bwa mem -t {bwa_t} {shlex.quote(marker_fasta)} - 2>/dev/null "
        f'| sambamba view -S -f sam '
        f'-F "not unmapped and mapping_quality >= {mapq}" '
        f"-t {sam_t} /dev/stdin 2>/dev/null "
        f"| awk {shlex.quote(AWK_COUNT)}"
    )

    # No pipefail: upstream SIGPIPE from head closing early is expected and
    # benign. Real upstream errors surface as empty counts -> caught by the
    # "<=3 mapped" check in classify_xy (same failure mode as the R script).
    result = subprocess.run(
        ["bash", "-c", pipeline],
        stdout=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        die(
            f"Pipeline failed (exit {result.returncode}). "
            f"bwa/sambamba stderr is silenced; to debug, remove `2>/dev/null` "
            f"from run_pipeline. If sambamba errored on SAM-on-stdin, try "
            f"`samtools view -F 4 -q {mapq}` in the source."
        )

    counts = {}
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        ref, n = parts
        try:
            counts[ref] = int(n)
        except ValueError:
            continue
    return counts


# ----------------------------------------------------------------------------
# Classification + R-compatible output

def split_xy_markers(marker_records):
    """Partition marker records into (x_markers, y_markers) by chr prefix.

    Dies if either group is empty. Called once before the pipeline so a
    malformed marker FASTA fails fast rather than after bwa has run.
    """
    x_markers = [(n, L) for n, L in marker_records if n.startswith("chrX")]
    y_markers = [(n, L) for n, L in marker_records if n.startswith("chrY")]
    if not x_markers or not y_markers:
        die("Marker FASTA must contain both chrX* and chrY* records")
    return x_markers, y_markers


def classify_xy(x_markers, y_markers, counts):
    x_total = sum(counts.get(n, 0) for n, _ in x_markers)
    y_total = sum(counts.get(n, 0) for n, _ in y_markers)

    if (x_total + y_total) <= 3:
        die("Total mapped marker read count is <= 3. "
            "Increase --max-reads or use higher-coverage input.")

    ratio = (y_total / x_total) if x_total > 0 else math.inf
    label = "F" if ratio < FEMALE_THRESHOLD else "M"
    return x_total, y_total, ratio, label


def write_r_compatible_output(path, x_markers, y_markers, counts, ratio, label,
                              fastq_display):
    """Write SEXCMD.R's exact output layout.

    Per R: chrX markers as rows; chrY counts paired by row index (the
    positional pairing is inherited from the R script's cbind-by-index and
    is only correct when the marker FASTA has paired X/Y markers in order).
    """
    with open(path, "wt") as out:
        out.write("\tLength\tXcount\tYcount\n")
        for i, (xname, xlen) in enumerate(x_markers):
            xc = counts.get(xname, 0)
            yc = counts.get(y_markers[i][0], 0) if i < len(y_markers) else 0
            out.write(f"{xname}\t{xlen}\t{xc}\t{yc}\n")
        out.write(
            f"Sex_Determination"
            f"\tRatio of X and Y counts is {ratio}"
            f"\tSex of this sample {fastq_display} is {label}"
            f"\t\n"
        )


# ----------------------------------------------------------------------------
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Human XY sex inference from FASTQ using SEXCMD marker "
                    "FASTA. Port of SEXCMD.R, XY+human only.",
        epilog=CFDNA_NOTE,
    )
    ap.add_argument("--marker", required=True,
                    help="SEXCMD marker FASTA (e.g. "
                         "sex_marker.hg38.filtered.final.fasta)")
    ap.add_argument("--fastq", nargs="+", required=True,
                    help="One or more FASTQ(.gz) files. R1 only is fine; "
                         "R1+R2 roughly doubles usable signal.")
    ap.add_argument("--seq-type", choices=["wes", "rna", "wgs"], required=True,
                    help="Sequencing type. Sets --max-reads default: "
                         "WES=1e6, RNA=5e6, WGS=1.5e8 (matches SEXCMD.R).")
    ap.add_argument("--max-reads", type=int, default=None,
                    help="Override --seq-type default. Max reads sampled "
                         "across all FASTQs combined.")
    ap.add_argument("--mapq", type=int, default=30,
                    help="Minimum MAPQ (default 30).")
    ap.add_argument("--threads", type=int, default=1,
                    help="Total thread budget. Split across pigz/bwa/sambamba "
                         "at runtime (see startup log).")
    ap.add_argument("--out", default=None,
                    help="Output path. Default: <first-fastq>.OUTPUT "
                         "(matches legacy R).")
    args = ap.parse_args()

    if args.threads < 1:
        die("--threads must be >= 1")

    # Tool checks
    require_cmd("bash")
    require_cmd("awk")
    require_cmd("head")
    require_cmd("bwa", "conda install -c bioconda bwa")
    require_cmd("sambamba", "conda install -c bioconda sambamba")
    require_cmd("pigz", "apt install pigz | conda install -c conda-forge pigz")

    # File checks
    ensure_file(args.marker)
    for fq in args.fastq:
        ensure_file(fq)

    # Seq-type resolution + RNA warning
    max_reads = (args.max_reads if args.max_reads is not None
                 else SEQ_TYPE_MAX_READS[args.seq_type])
    if args.seq_type == "rna":
        print(
            "WARNING: RNA-Seq sex inference via marker mapping is noisy. "
            "X-inactivation escape and Y-chromosome expression vary across "
            "tissues, and results should be interpreted with caution "
            "(especially for single-tissue or low-expression inputs).",
            file=sys.stderr,
        )

    # Thread split + run plan
    decomp_t, bwa_t, sam_t = split_threads(args.threads)
    log(f"threads total={args.threads} -> pigz={decomp_t} "
        f"bwa={bwa_t} sambamba={sam_t}")
    log(f"seq_type={args.seq_type} max_reads={max_reads} "
        f"trim_len={TRIM_LEN} mapq={args.mapq}")

    # Parse + split markers before the pipeline so a bad FASTA fails fast,
    # not after bwa has already run.
    ensure_bwa_index(args.marker)
    marker_records = parse_fasta_in_order(args.marker)
    x_markers, y_markers = split_xy_markers(marker_records)

    # legacy lovemun/SEXCMD R-compatible output file
    out_path = args.out if args.out else args.fastq[0] + ".OUTPUT"

    counts = run_pipeline(
        marker_fasta=args.marker,
        fastqs=args.fastq,
        max_reads=max_reads,
        trim_len=TRIM_LEN,
        mapq=args.mapq,
        decomp_t=decomp_t,
        bwa_t=bwa_t,
        sam_t=sam_t,
    )

    x_total, y_total, ratio, label = classify_xy(x_markers, y_markers, counts)
    log(f"chrX_mapped={x_total} chrY_mapped={y_total} y_over_x={ratio}")

    write_r_compatible_output(
        path=out_path,
        x_markers=x_markers,
        y_markers=y_markers,
        counts=counts,
        ratio=ratio,
        label=label,
        fastq_display=args.fastq[0],
    )
    log(f"SEXCMD result saved to {out_path}")

    # Stdout summary
    print(f"call={label}")
    print(f"chrX_reads={x_total}")
    print(f"chrY_reads={y_total}")
    print(f"y_over_x_ratio={ratio}")
    print(f"output={out_path}")


if __name__ == "__main__":
    main()
