# sexcmd-xy

Human XY sex inference from FASTQ via marker mapping. A standalone Python reimplementation of the SEXCMD classification step (XY + human only), using a `pigz | awk | bwa | sambamba` pipeline.

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

## Attribution

The marker-mapping approach and the reference hg38 marker FASTA are from SEXCMD by Seongmun Jeong et al. ([lovemun/SEXCMD](https://github.com/lovemun/SEXCMD)). This repository is an independent Python reimplementation of the classification step (originally `SEXCMD.R`) for the XY-human use case. It does not redistribute any upstream code. The marker FASTA must be obtained separately from the upstream repository.

## Requirements

- Python 3.8+
- `bwa`, `sambamba`, `pigz` on `PATH`

Install via conda:

```bash
conda install -c bioconda bwa sambamba
conda install -c conda-forge pigz
```

## Installation

```bash
git clone https://github.com/odinokov/sexcmd-xy
cd sexcmd-xy
chmod +x sexcmd_xy.py
```

## Usage

```bash
# One-time: fetch the hg38 marker FASTA from upstream SEXCMD
wget -O sex_marker.hg38.filtered.final.fasta \
    https://raw.githubusercontent.com/lovemun/SEXCMD/master/Examples/Human/sex_marker.hg38.filtered.final.fasta

# Run (bwa index is auto-built on first use)
./sexcmd_xy.py \
    --marker sex_marker.hg38.filtered.final.fasta \
    --fastq sample_R1.fq.gz \
    --seq-type wgs \
    --threads 16

# Multi-lane: pass all lanes together, reads are pooled up to --max-reads
./sexcmd_xy.py \
    --marker sex_marker.hg38.filtered.final.fasta \
    --fastq lane1_R1.fq.gz lane2_R1.fq.gz lane3_R1.fq.gz lane4_R1.fq.gz \
    --seq-type wgs \
    --threads 16 \
    --out sample_id.OUTPUT
```

`--seq-type` sets the read-sampling budget: `wes`=1M, `rna`=5M, `wgs`=150M. Override with `--max-reads`. See `./sexcmd_xy.py --help` for all options.

## Output

A TSV at `<first-fastq>.OUTPUT`, byte-compatible with the legacy `SEXCMD.R` output format. The final row contains the sex call (`F` or `M`) and the Y/X read-count ratio. A one-line summary is also printed to stdout.

## Limitations

- Reads are mapped single-end. Passing both R1 and R2 only helps when R1 contains fewer reads than `--max-reads` (e.g., low-coverage cfDNA); otherwise, R2 is not consumed. For a typical WGS sample, pass R1 only.
- The marker FASTA's chromosome build must match the sample (e.g., hg38 markers for hg38-derived reads). Mismatched builds produce silently wrong counts.
- XY human only. The upstream SEXCMD supports ZW and non-human species; this port does not.
- Per-row output pairs chrX and chrY markers by FASTA position (inherited from `SEXCMD.R`). The totals and classification are unaffected, but per-row X/Y pairs are only meaningful when the marker FASTA has matched X/Y markers in corresponding positions, as the official SEXCMD hg38 file does.

## License

[MIT](LICENSE)
