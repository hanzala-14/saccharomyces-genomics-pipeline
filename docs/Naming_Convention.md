# Module 1 — Input Contracts & Naming Convention

## Purpose
This document defines naming, input rules, and traceability conventions for Module 1 (ENA + local FASTQ hybrid ingestion).

---

## 1. Samplesheet schema

Required header:

```csv
strain_id,accession,r1,r2
```

### Row modes
- **ENA mode:** `accession` set, `r1/r2` empty
- **Local mode:** `r1` set, `r2` optional, `accession` empty

### Invalid rows
- both accession and r1 missing
- accession and local paths provided together in same row

### Path resolution
Local paths are resolved at channel creation time using `file(path, checkIfExists: true)`:
- Absolute paths (e.g., `/home/user/data/sample_R1.fastq.gz`) used as-is
- Relative paths resolved from the Nextflow launch directory
- Missing or empty files cause immediate pipeline abort with a clear error message

---

## 2. Naming conventions

### `strain_id`
- required
- sanitized to safe characters: `[A-Za-z0-9_.-]`
- underscore replaces any invalid character

### `run_id`
- **ENA:** accession (e.g., `SRR10047173`)
- **Local:** `LOCAL_<strain_id>_<r1_basename>` (sanitized)
  - Example: file `WLP830.masodikkk.R1.fastq.gz` → `run_id = LOCAL_WLP830_WLP830.masodikkk.R1.fastq`
  - The basename (without `.gz`) provides uniqueness per file

### Output file naming (per strain)
- PE merged: `<strain_id>_PE_R1.fastq.gz`, `<strain_id>_PE_R2.fastq.gz`
- SE merged: `<strain_id>_SE.fastq.gz`

### Internal task file naming (INGEST_LOCAL work directory)
- `<run_id>_R1.fastq.gz`
- `<run_id>_R2.fastq.gz` (empty file if no R2 provided)

---

## 3. Traceability

Pipeline traces runs via the combination of:
- `strain_id` — identifies the biological sample
- `run_id` — identifies the specific sequencing run or local file

This combination is globally unique within a pipeline execution and appears in:
- `run_qc.tsv` — per-run classification audit
- `strain_qc.tsv` — per-strain aggregate summary
- Process tags — visible in Nextflow log output

---

## 4. Output naming

Per strain:
- `<strain_id>_PE_R1.fastq.gz`
- `<strain_id>_PE_R2.fastq.gz`
- `<strain_id>_SE.fastq.gz` (if applicable)

Global:
- `run_qc.tsv`
- `strain_qc.tsv`
- `strains.txt`

---

## 5. Example mixed samplesheet

```csv
strain_id,accession,r1,r2
WLP830,SRR10047173,,
yHCT81-Peris,SRR2586163,,
yHCT81-Peris,SRR2586165,,
WLP830,,/home/user1/data/raw/WLP830.masodikkk.R1.fastq.gz,/home/user1/data/raw/WLP830.masodikkk.R2.fastq.gz
WLP830,,/home/user1/data/raw/WLP830.harmadikkk.R1.fastq.gz,/home/user1/data/raw/WLP830.harmadikkk.R2.fastq.gz
CBS1538-Langdon,,/home/user1/data/raw/CBS1538-Langdon.elsokkk.R1.fastq.gz,/home/user1/data/raw/CBS1538-Langdon.elsokkk.R2.fastq.gz
WLP820,,/home/user1/data/raw/WLP820.otodikkk.R1.fastq.gz,/home/user1/data/raw/WLP820.otodikkk.R2.fastq.gz
WLP820,,/home/user1/data/raw/WLP820.elsokkk.R1.fastq.gz,/home/user1/data/raw/WLP820.elsokkk.R2.fastq.gz
WY2007,,/home/user1/data/raw/WY2007.elsokkk.R1.fastq.gz,/home/user1/data/raw/WY2007.elsokkk.R2.fastq.gz
```

---

## 6. Merge behavior

When a strain has multiple runs (from any source combination):
- All `PE_PASS` runs are concatenated into a single `_PE_R1.fastq.gz` / `_PE_R2.fastq.gz`
- All `SE_INPUT` / `SE_FALLBACK` runs are concatenated into a single `_SE.fastq.gz`
- Concatenation order is deterministic: sorted by `run_id` alphabetically
- A strain can have both PE and SE outputs (final_mode = `MIXED`)

### Example: WLP830 with 3 runs (1 ENA + 2 local)

| run_id | source | classification |
|--------|--------|----------------|
| SRR10047173 | ENA | PE_PASS |
| LOCAL_WLP830_WLP830.harmadikkk.R1.fastq | local | PE_PASS |
| LOCAL_WLP830_WLP830.masodikkk.R1.fastq | local | PE_PASS |

Result: `WLP830_PE_R1.fastq.gz` and `WLP830_PE_R2.fastq.gz` containing all 3 runs merged in run_id-sorted order.