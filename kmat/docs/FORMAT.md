# kmat on-disk matrix formats

This document freezes the Phase 1 (v1) and Phase 2 (v2) `.kmat` layouts. Little-endian throughout.

## Common header (40 bytes)

| Offset | Type | Field |
|---|---|---|
| 0 | `char[4]` | Magic `KMAT` |
| 4 | `uint32` | `version` (1 or 2) |
| 8 | `uint32` | `kmer_size` |
| 12 | `uint32` | `num_accessions` |
| 16 | `uint32` | `num_stripes` (= `ceil(num_accessions / 64)`) |
| 20 | `uint64` | `num_rows` — **v1:** number of k-mer rows; **v2:** number of k-mer map entries |
| 28 | `uint64` | `reserved` — **v1:** unused (0); **v2:** number of unique PA patterns |

Presence bits: accession `g` → stripe `g / 64`, bit `g % 64` (LSB-first within the word).

---

## Version 1 — dense k-mer rows (Phase 1)

After the header, `num_rows` records:

```text
{ uint64 kmer_code; uint64 words[num_stripes]; }
```

Rows are sorted by ascending `kmer_code`. Used by early Phase 1 builds; readers still accept v1 and inflate to the in-memory pattern representation.

---

## Version 2 — pattern-compressed (Phase 2)

Stores each unique presence bit-vector **once**, then maps every k-mer to a `pattern_id`.

### Layout after header

1. **Pattern store** — `reserved` (= `num_patterns`) entries, each:

   ```text
   uint64 words[num_stripes]
   ```

   Pattern IDs are `0 .. num_patterns-1` in file order.

2. **K-mer map** — `num_rows` (= `num_kmers`) entries, sorted by `kmer_code`:

   ```text
   { uint64 kmer_code; uint32 pattern_id; uint32 pad; }  // 16 bytes
   ```

   `pattern_id` must be `< num_patterns`. `pad` is zero.

### Why

On wheat-scale data, many k-mers share the same PA vector (~2/3 duplicates observed). v2:

- Shrinks on-disk size toward unique-pattern cardinality.
- Lets **GWAS** score each unique `z` once, then expand hits to member k-mers.
- Lets **gene** lookup do binary search on the k-mer map, then fetch one pattern.

### Build

`kmat build` writes **v2** by default: online dedup of presence vectors while unioning accession k-mers. Inputs may be FASTA/FASTQ (optional `.gz`) or pre-counted `.kset` presence sets (see below).

### Codec choice (locked for v2)

Plain packed `uint64` bit-vectors per pattern (no roaring). Revisit only if a later bench shows clear wins.

---

## Presence set (`.kset`) — Phase 3

Per-accession filtered k-mer presence set produced by `kmat count` (or `kmat import-kmers`). Little-endian.

| Offset | Type | Field |
|---|---|---|
| 0 | `char[4]` | Magic `KSET` |
| 4 | `uint32` | `version` (= 1) |
| 8 | `uint32` | `kmer_size` |
| 12 | `uint32` | `min_count` (filter used at count time) |
| 16 | `uint32` | `reserved` |
| 20 | `uint64` | `num_kmers` |
| 28 | `uint64[num_kmers]` | Sorted unique k-mer codes |

`kmat build` detects `.kset` paths in the accession list and unions them into a v2 matrix (same list order as GWAS columns). Do not mix `.kset` and sequence paths in one list.

Production `.kset` files are usually produced by `kmat count --engine kmc` (KMC CLI → dump → encode). `--engine builtin` is for small fixtures only.

---

## In-memory model

Readers load either version into:

- `patterns[pattern_id][stripe]`
- `kmers[]` sorted `{kmer_code, pattern_id}`

v1 files are converted on read (each distinct word-vector becomes a pattern).
