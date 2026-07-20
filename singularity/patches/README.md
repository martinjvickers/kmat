# KMC patches for the kmat Singularity image

Stock KMC 3.2.4 release binaries vendor Cloudflare zlib and inflate `.gz` inputs
with a hand-rolled loop in `fastq_reader.cpp`. On many real Illumina / pigz
`.fq.gz` files that fails with:

```text
Error: Some error while reading gzip file in (kmc_core/fastq_reader.cpp: 1062)
```

even though `zcat` / system zlib read the same files fine.

| File | Purpose |
|---|---|
| `binary_reader.h` | Decompress `.gz` with `gzopen`/`gzread`, feed plain FASTQ packs to readers |
| `apply-kmc-system-zlib.sh` | Patch KMC `Makefile` + `fastq_reader.h` includes to use system `<zlib.h>` / `-lz` |

Applied during `singularity/kmat.def` `%post`.
