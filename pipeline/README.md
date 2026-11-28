# pipeline/

The steps that turn FASTQ into the depth tables that `data/` is derived from. You only
need these if you are starting from the reads; the published result reproduces from
`data/` alone.

All scripts assume a working tree at `/mnt/c/crc_strain_pilot` and are written for
bash on WSL/Linux. Comments and log output are in Japanese.

| Script | Step |
|---|---|
| `cn_download.py`, `de_download.py` | Resolve ENA URLs and fetch FASTQ for a cohort |
| `asm_get3.sh` | Fetch one run, taking file sizes from the ENA filereport API rather than from an HTTP HEAD (see [../CORRECTIONS.md](../CORRECTIONS.md) item 6 and the header comment) |
| `asm_get5.sh` | Fetch one run and **verify it against the `fastq_md5` the ENA manifest declares**; on a mismatch delete and refetch from zero rather than resuming, and fall back to `submitted_ftp` if the generated FASTQ never verifies |
| `pks_map.sh` | Map to the island reference AM229678.1 and write `samtools depth -Q 20` tables |
| `at_map.sh`, `cn_map.sh`, `de_map.sh` | Per-cohort drivers for the above, plus species abundance |
| `entero_ref.py` | Build the combined Enterobacteriaceae marker index and record fragment lengths |
| `asm_genomes.sh` | Map one sample to six complete Enterobacteriaceae genomes in a **single combined index**, so reads compete and MAPQ 20 removes cross-species bleed |
| `asm_measure.sh` | Map the same read set to the island and the marker sets, for the depth ratios |
| `asm_report.py` | Print the per-species depth / breadth / ratio table used in RESULTS.md |
| `asm_assemble.sh` | Assemble one sample with megahit, find the contigs carrying the island, and align them to the six genomes - the step that answers "whose island is it" rather than "not *E. coli*'s" |
| `asm_junction.sh` | Remap the reads to the island-bearing contig and count the evidence spanning the chromosome/island join, which is what rules out an assembler chimera. Produces `../data/assembly_junction_depth.tsv` |
| `asm_kleborate.sh` | Bin the contigs by species and run Kleborate (3.2.4, `-p kpsc`) on the *K. pneumoniae* bin. Its header records two earlier attempts that failed and why |
| `asm_linkage.sh` | Map the reads back to the **whole** assembly (21,759 contigs, one competitive index, MAPQ 20) for per-contig read depth and for read pairs spanning two contigs. The header says why the index must not be the bin alone and why proper-pair filtering would delete the signal. Produces `../data/kleborate/contig_read_depth.tsv` and `../data/kleborate/contig_links.tsv` |
| `asm_kleborate_plus.sh` | Add the island-bearing contigs to that bin and run Kleborate again. **Read this one before quoting the typing**: the ICEKp10 call comes only from here, and the header says why the addition was forced by the reference set and what it therefore cannot establish. See [../CORRECTIONS.md](../CORRECTIONS.md) item 16 |

## Conventions worth knowing before you reuse these

- **The reference set decides what the species assignment *can* say.** Of the six
  genomes, only *E. coli* IHE3034 and *C. koseri* BAA-895 carry a *clb* island;
  *K. pneumoniae* HS11286 shares 0.11% of the island's 40-mers and has none of its 18
  genes. So a contig made of *clb* can never be assigned to *Klebsiella* by this
  index, however much *Klebsiella* is in the sample. `../data/clb_in_six_genomes.tsv`
  records this for all six.
- **MAPQ 20 everywhere.** Lowering it lets the ~1,290 bp duplication back in and
  inflates breadth; raising it costs real coverage.
- **Reads, not fragments.** Mates are passed independently with `-U`, so all counts
  are reads. Published read counts in the source papers are often fragment counts, a
  factor of two away.
- **ENA FASTQ paths for 7-digit accessions** are
  `vol1/fastq/ERR101/00<last digit>/<accession>/`, not the last three digits.
  Do not build these by hand; ask the filereport API, as `asm_get3.sh` does.
- **The ENA filereport API always prepends `run_accession`** as the first column, even
  when you did not request it. Parse by header name, not column number.
- **Keep concurrent connections to ENA low.** It is publicly funded shared
  infrastructure. The earlier scripts take four parts of one file at a time;
  `asm_get5.sh` uses a single stream.
- **Never resume onto a partial file you have not verified.** `curl -C -` appends from
  the current file size. If the interrupted attempt left a gap or a half-written block,
  resuming completes a file that is **exactly the declared size and still corrupt**, so
  a size check can never catch it. Three runs reached full size this way and failed
  only at the integrity check. On a failed check, delete and refetch from zero.
- **`gzip -t` is a real integrity check**, not a formality: it inflates the whole
  stream and compares the stored CRC32 and length, so mid-file corruption fails it.
  It is what caught all three bad files here. What it cannot catch is a perfectly
  valid archive of the *wrong* content, which is why `asm_get5.sh` also checks the
  manifest's md5 — the point is to check against the archive's declared value rather
  than against the code that just did the download.
